from __future__ import annotations

import csv
import json
from datetime import datetime
from pathlib import Path
from typing import Any

from pcsec_pichia.core.paths import ProjectPaths


FACT_PACK_SCHEMA_VERSION = 1
TARGET_KEYS = ("hLF", "OPN")
SECRETION_UP_THRESHOLD = 1.01
SECRETION_DOWN_THRESHOLD = 0.99
GROWTH_RISK_THRESHOLD = 0.99
TOP_N_PER_BUCKET = 20


def build_screen_report_fact_pack(
    paths: ProjectPaths,
    *,
    run_names: tuple[str, ...] | None = None,
    csv_paths: tuple[Path | str, ...] | None = None,
) -> dict[str, Any]:
    """Build a strict JSON-compatible fact pack from existing local_runs screen outputs."""
    source_files = _discover_source_csvs(paths, run_names=run_names, csv_paths=csv_paths)
    source_runs: list[dict[str, Any]] = []
    evidence_items: list[dict[str, Any]] = []
    protected_results_detected = paths.results_dir.exists()

    for csv_path in source_files:
        source_run = _source_run_for_csv(paths, csv_path)
        source_runs.append(source_run)
        evidence_items.extend(_evidence_items_from_csv(paths, csv_path, source_run["run_name"]))

    evidence_items = sorted(
        evidence_items,
        key=lambda item: (
            TARGET_KEYS.index(str(item["target_key"])) if str(item["target_key"]) in TARGET_KEYS else len(TARGET_KEYS),
            str(item["intervention_type"]),
            _recommendation_rank(str(item["recommendation_tier"])),
            -_optional_float(item.get("secretion_ratio_vs_wildtype"), default=0.0),
            str(item.get("gene_id") or item.get("reaction_id") or ""),
        ),
    )
    _assign_evidence_ids(evidence_items)
    targets = {target_key: _target_payload(target_key, evidence_items) for target_key in TARGET_KEYS}
    warnings = []
    if protected_results_detected:
        warnings.append("Results/ exists and is treated as protected legacy MATLAB reference; fact pack did not read it.")
    return {
        "schema_version": FACT_PACK_SCHEMA_VERSION,
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "source_runs": source_runs,
        "targets": targets,
        "evidence_items": [_public_evidence_item(item) for item in evidence_items],
        "warnings": warnings,
    }


def summarize_fact_pack(fact_pack: dict[str, Any]) -> dict[str, Any]:
    return {
        "schema_version": fact_pack.get("schema_version"),
        "source_run_count": len(fact_pack.get("source_runs") or []),
        "evidence_item_count": len(fact_pack.get("evidence_items") or []),
        "targets": {
            target_key: {
                "candidate_count": target.get("candidate_count", 0),
                "useful_ko": len(target.get("useful_ko_candidates") or []),
                "useful_oe": len(target.get("useful_oe_candidates") or []),
                "manual_review": len(target.get("manual_review_candidates") or []),
                "growth_risk": len(target.get("growth_risk_candidates") or []),
            }
            for target_key, target in (fact_pack.get("targets") or {}).items()
        },
        "warnings": list(fact_pack.get("warnings") or []),
    }


def _discover_source_csvs(
    paths: ProjectPaths,
    *,
    run_names: tuple[str, ...] | None,
    csv_paths: tuple[Path | str, ...] | None,
) -> list[Path]:
    if csv_paths:
        return [Path(path).resolve() for path in csv_paths if Path(path).exists()]
    if run_names:
        candidates: list[Path] = []
        for run_name in run_names:
            run_dir = paths.local_runs_dir / run_name
            status_path = run_dir / "status.json"
            csv_path = _csv_path_from_status(status_path)
            candidates.append(csv_path if csv_path else run_dir / "gene_tradeoff_rows.csv")
        return [path.resolve() for path in candidates if path.exists()]
    return sorted(
        paths.local_runs_dir.glob("*/gene_tradeoff_rows.csv"),
        key=lambda path: path.stat().st_mtime,
        reverse=True,
    )


def _source_run_for_csv(paths: ProjectPaths, csv_path: Path) -> dict[str, Any]:
    run_dir = csv_path.parent
    status_path = run_dir / "status.json"
    status = _read_json(status_path)
    rel_csv = _relative_or_absolute(paths.repo_root, csv_path)
    return {
        "run_name": run_dir.name,
        "source_file": rel_csv,
        "status_file": _relative_or_absolute(paths.repo_root, status_path) if status_path.exists() else "",
        "status": status.get("status", "unknown") if isinstance(status, dict) else "unknown",
        "targets": tuple(status.get("targets") or ()) if isinstance(status, dict) else (),
        "scope": status.get("scope", "") if isinstance(status, dict) else "",
        "updated_at": status.get("updated_at", "") if isinstance(status, dict) else "",
    }


def _evidence_items_from_csv(paths: ProjectPaths, csv_path: Path, source_run: str) -> list[dict[str, Any]]:
    items: list[dict[str, Any]] = []
    with csv_path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        for row_index, row in enumerate(reader, start=1):
            target_key = _target_key(str(row.get("target_id") or ""))
            if target_key is None:
                continue
            items.append(_normalize_evidence_item(paths, csv_path, source_run, row_index, row, target_key))
    return items


def _normalize_evidence_item(
    paths: ProjectPaths,
    csv_path: Path,
    source_run: str,
    row_index: int,
    row: dict[str, str],
    target_key: str,
) -> dict[str, Any]:
    numeric_fields = _numeric_fields(row)
    intervention_type = str(row.get("intervention_type") or "").strip()
    candidate_kind = str(row.get("candidate_kind") or "gene").strip() or "gene"
    secretion_ratio = _optional_float(row.get("secretion_ratio_vs_wildtype"))
    growth_retention = _optional_float(row.get("growth_retention_ratio"))
    max_feasible_mu = _optional_float(row.get("max_feasible_mu"))
    support_status = str(row.get("support_status") or row.get("ko_support_status") or row.get("oe_support_status") or "")
    explicit_model_gpr = _optional_bool(row.get("model_gpr_executable"))
    explicit_oe_proxy = _optional_bool(row.get("oe_reaction_proxy"))
    recommendation_tier = str(row.get("recommendation_tier") or "").strip() or _derive_recommendation_tier(
        row,
        intervention_type=intervention_type,
        candidate_kind=candidate_kind,
        secretion_ratio=secretion_ratio,
        growth_retention=growth_retention,
        max_feasible_mu=max_feasible_mu,
        support_status=support_status,
    )
    return {
        "target_key": target_key,
        "target_id": str(row.get("target_id") or ""),
        "source_run": source_run,
        "source_file": _relative_or_absolute(paths.repo_root, csv_path),
        "row_index": row_index,
        "gene_id": str(row.get("gene_id") or ""),
        "reaction_id": str(row.get("reaction_id") or row.get("affected_reactions") or ""),
        "canonical_gene_id": str(row.get("canonical_gene_id") or row.get("gene_id") or ""),
        "gene_display_name": str(row.get("gene_display_name") or row.get("common_name") or row.get("gene_id") or ""),
        "standard_symbol": str(row.get("standard_symbol") or ""),
        "protein_name": str(row.get("protein_name") or ""),
        "intervention_type": intervention_type,
        "candidate_kind": candidate_kind,
        "effect_label": str(row.get("effect_label") or _effect_label(secretion_ratio, max_feasible_mu)),
        "recommendation_tier": recommendation_tier,
        "recommendation_tier_reason": str(row.get("recommendation_tier_reason") or _recommendation_reason(row, recommendation_tier)),
        "delta_objective": _optional_float(row.get("delta_objective")),
        "growth_retention_ratio": growth_retention,
        "secretion_ratio_vs_wildtype": secretion_ratio,
        "max_feasible_mu": max_feasible_mu,
        "support_status": support_status,
        "database_annotation_sources": _json_list(row.get("database_annotation_sources") or row.get("external_gene_function_sources")),
        "database_annotation_confidence": str(row.get("database_annotation_confidence") or ""),
        "model_operable": _model_operable(row, support_status),
        "model_gpr_executable": (
            explicit_model_gpr
            if explicit_model_gpr is not None
            else ("gpr" in support_status.lower() or support_status.endswith("gene_deletion"))
        ),
        "oe_reaction_proxy": (
            explicit_oe_proxy
            if explicit_oe_proxy is not None
            else (intervention_type == "OE" or "OE" in intervention_type)
        ),
        "phenotype_evidence": _json_value(row.get("phenotype_evidence") or row.get("phenotype_evidence_tier") or ""),
        "external_gene_function_evidence": _json_list(row.get("external_gene_function_evidence")),
        "external_gene_function_confidence": _json_list(row.get("external_gene_function_confidence")),
        "external_gene_function_sources": _json_list(row.get("external_gene_function_sources")),
        "external_gene_function_warnings": _json_list(row.get("external_gene_function_warnings")),
        "external_gpr_candidate_evidence": _json_list(
            row.get("external_gpr_candidate_evidence") or row.get("external_gpr_candidates")
        ),
        "external_model_sources": _json_list(row.get("external_model_sources")),
        "gpr_source_priority": _json_value(row.get("gpr_source_priority")),
        "external_gpr_candidate_count": _optional_int(row.get("external_gpr_candidate_count"), default=0),
        "best_external_gpr_source": str(row.get("best_external_gpr_source") or ""),
        "external_gpr_mapping_status": _json_value(row.get("external_gpr_mapping_status")),
        "external_gpr_conflict_warnings": _json_list(row.get("external_gpr_conflict_warnings")),
        "manual_review_reasons": _json_list(row.get("manual_review_reasons")),
        "ko_oe_external_gene_evidence": _json_object(row.get("ko_oe_external_gene_evidence")),
        "evidence_type": str(row.get("evidence_type") or row.get("recommendation_tier") or ""),
        "homology_review_status": str(row.get("homology_review_status") or row.get("standard_name_status") or ""),
        "rule_transfer_status": str(row.get("rule_transfer_status") or ""),
        "warnings": _row_warnings(row, recommendation_tier, intervention_type),
        "numeric_fields": numeric_fields,
    }


def _assign_evidence_ids(items: list[dict[str, Any]]) -> None:
    counters: dict[tuple[str, str], int] = {}
    for item in items:
        key = (str(item["target_key"]), str(item["intervention_type"]) or "UNK")
        counters[key] = counters.get(key, 0) + 1
        item["evidence_id"] = f"{key[0]}-{key[1] or 'UNK'}-{counters[key]:04d}"


def _target_payload(target_key: str, items: list[dict[str, Any]]) -> dict[str, Any]:
    target_items = [item for item in items if item["target_key"] == target_key]
    useful_ko = [item for item in target_items if item["intervention_type"] == "KO" and _is_useful_executable(item)]
    useful_oe = [item for item in target_items if item["intervention_type"] == "OE" and _is_useful_executable(item)]
    growth_risk = [item for item in target_items if _is_growth_risk(item)]
    manual_review = [item for item in target_items if _is_manual_review(item)]
    model_external = [item for item in target_items if _is_model_external_or_homology(item)]
    return {
        "candidate_count": len(target_items),
        "top_candidates": _brief_items(sorted(target_items, key=_item_sort_key)[:TOP_N_PER_BUCKET]),
        "useful_ko_candidates": _brief_items(sorted(useful_ko, key=_item_sort_key)[:TOP_N_PER_BUCKET]),
        "useful_oe_candidates": _brief_items(sorted(useful_oe, key=_item_sort_key)[:TOP_N_PER_BUCKET]),
        "growth_risk_candidates": _brief_items(sorted(growth_risk, key=_item_sort_key)[:TOP_N_PER_BUCKET]),
        "manual_review_candidates": _brief_items(sorted(manual_review, key=_item_sort_key)[:TOP_N_PER_BUCKET]),
        "model_external_candidates": _brief_items(sorted(model_external, key=_item_sort_key)[:TOP_N_PER_BUCKET]),
        "warnings": _target_warnings(target_key, target_items),
    }


def _brief_items(items: list[dict[str, Any]]) -> list[dict[str, Any]]:
    fields = (
        "evidence_id",
        "target_id",
        "gene_id",
        "reaction_id",
        "standard_symbol",
        "protein_name",
        "intervention_type",
        "effect_label",
        "recommendation_tier",
        "database_annotation_sources",
        "database_annotation_confidence",
        "model_gpr_executable",
        "oe_reaction_proxy",
        "phenotype_evidence",
        "external_gene_function_confidence",
        "external_gene_function_sources",
        "external_gene_function_evidence",
        "external_gpr_candidate_evidence",
        "external_model_sources",
        "gpr_source_priority",
        "external_gpr_candidate_count",
        "best_external_gpr_source",
        "external_gpr_mapping_status",
        "external_gpr_conflict_warnings",
        "manual_review_reasons",
        "ko_oe_external_gene_evidence",
        "secretion_ratio_vs_wildtype",
        "growth_retention_ratio",
        "max_feasible_mu",
        "warnings",
    )
    return [{field: item.get(field) for field in fields} for item in items]


def _public_evidence_item(item: dict[str, Any]) -> dict[str, Any]:
    return {key: value for key, value in item.items() if key != "target_key"}


def _target_warnings(target_key: str, items: list[dict[str, Any]]) -> list[str]:
    warnings: list[str] = []
    if not items:
        warnings.append(f"No existing local_runs screen rows were found for {target_key}.")
    if any(item.get("oe_reaction_proxy") for item in items):
        warnings.append("OE rows are reaction-level proxy evidence unless a row explicitly proves gene-level OE.")
    if any(_is_model_external_or_homology(item) for item in items):
        warnings.append("Homology/model-external rows require manual review and are not executable model recommendations.")
    return warnings


def _derive_recommendation_tier(
    row: dict[str, str],
    *,
    intervention_type: str,
    candidate_kind: str,
    secretion_ratio: float | None,
    growth_retention: float | None,
    max_feasible_mu: float | None,
    support_status: str,
) -> str:
    if _is_model_external_or_homology_row(row, candidate_kind):
        return "manual_review_homology_auxiliary"
    if max_feasible_mu is None and intervention_type == "KO":
        return "not_recommended_growth_risk"
    if secretion_ratio is None:
        return "manual_review_solver_inconclusive"
    if secretion_ratio <= SECRETION_UP_THRESHOLD:
        return "not_recommended_no_model_gain"
    if growth_retention is not None and growth_retention < GROWTH_RISK_THRESHOLD:
        return "model_executable_growth_risk"
    if intervention_type == "OE":
        return "promising_but_oe_proxy"
    if "runnable" in support_status or "gpr" in support_status.lower() or candidate_kind == "gene":
        return "strong_model_candidate"
    return "biology_interesting_manual_review"


def _recommendation_reason(row: dict[str, str], tier: str) -> str:
    if tier == "strong_model_candidate":
        return "Model-executable candidate with secretion increase and manageable growth cost in existing screen rows."
    if tier == "promising_but_oe_proxy":
        return "Secretion increases in an OE row, but OE is treated as reaction-level proxy evidence."
    if "growth_risk" in tier:
        return "Candidate has model gain but growth feasibility or retention needs review."
    if "homology" in tier or "manual" in tier:
        return "Evidence requires manual review before wet-lab actionability."
    return "Existing screen row does not support prioritization."


def _row_warnings(row: dict[str, str], recommendation_tier: str, intervention_type: str) -> list[str]:
    warnings: list[str] = []
    if intervention_type == "OE":
        warnings.append("OE is reaction-level proxy evidence unless explicitly stated otherwise.")
    if "homology" in recommendation_tier or "auxiliary" in recommendation_tier:
        warnings.append("Homology/annotation evidence is not phenotype validation.")
    if str(row.get("has_timeout") or "").lower() == "true":
        warnings.append("Solver timeout evidence present; interpret cautiously.")
    if str(row.get("skipped_reason") or "").strip():
        warnings.append(f"Skipped row: {row.get('skipped_reason')}")
    warnings.extend(str(item) for item in _json_list(row.get("external_gpr_conflict_warnings")) if str(item).strip())
    return warnings


def _effect_label(secretion_ratio: float | None, max_feasible_mu: float | None) -> str:
    if max_feasible_mu is None:
        return "求解失败或不可行"
    if secretion_ratio is None:
        return "无可比较分泌结果"
    if secretion_ratio > SECRETION_UP_THRESHOLD:
        return "提升分泌"
    if secretion_ratio < SECRETION_DOWN_THRESHOLD:
        return "降低分泌"
    return "无明显变化"


def _model_operable(row: dict[str, str], support_status: str) -> bool:
    text = " ".join(str(value or "") for value in (support_status, row.get("candidate_kind"), row.get("mapping_confidence")))
    lowered = text.lower()
    return "runnable" in lowered or "gpr" in lowered or "gene" in lowered or "proxy" in lowered


def _is_useful_executable(item: dict[str, Any]) -> bool:
    return (
        item.get("secretion_ratio_vs_wildtype") is not None
        and float(item["secretion_ratio_vs_wildtype"]) > SECRETION_UP_THRESHOLD
        and bool(item.get("model_operable"))
        and not _is_model_external_or_homology(item)
    )


def _is_growth_risk(item: dict[str, Any]) -> bool:
    growth = item.get("growth_retention_ratio")
    return item.get("max_feasible_mu") is None or (growth is not None and float(growth) < GROWTH_RISK_THRESHOLD)


def _is_manual_review(item: dict[str, Any]) -> bool:
    tier = str(item.get("recommendation_tier") or "").lower()
    return (
        "manual" in tier
        or "review" in tier
        or bool(item.get("manual_review_reasons"))
        or bool(item.get("external_gpr_conflict_warnings"))
        or _is_model_external_or_homology(item)
    )


def _is_model_external_or_homology(item: dict[str, Any]) -> bool:
    return _is_model_external_or_homology_row({key: str(value or "") for key, value in item.items()}, str(item.get("candidate_kind") or ""))


def _is_model_external_or_homology_row(row: dict[str, str], candidate_kind: str) -> bool:
    text = " ".join(
        str(value or "")
        for value in (
            candidate_kind,
            row.get("recommendation_tier"),
            row.get("evidence_type"),
            row.get("homology_review_status"),
            row.get("support_status"),
            row.get("standard_name_status"),
        )
    ).lower()
    return "model_external" in text or "homology" in text or "auxiliary" in text


def _item_sort_key(item: dict[str, Any]) -> tuple[int, int, float, float, str]:
    secretion = _optional_float(item.get("secretion_ratio_vs_wildtype"), default=0.0)
    growth = _optional_float(item.get("growth_retention_ratio"), default=0.0)
    return (
        _recommendation_rank(str(item.get("recommendation_tier") or "")),
        1 if item.get("model_operable") else 2,
        -secretion,
        -growth,
        str(item.get("evidence_id") or item.get("gene_id") or ""),
    )


def _recommendation_rank(tier: str) -> int:
    lowered = tier.lower()
    if "strong" in lowered:
        return 0
    if "promising" in lowered:
        return 1
    if "growth_risk" in lowered:
        return 2
    if "manual" in lowered or "review" in lowered:
        return 3
    if "not_recommended" in lowered:
        return 4
    return 5


def _numeric_fields(row: dict[str, str]) -> dict[str, float]:
    numeric: dict[str, float] = {}
    for key, value in row.items():
        parsed = _optional_float(value)
        if parsed is not None:
            numeric[key] = parsed
    return numeric


def _json_value(value: Any) -> Any:
    if value is None:
        return ""
    if isinstance(value, (dict, list, tuple)):
        return value
    text = str(value).strip()
    if not text:
        return ""
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        return text


def _json_list(value: Any) -> list[Any]:
    parsed = _json_value(value)
    if parsed in ("", None):
        return []
    if isinstance(parsed, list):
        return parsed
    if isinstance(parsed, tuple):
        return list(parsed)
    return [parsed]


def _json_object(value: Any) -> dict[str, Any]:
    parsed = _json_value(value)
    return parsed if isinstance(parsed, dict) else {}


def _optional_bool(value: Any) -> bool | None:
    if isinstance(value, bool):
        return value
    if value is None:
        return None
    text = str(value).strip().lower()
    if text in {"true", "1", "yes", "y"}:
        return True
    if text in {"false", "0", "no", "n"}:
        return False
    return None


def _optional_float(value: Any, default: float | None = None) -> float | None:
    if value is None:
        return default
    text = str(value).strip()
    if not text or text.lower() in {"nan", "none", "null"}:
        return default
    try:
        return float(text)
    except ValueError:
        return default


def _optional_int(value: Any, default: int = 0) -> int:
    parsed = _optional_float(value)
    return default if parsed is None else int(parsed)


def _target_key(target_id: str) -> str | None:
    upper = target_id.upper()
    if "HLF" in upper:
        return "hLF"
    if "OPN" in upper:
        return "OPN"
    return None


def _read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError, UnicodeDecodeError):
        return {}
    return data if isinstance(data, dict) else {}


def _csv_path_from_status(status_path: Path) -> Path | None:
    status = _read_json(status_path)
    csv_path = status.get("csv_path")
    return Path(csv_path) if csv_path else None


def _relative_or_absolute(root: Path, path: Path) -> str:
    try:
        return path.resolve().relative_to(root.resolve()).as_posix()
    except ValueError:
        return str(path)


__all__ = [
    "FACT_PACK_SCHEMA_VERSION",
    "build_screen_report_fact_pack",
    "summarize_fact_pack",
]
