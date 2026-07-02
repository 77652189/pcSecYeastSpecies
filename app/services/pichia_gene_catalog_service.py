from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from app import ensure_python_pichia_on_path

ensure_python_pichia_on_path()

from pcsec_pichia.core.paths import ProjectPaths

GENE_CATALOG_CACHE_SCHEMA_VERSION = 1
GENE_CATALOG_CACHE_DIR = Path("local_runs") / "streamlit_pichia_runs" / "gene_catalog_cache"
GENE_CATALOG_CACHE_FILE = "full_model_gene_catalog.json"
SECRETION_GENE_EVIDENCE_CACHE_FILE = "secretion_gene_evidence.json"


def list_pichia_gene_options(
    query: str = "",
    limit: int = 50,
    paths: ProjectPaths | None = None,
) -> list[dict[str, Any]]:
    resolved_paths = paths or ProjectPaths.discover(Path(__file__))
    ensure_python_pichia_on_path()

    from pcsec_pichia.loading import load_pcsec_pichia_inputs
    from pcsec_pichia.screens import default_ko_genes
    from pcsec_pichia.services.gene_catalog import search_full_catalog

    inputs = load_pcsec_pichia_inputs(resolved_paths.repo_root)
    model = inputs.prepared_model
    max_rows = max(1, min(int(limit), 200))
    query_text = str(query or "").strip().lower()
    preferred = set(default_ko_genes(model, 20)) if not query_text else set()

    rows = [
        _gene_option_row(row)
        for row in search_full_catalog(query_text, model=model)
        if query_text or row.get("gene_id") in preferred or row.get("primary_category") == "分泌相关"
    ]
    rows.sort(key=lambda row: _gene_lookup_sort_key(row, query_text, preferred))
    return rows[:max_rows]


def list_curated_pichia_gene_catalog(query: str = "") -> list[Any]:
    """Return curated secretion-pathway KO/OE entries through the app facade."""
    ensure_python_pichia_on_path()

    from pcsec_pichia.services.gene_catalog import search_catalog

    return search_catalog(query)


def list_curated_pichia_gene_catalog_by_category() -> dict[str, list[Any]]:
    """Return curated secretion-pathway KO/OE entries grouped by category."""
    ensure_python_pichia_on_path()

    from pcsec_pichia.services.gene_catalog import get_catalog_by_category

    return get_catalog_by_category()


def list_pichia_secretion_gene_evidence(
    query: str = "",
    *,
    force_refresh: bool = False,
    paths: ProjectPaths | None = None,
) -> list[dict[str, Any]]:
    """Return curated secretion gene names mapped to model GPR genes and reaction proxies."""
    ensure_python_pichia_on_path()

    cache_path = pichia_secretion_gene_evidence_cache_path(paths)
    rows = None if force_refresh else _read_gene_catalog_cache(cache_path)
    if rows is None:
        if force_refresh:
            from pcsec_pichia.services.gene_catalog import search_secretion_gene_evidence

            rows = search_secretion_gene_evidence("")
            _write_gene_catalog_cache(cache_path, rows)
        else:
            rows = _static_secretion_gene_evidence_rows()
    return _filter_secretion_gene_evidence_rows(rows, query)


def list_pichia_gene_rule_evidence(
    query: str = "",
    paths: ProjectPaths | None = None,
) -> list[dict[str, Any]]:
    """Return external-evidence GPR overlay rows for secretion-engineering names."""
    ensure_python_pichia_on_path()

    from pcsec_pichia.services.gene_rule_overlay import (
        DEFAULT_GENE_RULE_EVIDENCE_CACHE,
        load_gene_rule_evidence_cache,
    )

    resolved_paths = paths or ProjectPaths.discover(Path(__file__))
    query_text = str(query or "").strip().lower()
    rows = [
        record.to_dict()
        for record in load_gene_rule_evidence_cache(resolved_paths.repo_root / DEFAULT_GENE_RULE_EVIDENCE_CACHE).values()
    ]
    if not query_text:
        return rows
    return [
        row
        for row in rows
        if query_text in str(row.get("common_name") or "").lower()
        or query_text in str(row.get("candidate_locus_tag") or "").lower()
        or query_text in str(row.get("protein_name") or "").lower()
        or query_text in " ".join(str(item) for item in row.get("target_reaction_ids") or []).lower()
    ]


def list_verified_secretion_gene_library(
    query: str = "",
    *,
    paths: ProjectPaths | None = None,
) -> list[dict[str, Any]]:
    """Return a compact, evidence-aware secretion-engineering gene library for UI display.

    The rows are presentation records only. They do not add executable GPR rules and do not
    change KO/OE simulation semantics.
    """
    curated_rows = list_pichia_secretion_gene_evidence("", paths=paths)
    rule_rows = {
        str(row.get("common_name") or "").strip().lower(): row
        for row in list_pichia_gene_rule_evidence("", paths=paths)
    }
    rows = [
        _verified_secretion_gene_row(row, rule_rows.get(str(row.get("common_name") or "").strip().lower()))
        for row in curated_rows
    ]
    query_text = str(query or "").strip().lower()
    if not query_text:
        return rows
    return [row for row in rows if _verified_secretion_gene_row_matches(row, query_text)]


def load_pichia_full_model_gene_catalog(
    *,
    force_refresh: bool = False,
    paths: ProjectPaths | None = None,
) -> list[dict[str, object]]:
    """Return all model genes with lightweight reaction/process metadata."""
    ensure_python_pichia_on_path()

    cache_path = pichia_full_model_gene_catalog_cache_path(paths)
    if not force_refresh:
        cached_rows = _read_gene_catalog_cache(cache_path)
        if cached_rows is not None:
            return cached_rows

    from pcsec_pichia.services.gene_catalog import load_full_model_genes

    rows = load_full_model_genes()
    _write_gene_catalog_cache(cache_path, rows)
    return rows


def build_pichia_gene_evidence_cache(
    paths: ProjectPaths | None = None,
    progress=None,
    refresh_full_catalog: bool = True,
) -> dict[str, Any]:
    """Build wet-lab gene annotation evidence for all current model genes."""
    resolved_paths = paths or ProjectPaths.discover(Path(__file__))
    ensure_python_pichia_on_path()

    from pcsec_pichia.loading import load_pcsec_pichia_inputs
    from pcsec_pichia.services.gene_evidence import (
        DEFAULT_GENE_EVIDENCE_CACHE,
        DEFAULT_GENE_EVIDENCE_SUMMARY,
        build_gene_evidence_cache,
        load_gene_evidence_cache,
    )

    if progress:
        progress("Loading pcSecPichia model genes...")
    inputs = load_pcsec_pichia_inputs(resolved_paths.repo_root)
    gene_ids = [str(gene_id) for gene_id in inputs.prepared_model.genes]
    summary = build_gene_evidence_cache(
        gene_ids,
        output_path=resolved_paths.repo_root / DEFAULT_GENE_EVIDENCE_CACHE,
        summary_path=resolved_paths.repo_root / DEFAULT_GENE_EVIDENCE_SUMMARY,
        progress=progress,
    )
    if refresh_full_catalog:
        if progress:
            progress("Merging wet-lab evidence into full model gene catalog cache...")
        _merge_gene_evidence_into_catalog_cache(
            pichia_full_model_gene_catalog_cache_path(resolved_paths),
            load_gene_evidence_cache(resolved_paths.repo_root / DEFAULT_GENE_EVIDENCE_CACHE),
        )
    return summary


def pichia_full_model_gene_catalog_cache_path(paths: ProjectPaths | None = None) -> Path:
    resolved_paths = paths or ProjectPaths.discover(Path(__file__))
    return resolved_paths.repo_root / GENE_CATALOG_CACHE_DIR / GENE_CATALOG_CACHE_FILE


def pichia_secretion_gene_evidence_cache_path(paths: ProjectPaths | None = None) -> Path:
    resolved_paths = paths or ProjectPaths.discover(Path(__file__))
    return resolved_paths.repo_root / GENE_CATALOG_CACHE_DIR / SECRETION_GENE_EVIDENCE_CACHE_FILE


def _verified_secretion_gene_row(
    curated_row: dict[str, object],
    rule_evidence: dict[str, object] | None = None,
) -> dict[str, Any]:
    common_name = str(curated_row.get("common_name") or "")
    model_gene_id = str(curated_row.get("mapped_model_gene_id") or curated_row.get("declared_model_gene_id") or "")
    candidate_locus = str((rule_evidence or {}).get("candidate_locus_tag") or "")
    ko_reaction = str(curated_row.get("ko_reaction_id") or "")
    oe_reaction = str(curated_row.get("oe_reaction_id") or "")
    gene_level_ready = bool(curated_row.get("gene_level_ready") or model_gene_id)
    reaction_proxy_ready = bool(curated_row.get("reaction_proxy_ready") or ko_reaction or oe_reaction)
    operation_status = _verified_operation_status(
        gene_level_ready=gene_level_ready,
        ko_reaction=ko_reaction,
        oe_reaction=oe_reaction,
        reaction_proxy_ready=reaction_proxy_ready,
    )
    evidence_tier = _verified_evidence_tier(curated_row, rule_evidence, gene_level_ready, reaction_proxy_ready)
    recommended_use = _verified_recommended_use(curated_row, gene_level_ready, reaction_proxy_ready)
    detail_payload = {
        "curated": dict(curated_row),
        "rule_evidence": dict(rule_evidence or {}),
    }
    return {
        "display_name": common_name,
        "model_gene_id": model_gene_id,
        "locus_tag": candidate_locus or model_gene_id,
        "function_annotation": curated_row.get("description") or "",
        "operation_status": operation_status,
        "evidence_tier": evidence_tier,
        "recommended_use": recommended_use,
        "source_summary": _verified_source_summary(curated_row, rule_evidence),
        "category": curated_row.get("category") or "",
        "ko_reaction_id": ko_reaction,
        "oe_reaction_id": oe_reaction,
        "mapping_status": curated_row.get("mapping_status") or "",
        "rule_status": (rule_evidence or {}).get("rule_status") or "",
        "rule_confidence": (rule_evidence or {}).get("confidence") or "",
        "detail_payload": detail_payload,
    }


def _verified_operation_status(
    *,
    gene_level_ready: bool,
    ko_reaction: str,
    oe_reaction: str,
    reaction_proxy_ready: bool,
) -> str:
    actions: list[str] = []
    if gene_level_ready:
        actions.append("基因级 KO")
    proxy_actions: list[str] = []
    if ko_reaction:
        proxy_actions.append("KO")
    if oe_reaction:
        proxy_actions.append("OE")
    if proxy_actions:
        actions.append(f"反应级 {'/'.join(proxy_actions)} proxy")
    elif reaction_proxy_ready:
        actions.append("反应级 proxy")
    if not actions:
        return "仅证据展示 / 需人工确认"
    return "；".join(actions)


def _verified_evidence_tier(
    curated_row: dict[str, object],
    rule_evidence: dict[str, object] | None,
    gene_level_ready: bool,
    reaction_proxy_ready: bool,
) -> str:
    confidence = str((rule_evidence or {}).get("confidence") or "")
    rule_status = str((rule_evidence or {}).get("rule_status") or "")
    if gene_level_ready:
        return "模型可执行 GPR + curated 证据"
    if confidence == "high_exact_multi_source" and rule_status != "overlay_executable":
        return "数据库高置信 locus 候选；GPR 未执行"
    if reaction_proxy_ready:
        return "模型反应代理 + curated 证据"
    if confidence:
        return f"外部证据待复核：{confidence}"
    if curated_row.get("curated_evidence"):
        return "curated 文献/路径证据；需人工确认"
    return "仅模型/名称证据；需人工确认"


def _verified_recommended_use(
    curated_row: dict[str, object],
    gene_level_ready: bool,
    reaction_proxy_ready: bool,
) -> str:
    if gene_level_ready:
        return "可直接作为小批量 gene-level KO 输入；实验前仍需复核 locus 注释"
    if reaction_proxy_ready:
        return "可用于 reaction-level proxy 解释/筛选；湿实验前需确认 K. phaffii locus"
    label = str(curated_row.get("recommended_use") or "")
    if label == "manual_review_required":
        return "仅作候选名展示；需人工确认后再进入 KO/OE"
    return "需人工确认后使用"


def _verified_source_summary(
    curated_row: dict[str, object],
    rule_evidence: dict[str, object] | None,
) -> str:
    parts: list[str] = []
    if curated_row.get("curated_evidence"):
        parts.append(str(curated_row.get("curated_evidence")))
    if rule_evidence:
        confidence = str(rule_evidence.get("confidence") or "")
        sources = ", ".join(str(item) for item in rule_evidence.get("evidence_sources") or [])
        locus = str(rule_evidence.get("candidate_locus_tag") or "")
        rule_status = str(rule_evidence.get("rule_status") or "")
        evidence_bits = [bit for bit in (confidence, sources, locus, rule_status) if bit]
        if evidence_bits:
            parts.append("外部证据：" + " | ".join(evidence_bits))
    if curated_row.get("mapped_model_gene_id") or curated_row.get("declared_model_gene_id"):
        parts.append("已映射模型 GPR gene ID")
    elif curated_row.get("reaction_proxy_ready"):
        parts.append("模型中存在反应级代理")
    return "；".join(parts)


def _verified_secretion_gene_row_matches(row: dict[str, Any], query_text: str) -> bool:
    detail = row.get("detail_payload") if isinstance(row.get("detail_payload"), dict) else {}
    curated = detail.get("curated") if isinstance(detail.get("curated"), dict) else {}
    rule = detail.get("rule_evidence") if isinstance(detail.get("rule_evidence"), dict) else {}
    values: list[object] = [
        row.get("display_name"),
        row.get("model_gene_id"),
        row.get("locus_tag"),
        row.get("function_annotation"),
        row.get("operation_status"),
        row.get("evidence_tier"),
        row.get("recommended_use"),
        row.get("source_summary"),
        row.get("category"),
        row.get("ko_reaction_id"),
        row.get("oe_reaction_id"),
        row.get("mapping_status"),
        row.get("rule_status"),
        row.get("rule_confidence"),
        curated.get("homolog_note"),
        curated.get("mapped_display_name"),
        " ".join(str(item) for item in curated.get("mapped_aliases") or []),
        rule.get("protein_name"),
        " ".join(str(item) for item in rule.get("target_reaction_ids") or []),
    ]
    return any(query_text in str(value or "").lower() for value in values)


def get_pichia_ko_genes_for_selection(selected_names: list[str]) -> list[str]:
    ensure_python_pichia_on_path()

    from pcsec_pichia.services.gene_catalog import get_ko_genes_for_selection

    return get_ko_genes_for_selection(selected_names)


def get_pichia_ko_reactions_for_selection(selected_names: list[str]) -> list[str]:
    ensure_python_pichia_on_path()

    from pcsec_pichia.services.gene_catalog import get_ko_reactions_for_selection

    return get_ko_reactions_for_selection(selected_names)


def get_pichia_oe_reactions_for_selection(selected_names: list[str]) -> list[str]:
    ensure_python_pichia_on_path()

    from pcsec_pichia.services.gene_catalog import get_oe_reactions_for_selection

    return get_oe_reactions_for_selection(selected_names)


def _gene_option_row(row: dict[str, object]) -> dict[str, Any]:
    sample_reactions = [str(item) for item in row.get("sample_reactions") or []]
    processes = [item.strip() for item in str(row.get("processes") or "").split(",") if item.strip()]
    return {
        "gene_id": str(row.get("gene_id") or ""),
        "canonical_gene_id": str(row.get("canonical_gene_id") or row.get("gene_id") or ""),
        "aliases": [str(item) for item in row.get("aliases") or []],
        "reaction_count": int(row.get("n_reactions") or 0),
        "reactions_preview": sample_reactions[:5],
        "secretory_processes": processes,
        "gr_rule_preview": "",
        "suggested_use": "KO / OE gene proxy",
        "primary_category": row.get("primary_category"),
        "ko_support_status": row.get("ko_support_status") or "",
        "oe_support_status": row.get("oe_support_status") or "",
        "gpr_role": row.get("gpr_role") or "",
        "support_reason": row.get("support_reason") or "",
        "missing_information": list(row.get("missing_information") or []),
        "confidence": row.get("confidence") or "",
        "standard_gene_symbol": row.get("standard_gene_symbol") or "",
        "display_name": row.get("display_name") or "",
        "protein_name": row.get("protein_name") or "",
        "function_annotation": row.get("function_annotation") or "",
        "external_ids": dict(row.get("external_ids") or {}),
        "ec_numbers": list(row.get("ec_numbers") or []),
        "go_terms": list(row.get("go_terms") or []),
        "ortholog_symbol": row.get("ortholog_symbol") or "",
        "wet_lab_readiness": row.get("wet_lab_readiness") or "model_only_not_experiment_ready",
        "evidence_sources": list(row.get("evidence_sources") or []),
        "evidence_confidence": row.get("evidence_confidence") or "",
    }


def _read_gene_catalog_cache(cache_path: Path) -> list[dict[str, object]] | None:
    try:
        payload = json.loads(cache_path.read_text(encoding="utf-8"))
    except (FileNotFoundError, OSError, json.JSONDecodeError, UnicodeDecodeError):
        return None
    if not isinstance(payload, dict):
        return None
    if payload.get("schema_version") != GENE_CATALOG_CACHE_SCHEMA_VERSION:
        return None
    rows = payload.get("rows")
    if not isinstance(rows, list) or not all(isinstance(row, dict) for row in rows):
        return None
    return rows


def _write_gene_catalog_cache(cache_path: Path, rows: list[dict[str, object]]) -> None:
    payload = {
        "schema_version": GENE_CATALOG_CACHE_SCHEMA_VERSION,
        "rows": rows,
    }
    try:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        cache_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    except OSError:
        return


def _filter_secretion_gene_evidence_rows(
    rows: list[dict[str, object]],
    query: str = "",
) -> list[dict[str, Any]]:
    query_text = str(query or "").strip().lower()
    normalized_rows = [dict(row) for row in rows]
    if not query_text:
        return normalized_rows
    return [
        row
        for row in normalized_rows
        if _secretion_gene_evidence_row_matches(row, query_text)
    ]


def _secretion_gene_evidence_row_matches(row: dict[str, object], query_text: str) -> bool:
    reaction_text = " ".join(
        str(item.get("reaction_id") or "")
        for item in row.get("reaction_evidence") or []
        if isinstance(item, dict)
    )
    fields = (
        row.get("common_name"),
        row.get("category"),
        row.get("description"),
        row.get("curated_evidence"),
        row.get("homolog_note"),
        row.get("declared_model_gene_id"),
        row.get("mapped_model_gene_id"),
        row.get("mapped_display_name"),
        " ".join(str(alias) for alias in row.get("mapped_aliases") or []),
        row.get("mapping_status"),
        row.get("recommended_use"),
        row.get("ko_reaction_id"),
        row.get("oe_reaction_id"),
        reaction_text,
    )
    return any(query_text in str(value or "").lower() for value in fields)


def _static_secretion_gene_evidence_rows() -> list[dict[str, object]]:
    """Return a fast curated gene/proxy table without loading the pcSec model."""
    from pcsec_pichia.services.gene_catalog import SECRETION_GENE_CATALOG

    rows: list[dict[str, object]] = []
    for entry in SECRETION_GENE_CATALOG:
        proxy_reactions = [value for value in (entry.ko_reaction_id, entry.oe_reaction_id) if value]
        if entry.gene_id:
            mapping_status = "model_gpr_gene_available"
            recommended_use = "gene_level_gpr_perturbation"
        elif proxy_reactions:
            mapping_status = "reaction_proxy_only"
            recommended_use = "reaction_level_proxy_requires_locus_review"
        else:
            mapping_status = "literature_name_only"
            recommended_use = "manual_review_required"
        rows.append(
            {
                "common_name": entry.common_name,
                "category": entry.category,
                "description": entry.description,
                "curated_evidence": entry.evidence,
                "homolog_note": entry.homolog_note,
                "declared_model_gene_id": entry.gene_id,
                "mapped_model_gene_id": entry.gene_id,
                "mapped_display_name": "",
                "mapped_aliases": [],
                "mapping_status": mapping_status,
                "recommended_use": recommended_use,
                "ko_reaction_id": entry.ko_reaction_id,
                "oe_reaction_id": entry.oe_reaction_id,
                "reaction_evidence": [
                    {
                        "reaction_id": reaction_id,
                        "exists_in_model": True,
                        "reaction_index_1based": None,
                        "has_gpr_rule": False,
                        "rule": "",
                        "gr_rule": "",
                        "source": "static_curated_catalog",
                    }
                    for reaction_id in dict.fromkeys(proxy_reactions)
                ],
                "proxy_exists_in_model": bool(proxy_reactions),
                "proxy_has_gpr_rule": False,
                "gene_level_ready": bool(entry.gene_id),
                "reaction_proxy_ready": bool(proxy_reactions),
                "evidence_source": "static_curated_catalog",
            }
        )
    return rows


def _merge_gene_evidence_into_catalog_cache(cache_path: Path, evidence_by_gene: dict[str, Any]) -> None:
    cached_rows = _read_gene_catalog_cache(cache_path)
    if cached_rows is None:
        return
    merged_rows: list[dict[str, object]] = []
    for row in cached_rows:
        gene_id = str(row.get("gene_id") or row.get("canonical_gene_id") or "")
        evidence = evidence_by_gene.get(gene_id)
        if evidence is None:
            merged_rows.append({**row, "wet_lab_readiness": "model_only_not_experiment_ready"})
            continue
        merged_rows.append(
            {
                **row,
                "standard_gene_symbol": evidence.standard_gene_symbol,
                "display_name": evidence.display_name,
                "aliases": list(evidence.aliases),
                "external_ids": dict(evidence.external_ids or {}),
                "protein_name": evidence.protein_name,
                "function_annotation": evidence.function_annotation,
                "subcellular_location": evidence.subcellular_location,
                "ec_numbers": list(evidence.ec_numbers),
                "go_terms": list(evidence.go_terms),
                "ortholog_symbol": evidence.ortholog_symbol,
                "wet_lab_readiness": evidence.wet_lab_readiness,
                "evidence_sources": list(evidence.evidence_sources),
                "evidence_confidence": evidence.evidence_confidence,
                "last_refreshed": evidence.last_refreshed,
            }
        )
    _write_gene_catalog_cache(cache_path, merged_rows)


def _gene_lookup_sort_key(row: dict[str, Any], query: str, preferred: set[str]) -> tuple[int, int, str]:
    gene_id = str(row.get("gene_id") or "")
    processes = " ".join(str(item) for item in row.get("secretory_processes") or [])
    exact_rank = 0 if query and gene_id.lower().startswith(query) else 1
    preferred_rank = 0 if gene_id in preferred or _is_secretory_lookup_process(processes) else 1
    return (exact_rank, preferred_rank, gene_id)


def _is_secretory_lookup_process(process: str) -> bool:
    return any(token in process for token in ("翻译", "ER", "DSB", "糖基化", "错误折叠", "Golgi", "分泌"))


__all__ = [
    "GENE_CATALOG_CACHE_DIR",
    "SECRETION_GENE_EVIDENCE_CACHE_FILE",
    "build_pichia_gene_evidence_cache",
    "get_pichia_ko_genes_for_selection",
    "get_pichia_oe_reactions_for_selection",
    "get_pichia_ko_reactions_for_selection",
    "list_curated_pichia_gene_catalog",
    "list_curated_pichia_gene_catalog_by_category",
    "list_pichia_secretion_gene_evidence",
    "list_pichia_gene_rule_evidence",
    "list_pichia_gene_options",
    "list_verified_secretion_gene_library",
    "load_pichia_full_model_gene_catalog",
    "pichia_full_model_gene_catalog_cache_path",
    "pichia_secretion_gene_evidence_cache_path",
]
