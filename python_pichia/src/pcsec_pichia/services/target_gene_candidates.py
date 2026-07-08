from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

from pcsec_pichia.services.homology_evidence import GeneHomologyEvidence, homology_evidence_for_gene


HLF_OPN_CANDIDATE_SCHEMA_VERSION = 1

TARGET_HLF = "hLF"
TARGET_OPN = "OPN"
TARGET_SHARED = "shared"

MODEL_KO_EXECUTABLE = "model_ko_executable"
MODEL_OE_PROXY_EXECUTABLE = "model_oe_proxy_executable"
MODEL_EXPLAIN_ONLY_COMPLEX_SUBUNIT = "model_explain_only_complex_subunit"
NOT_IN_MODEL = "not_in_model"
UNRESOLVED_NAME = "unresolved_name"
MANUAL_REVIEW_REQUIRED = "manual_review_required"

MODEL_OPERABLE = "model_operable"
NOT_MODEL_OPERABLE = "not_model_operable"


@dataclass(frozen=True)
class HlfOpnCandidateGene:
    target_context: str
    gene_id: str
    candidate_role: str
    evidence_type: str
    evidence_confidence: str
    model_operable: bool
    recommended_intervention: str
    reason: str
    warnings: tuple[str, ...] = ()
    operability_status: str = MANUAL_REVIEW_REQUIRED
    model_operability_label: str = NOT_MODEL_OPERABLE
    source_common_name: str = ""
    source_category: str = ""
    display_name: str = ""
    standard_symbol: str = ""
    protein_name: str = ""
    external_ids: dict[str, str] = field(default_factory=dict)
    annotation_sources: tuple[str, ...] = ()
    source_reaction_ids: tuple[str, ...] = ()
    executable_ko_reactions: tuple[str, ...] = ()
    executable_oe_proxy_reactions: tuple[str, ...] = ()
    review_reactions: tuple[str, ...] = ()
    homology_review_status: str = ""
    rule_transfer_status: str = ""
    homology_query_symbol: str = ""
    source_gene_resolution: str = ""

    def to_dict(self) -> dict[str, Any]:
        payload = asdict(self)
        for key in (
            "warnings",
            "annotation_sources",
            "source_reaction_ids",
            "executable_ko_reactions",
            "executable_oe_proxy_reactions",
            "review_reactions",
        ):
            payload[key] = list(payload[key])
        payload["external_ids"] = dict(self.external_ids)
        return payload


def build_hlf_opn_candidate_gene_rows(
    *,
    full_model_rows: list[dict[str, object]] | tuple[dict[str, object], ...],
    standard_name_rows: list[dict[str, object]] | tuple[dict[str, object], ...],
    secretion_gene_catalog: tuple[Any, ...] | list[Any],
    homology_evidence_by_gene: dict[str, GeneHomologyEvidence] | None = None,
) -> tuple[HlfOpnCandidateGene, ...]:
    """Build target-context candidates without expanding beyond the curated secretion catalog."""

    full_by_gene = {_text(row.get("gene_id")): row for row in full_model_rows if _text(row.get("gene_id"))}
    standard_by_gene = {_text(row.get("gene_id")): row for row in standard_name_rows if _text(row.get("gene_id"))}
    rows: list[HlfOpnCandidateGene] = []
    seen: set[tuple[str, str, str, str]] = set()
    for entry in secretion_gene_catalog:
        contexts = _target_contexts_for_category(_entry_value(entry, "category"))
        homology = _homology_for_entry(entry, homology_evidence_by_gene)
        gene_id, resolution = _resolve_candidate_gene_id(entry, homology)
        for target_context in contexts:
            candidate = _candidate_row(
                entry=entry,
                target_context=target_context,
                gene_id=gene_id,
                source_gene_resolution=resolution,
                full_model_row=full_by_gene.get(gene_id),
                standard_name_row=standard_by_gene.get(gene_id),
                homology=homology,
            )
            key = (
                candidate.target_context,
                candidate.source_common_name,
                candidate.gene_id,
                "|".join(candidate.source_reaction_ids),
            )
            if key in seen:
                continue
            seen.add(key)
            rows.append(candidate)
    return tuple(sorted(rows, key=_candidate_sort_key))


def filter_hlf_opn_candidate_gene_rows(
    rows: tuple[HlfOpnCandidateGene, ...] | list[HlfOpnCandidateGene],
    *,
    target_context: str | None = None,
    include_shared: bool = True,
) -> tuple[HlfOpnCandidateGene, ...]:
    context = _normalize_target_context(target_context)
    if context is None:
        return tuple(rows)
    allowed = {context}
    if include_shared and context in {TARGET_HLF, TARGET_OPN}:
        allowed.add(TARGET_SHARED)
    return tuple(row for row in rows if row.target_context in allowed)


def executable_inputs_for_hlf_opn_candidates(
    rows: tuple[HlfOpnCandidateGene, ...] | list[HlfOpnCandidateGene],
    *,
    target_context: str,
    include_shared: bool = True,
) -> dict[str, Any]:
    filtered = filter_hlf_opn_candidate_gene_rows(
        rows,
        target_context=target_context,
        include_shared=include_shared,
    )
    ko_gene_ids = _dedupe(
        row.gene_id
        for row in filtered
        if row.recommended_intervention == "KO" and row.operability_status == MODEL_KO_EXECUTABLE
    )
    oe_gene_ids = _dedupe(
        row.gene_id
        for row in filtered
        if row.recommended_intervention == "OE" and row.operability_status == MODEL_OE_PROXY_EXECUTABLE
    )
    excluded = [
        row
        for row in filtered
        if row.operability_status not in {MODEL_KO_EXECUTABLE, MODEL_OE_PROXY_EXECUTABLE}
        or row.recommended_intervention == "review_only"
    ]
    return {
        "target_context": _normalize_target_context(target_context) or target_context,
        "ko_gene_ids": ko_gene_ids,
        "oe_gene_ids": oe_gene_ids,
        "excluded_count": len(excluded),
        "warnings": [
            "OE candidates are reaction-level proxies in the current model, not gene-level overexpression simulations.",
            "not_in_model and unresolved_name candidates are excluded from executable KO/OE inputs.",
        ],
    }


def summarize_hlf_opn_candidate_gene_rows(
    rows: tuple[HlfOpnCandidateGene, ...] | list[HlfOpnCandidateGene],
) -> dict[str, Any]:
    all_rows = tuple(rows)
    return {
        "schema_version": HLF_OPN_CANDIDATE_SCHEMA_VERSION,
        "total_candidates": len(all_rows),
        "target_context_counts": _count_by(all_rows, "target_context"),
        "target_candidate_counts": {
            TARGET_HLF: len(filter_hlf_opn_candidate_gene_rows(all_rows, target_context=TARGET_HLF)),
            TARGET_OPN: len(filter_hlf_opn_candidate_gene_rows(all_rows, target_context=TARGET_OPN)),
        },
        "operability_status_counts": _count_by(all_rows, "operability_status"),
        "recommended_intervention_counts": _count_by(all_rows, "recommended_intervention"),
        "evidence_type_counts": _count_by(all_rows, "evidence_type"),
        "model_operable_count": sum(1 for row in all_rows if row.model_operable),
        "not_model_operable_count": sum(1 for row in all_rows if not row.model_operable),
    }


def write_hlf_opn_candidate_gene_cache(
    rows: tuple[HlfOpnCandidateGene, ...] | list[HlfOpnCandidateGene],
    output_path: Path,
) -> None:
    payload = {
        "schema_version": HLF_OPN_CANDIDATE_SCHEMA_VERSION,
        "summary": summarize_hlf_opn_candidate_gene_rows(tuple(rows)),
        "rows": [row.to_dict() for row in rows],
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")


def load_hlf_opn_candidate_gene_cache(path: Path) -> tuple[HlfOpnCandidateGene, ...]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (FileNotFoundError, OSError, json.JSONDecodeError, UnicodeDecodeError):
        return ()
    if not isinstance(payload, dict) or payload.get("schema_version") != HLF_OPN_CANDIDATE_SCHEMA_VERSION:
        return ()
    rows = payload.get("rows")
    if not isinstance(rows, list):
        return ()
    return tuple(_row_from_payload(row) for row in rows if isinstance(row, dict))


def _candidate_row(
    *,
    entry: Any,
    target_context: str,
    gene_id: str,
    source_gene_resolution: str,
    full_model_row: dict[str, object] | None,
    standard_name_row: dict[str, object] | None,
    homology: GeneHomologyEvidence | None,
) -> HlfOpnCandidateGene:
    source_reactions = _source_reactions(entry)
    recommended = _recommended_intervention(entry)
    status = _operability_status(full_model_row, recommended, gene_id)
    model_operable = status in {MODEL_KO_EXECUTABLE, MODEL_OE_PROXY_EXECUTABLE}
    standard = standard_name_row or {}
    full = full_model_row or {}
    warnings = _warnings(
        target_context=target_context,
        gene_id=gene_id,
        entry=entry,
        status=status,
        recommended_intervention=recommended,
        source_reactions=source_reactions,
        homology=homology,
    )
    return HlfOpnCandidateGene(
        target_context=target_context,
        gene_id=gene_id,
        candidate_role=_candidate_role(_entry_value(entry, "category")),
        evidence_type=_evidence_type(entry, homology),
        evidence_confidence=_evidence_confidence(entry, full, homology, model_operable),
        model_operable=model_operable,
        recommended_intervention=recommended,
        reason=_reason(entry, target_context, gene_id, status),
        warnings=warnings,
        operability_status=status,
        model_operability_label=MODEL_OPERABLE if model_operable else NOT_MODEL_OPERABLE,
        source_common_name=_entry_value(entry, "common_name"),
        source_category=_entry_value(entry, "category"),
        display_name=_text(standard.get("display_name") or full.get("display_name") or gene_id),
        standard_symbol=_text(standard.get("standard_symbol") or full.get("standard_gene_symbol")),
        protein_name=_text(standard.get("protein_name") or full.get("protein_name")),
        external_ids=_dict_strings(standard.get("external_ids") or full.get("external_ids")),
        annotation_sources=_tuple_strings(standard.get("annotation_sources") or full.get("evidence_sources")),
        source_reaction_ids=source_reactions,
        executable_ko_reactions=tuple(_tuple_strings(full.get("inactive_reactions_if_ko")) if status == MODEL_KO_EXECUTABLE else ()),
        executable_oe_proxy_reactions=tuple(
            _tuple_strings(full.get("oe_executable_reactions") or full.get("affected_reactions"))
            if status == MODEL_OE_PROXY_EXECUTABLE
            else ()
        ),
        review_reactions=source_reactions if not model_operable else (),
        homology_review_status=_text(getattr(homology, "homology_review_status", "")),
        rule_transfer_status=_text(getattr(homology, "rule_transfer_status", "")),
        homology_query_symbol=_text(getattr(homology, "query_symbol", "")),
        source_gene_resolution=source_gene_resolution,
    )


def _target_contexts_for_category(category: str) -> tuple[str, ...]:
    if category in {"二硫键 (DSB)", "N-糖基化"}:
        return (TARGET_HLF,)
    if category == "O-糖基化":
        return (TARGET_OPN,)
    return (TARGET_SHARED,)


def _candidate_role(category: str) -> str:
    role_by_category = {
        "ER 转运": "er_translocation",
        "ER 折叠与分子伴侣": "er_folding_chaperone",
        "二硫键 (DSB)": "disulfide_bond_folding",
        "N-糖基化": "n_glycosylation_processing",
        "O-糖基化": "o_glycosylation_processing",
        "COPII 囊泡转运": "anterograde_vesicle_trafficking",
        "COPI 逆向转运": "retrograde_vesicle_trafficking",
        "错误折叠与 ERAD": "erad_quality_control",
        "蛋白酶体与降解": "protease_degradation_control",
        "液泡/内体分选（竞争性分流）": "vacuolar_sorting_review",
        "胞吐与分泌": "exocytosis",
        "GPI 锚定加工": "gpi_processing_review",
        "通用/其他": "global_capacity_proxy",
    }
    return role_by_category.get(category, "manual_review")


def _recommended_intervention(entry: Any) -> str:
    value = _entry_value(entry, "intervention").upper()
    if value in {"KO", "OE"}:
        return value
    return "review_only"


def _operability_status(full_model_row: dict[str, object] | None, recommended: str, gene_id: str) -> str:
    if not gene_id:
        return UNRESOLVED_NAME
    if full_model_row is None:
        return NOT_IN_MODEL
    ko_status = _text(full_model_row.get("ko_support_status"))
    oe_status = _text(full_model_row.get("oe_support_status"))
    if recommended == "KO" and ko_status == "ko_runnable_gpr_gene_deletion":
        return MODEL_KO_EXECUTABLE
    if recommended == "OE" and oe_status == "oe_runnable_reaction_proxy":
        return MODEL_OE_PROXY_EXECUTABLE
    if recommended == "OE" and oe_status == "oe_explain_only_complex_subunit":
        return MODEL_EXPLAIN_ONLY_COMPLEX_SUBUNIT
    if _tuple_strings(full_model_row.get("affected_reactions")):
        return MANUAL_REVIEW_REQUIRED
    return MANUAL_REVIEW_REQUIRED


def _evidence_type(entry: Any, homology: GeneHomologyEvidence | None) -> str:
    evidence = _entry_value(entry, "evidence")
    if any(token in evidence for token in ("已报道", "PMID", "DOI", "常用 KO", "默认 KO")):
        return "phenotype_or_literature_evidence"
    if homology is not None:
        return "homology_auxiliary"
    if _source_reactions(entry):
        return "model_reaction_proxy"
    return "curated_review"


def _evidence_confidence(
    entry: Any,
    full_model_row: dict[str, object],
    homology: GeneHomologyEvidence | None,
    model_operable: bool,
) -> str:
    evidence = _entry_value(entry, "evidence")
    if model_operable and any(token in evidence for token in ("已报道", "PMID", "DOI", "常用 KO", "默认 KO")):
        return "curated_evidence_model_executable"
    if model_operable:
        return "model_executable_annotation_supported"
    if homology is not None and getattr(homology, "is_rbh", False):
        return "homology_supported_review_required"
    if _source_reactions(entry):
        return "reaction_proxy_review_required"
    return _text(full_model_row.get("evidence_confidence")) or "manual_review_required"


def _reason(entry: Any, target_context: str, gene_id: str, status: str) -> str:
    common_name = _entry_value(entry, "common_name")
    category = _entry_value(entry, "category")
    evidence = _entry_value(entry, "evidence")
    gene_part = gene_id if gene_id else "unresolved Pichia gene_id"
    return f"{target_context} candidate from {category}: {common_name} -> {gene_part}; {status}; {evidence}"


def _warnings(
    *,
    target_context: str,
    gene_id: str,
    entry: Any,
    status: str,
    recommended_intervention: str,
    source_reactions: tuple[str, ...],
    homology: GeneHomologyEvidence | None,
) -> tuple[str, ...]:
    warnings: list[str] = []
    if target_context == TARGET_HLF:
        warnings.append(
            "hLF target parameters are draft; candidate rationale must be reviewed after sequence/PTM confirmation."
        )
    if homology is not None:
        warnings.append("SCE homology/RBH is auxiliary evidence only and is not phenotype evidence.")
    curated_gene_id = _entry_value(entry, "gene_id")
    homology_gene_id = _text(getattr(homology, "pichia_gene_id", "") if homology else "")
    if curated_gene_id and homology_gene_id and curated_gene_id != homology_gene_id:
        warnings.append(
            f"Curated gene_id {curated_gene_id} differs from homology candidate {homology_gene_id}; manual review required."
        )
    if status == MODEL_OE_PROXY_EXECUTABLE and recommended_intervention == "OE":
        warnings.append("OE is represented as a reaction-level proxy, not a gene-level overexpression simulation.")
    if status in {NOT_IN_MODEL, UNRESOLVED_NAME, MANUAL_REVIEW_REQUIRED, MODEL_EXPLAIN_ONLY_COMPLEX_SUBUNIT}:
        warnings.append(f"{NOT_MODEL_OPERABLE}: {status}.")
    if source_reactions and status not in {MODEL_KO_EXECUTABLE, MODEL_OE_PROXY_EXECUTABLE}:
        warnings.append("Curated secretion reaction proxy is retained for review and is not executable gene evidence.")
    if not gene_id:
        warnings.append("No Pichia gene_id could be resolved from curated gene_id or offline homology cache.")
    return tuple(dict.fromkeys(warnings))


def _homology_for_entry(
    entry: Any,
    homology_evidence_by_gene: dict[str, GeneHomologyEvidence] | None,
) -> GeneHomologyEvidence | None:
    if not homology_evidence_by_gene:
        return None
    common_name = _entry_value(entry, "common_name")
    return homology_evidence_for_gene(common_name, homology_evidence_by_gene, aliases=_common_name_aliases(common_name))


def _resolve_candidate_gene_id(entry: Any, homology: GeneHomologyEvidence | None) -> tuple[str, str]:
    curated_gene_id = _entry_value(entry, "gene_id")
    if curated_gene_id:
        return curated_gene_id, "curated_gene_id"
    if homology is not None:
        model_gene_id = _text(getattr(homology, "pichia_model_gene_id", ""))
        if model_gene_id:
            return model_gene_id, "homology_model_gene_id"
        pichia_gene_id = _text(getattr(homology, "pichia_gene_id", ""))
        if pichia_gene_id:
            return pichia_gene_id, "homology_pichia_gene_id"
    return "", "unresolved_name"


def _source_reactions(entry: Any) -> tuple[str, ...]:
    return tuple(
        dict.fromkeys(
            value
            for value in (
                _entry_value(entry, "ko_reaction_id"),
                _entry_value(entry, "oe_reaction_id"),
            )
            if value
        )
    )


def _common_name_aliases(common_name: str) -> tuple[str, ...]:
    separators = ("/", "（", "(", "+", " ")
    aliases = {common_name}
    current = [common_name]
    for separator in separators:
        next_values: list[str] = []
        for value in current:
            parts = [part.strip(" ）)") for part in value.split(separator) if part.strip(" ）)")]
            next_values.extend(parts)
            aliases.update(parts)
        current = next_values or current
    return tuple(alias for alias in aliases if alias)


def _normalize_target_context(value: str | None) -> str | None:
    if value is None or str(value).strip() == "":
        return None
    text = str(value).strip()
    normalized = text.lower()
    if normalized == "hlf":
        return TARGET_HLF
    if normalized == "opn":
        return TARGET_OPN
    if normalized == "shared":
        return TARGET_SHARED
    raise ValueError(f"Unknown target_context: {value}")


def _row_from_payload(payload: dict[str, object]) -> HlfOpnCandidateGene:
    return HlfOpnCandidateGene(
        target_context=_text(payload.get("target_context")),
        gene_id=_text(payload.get("gene_id")),
        candidate_role=_text(payload.get("candidate_role")),
        evidence_type=_text(payload.get("evidence_type")),
        evidence_confidence=_text(payload.get("evidence_confidence")),
        model_operable=bool(payload.get("model_operable", False)),
        recommended_intervention=_text(payload.get("recommended_intervention")),
        reason=_text(payload.get("reason")),
        warnings=_tuple_strings(payload.get("warnings")),
        operability_status=_text(payload.get("operability_status")),
        model_operability_label=_text(payload.get("model_operability_label")),
        source_common_name=_text(payload.get("source_common_name")),
        source_category=_text(payload.get("source_category")),
        display_name=_text(payload.get("display_name")),
        standard_symbol=_text(payload.get("standard_symbol")),
        protein_name=_text(payload.get("protein_name")),
        external_ids=_dict_strings(payload.get("external_ids")),
        annotation_sources=_tuple_strings(payload.get("annotation_sources")),
        source_reaction_ids=_tuple_strings(payload.get("source_reaction_ids")),
        executable_ko_reactions=_tuple_strings(payload.get("executable_ko_reactions")),
        executable_oe_proxy_reactions=_tuple_strings(payload.get("executable_oe_proxy_reactions")),
        review_reactions=_tuple_strings(payload.get("review_reactions")),
        homology_review_status=_text(payload.get("homology_review_status")),
        rule_transfer_status=_text(payload.get("rule_transfer_status")),
        homology_query_symbol=_text(payload.get("homology_query_symbol")),
        source_gene_resolution=_text(payload.get("source_gene_resolution")),
    )


def _candidate_sort_key(row: HlfOpnCandidateGene) -> tuple[int, str, str, str]:
    target_rank = {TARGET_HLF: 0, TARGET_OPN: 1, TARGET_SHARED: 2}.get(row.target_context, 3)
    operable_rank = "0" if row.model_operable else "1"
    return (target_rank, operable_rank, row.source_category, row.source_common_name)


def _count_by(rows: tuple[HlfOpnCandidateGene, ...], field_name: str) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in rows:
        value = _text(getattr(row, field_name))
        counts[value] = counts.get(value, 0) + 1
    return dict(sorted(counts.items()))


def _dedupe(values: Any) -> list[str]:
    return [item for item in dict.fromkeys(_text(value) for value in values if _text(value))]


def _entry_value(entry: Any, name: str) -> str:
    if isinstance(entry, dict):
        return _text(entry.get(name))
    return _text(getattr(entry, name, ""))


def _text(value: object) -> str:
    return str(value or "").strip()


def _dict_strings(value: object) -> dict[str, str]:
    if not isinstance(value, dict):
        return {}
    return {str(key): str(item) for key, item in value.items() if str(key).strip() and str(item).strip()}


def _tuple_strings(value: object) -> tuple[str, ...]:
    if value is None:
        return ()
    if isinstance(value, str):
        return (value,) if value.strip() else ()
    if isinstance(value, (list, tuple, set)):
        return tuple(str(item).strip() for item in value if str(item).strip())
    return ()


__all__ = [
    "HLF_OPN_CANDIDATE_SCHEMA_VERSION",
    "HlfOpnCandidateGene",
    "MODEL_EXPLAIN_ONLY_COMPLEX_SUBUNIT",
    "MODEL_KO_EXECUTABLE",
    "MODEL_OE_PROXY_EXECUTABLE",
    "NOT_IN_MODEL",
    "UNRESOLVED_NAME",
    "MANUAL_REVIEW_REQUIRED",
    "build_hlf_opn_candidate_gene_rows",
    "executable_inputs_for_hlf_opn_candidates",
    "filter_hlf_opn_candidate_gene_rows",
    "load_hlf_opn_candidate_gene_cache",
    "summarize_hlf_opn_candidate_gene_rows",
    "write_hlf_opn_candidate_gene_cache",
]
