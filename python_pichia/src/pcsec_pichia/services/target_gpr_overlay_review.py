from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

from pcsec_pichia.services.gene_rule_overlay import GeneRuleEvidence, HIGH_CONFIDENCE


HLF_OPN_OVERLAY_REVIEW_SCHEMA_VERSION = 1

CANDIDATE_GPR_OVERLAY_REVIEW = "candidate_gpr_overlay_review"
MODEL_EXPANSION_REQUIRED = "model_expansion_required"
MANUAL_REVIEW_REQUIRED = "manual_review_required"


@dataclass(frozen=True)
class HlfOpnGprOverlayReview:
    target_context: str
    gene_id: str
    source_common_name: str
    source_candidate_gene_id: str
    candidate_role: str
    recommended_intervention: str
    review_status: str
    evidence_confidence: str
    evidence_sources: tuple[str, ...]
    external_ids: dict[str, str]
    protein_name: str
    target_reaction_ids: tuple[str, ...]
    existing_model_reaction_ids: tuple[str, ...]
    missing_model_reaction_ids: tuple[str, ...]
    rule_status: str
    recommended_action: str
    risk: str
    warnings: tuple[str, ...] = field(default_factory=tuple)
    homology_review_status: str = ""
    rule_transfer_status: str = ""
    source_cache: str = "gene_rule_evidence_cache"

    def to_dict(self) -> dict[str, Any]:
        payload = asdict(self)
        for key in (
            "evidence_sources",
            "target_reaction_ids",
            "existing_model_reaction_ids",
            "missing_model_reaction_ids",
            "warnings",
        ):
            payload[key] = list(payload[key])
        payload["external_ids"] = dict(self.external_ids)
        return payload


def build_hlf_opn_gpr_overlay_review_rows(
    *,
    candidate_rows: list[dict[str, object]] | tuple[dict[str, object], ...],
    gene_rule_evidence_by_name: dict[str, GeneRuleEvidence],
    model_reaction_ids: set[str] | frozenset[str] | tuple[str, ...] | list[str],
) -> tuple[HlfOpnGprOverlayReview, ...]:
    """Review only high-value hLF/OPN model-external candidates for GPR overlay potential."""

    evidence_by_name = _normalise_evidence(gene_rule_evidence_by_name)
    reaction_ids = {str(item) for item in model_reaction_ids}
    rows: list[HlfOpnGprOverlayReview] = []
    seen: set[tuple[str, str, str]] = set()
    for candidate in candidate_rows:
        if bool(candidate.get("model_operable")):
            continue
        evidence = _evidence_for_candidate(candidate, evidence_by_name)
        if evidence is None:
            continue
        row = _review_row(candidate, evidence, reaction_ids)
        key = (row.target_context, row.gene_id, row.source_common_name)
        if key in seen:
            continue
        seen.add(key)
        rows.append(row)
    return tuple(sorted(rows, key=_sort_key))


def filter_hlf_opn_gpr_overlay_review_rows(
    rows: tuple[HlfOpnGprOverlayReview, ...] | list[HlfOpnGprOverlayReview],
    *,
    target_context: str | None = None,
    include_shared: bool = True,
) -> tuple[HlfOpnGprOverlayReview, ...]:
    context = _normalise_target_context(target_context)
    if context is None:
        return tuple(rows)
    allowed = {context}
    if include_shared and context in {"hLF", "OPN"}:
        allowed.add("shared")
    return tuple(row for row in rows if row.target_context in allowed)


def summarize_hlf_opn_gpr_overlay_review_rows(
    rows: tuple[HlfOpnGprOverlayReview, ...] | list[HlfOpnGprOverlayReview],
) -> dict[str, Any]:
    all_rows = tuple(rows)
    return {
        "schema_version": HLF_OPN_OVERLAY_REVIEW_SCHEMA_VERSION,
        "total_review_rows": len(all_rows),
        "review_status_counts": _count_by(all_rows, "review_status"),
        "target_context_counts": _count_by(all_rows, "target_context"),
        "high_confidence_review_count": sum(1 for row in all_rows if row.evidence_confidence == HIGH_CONFIDENCE),
        "candidate_gpr_overlay_review_count": sum(
            1 for row in all_rows if row.review_status == CANDIDATE_GPR_OVERLAY_REVIEW
        ),
        "model_expansion_required_count": sum(1 for row in all_rows if row.review_status == MODEL_EXPANSION_REQUIRED),
        "manual_review_required_count": sum(1 for row in all_rows if row.review_status == MANUAL_REVIEW_REQUIRED),
    }


def write_hlf_opn_gpr_overlay_review_cache(
    rows: tuple[HlfOpnGprOverlayReview, ...] | list[HlfOpnGprOverlayReview],
    output_path: Path,
) -> None:
    payload = {
        "schema_version": HLF_OPN_OVERLAY_REVIEW_SCHEMA_VERSION,
        "summary": summarize_hlf_opn_gpr_overlay_review_rows(tuple(rows)),
        "rows": [row.to_dict() for row in rows],
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")


def load_hlf_opn_gpr_overlay_review_cache(path: Path) -> tuple[HlfOpnGprOverlayReview, ...]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (FileNotFoundError, OSError, json.JSONDecodeError, UnicodeDecodeError):
        return ()
    if not isinstance(payload, dict) or payload.get("schema_version") != HLF_OPN_OVERLAY_REVIEW_SCHEMA_VERSION:
        return ()
    rows = payload.get("rows")
    if not isinstance(rows, list):
        return ()
    return tuple(_row_from_payload(row) for row in rows if isinstance(row, dict))


def _review_row(
    candidate: dict[str, object],
    evidence: GeneRuleEvidence,
    model_reaction_ids: set[str],
) -> HlfOpnGprOverlayReview:
    target_reactions = tuple(str(item) for item in evidence.target_reaction_ids if str(item))
    existing = tuple(reaction_id for reaction_id in target_reactions if reaction_id in model_reaction_ids)
    missing = tuple(reaction_id for reaction_id in target_reactions if reaction_id not in model_reaction_ids)
    review_status = _review_status(evidence, existing, missing)
    gene_id = str(evidence.candidate_locus_tag or candidate.get("gene_id") or "").strip()
    warnings = _warnings(candidate, evidence, review_status, existing, missing, gene_id)
    return HlfOpnGprOverlayReview(
        target_context=str(candidate.get("target_context") or ""),
        gene_id=gene_id,
        source_common_name=str(candidate.get("source_common_name") or evidence.common_name or ""),
        source_candidate_gene_id=str(candidate.get("gene_id") or ""),
        candidate_role=str(candidate.get("candidate_role") or ""),
        recommended_intervention=str(candidate.get("recommended_intervention") or "review_only"),
        review_status=review_status,
        evidence_confidence=evidence.confidence,
        evidence_sources=tuple(evidence.evidence_sources),
        external_ids=dict(evidence.external_ids or {}),
        protein_name=evidence.protein_name,
        target_reaction_ids=target_reactions,
        existing_model_reaction_ids=existing,
        missing_model_reaction_ids=missing,
        rule_status=evidence.rule_status,
        recommended_action=evidence.recommended_action,
        risk=_risk_text(evidence, review_status, existing, missing),
        warnings=warnings,
        homology_review_status=str(candidate.get("homology_review_status") or ""),
        rule_transfer_status=str(candidate.get("rule_transfer_status") or ""),
    )


def _review_status(evidence: GeneRuleEvidence, existing: tuple[str, ...], missing: tuple[str, ...]) -> str:
    if evidence.confidence != HIGH_CONFIDENCE:
        return MANUAL_REVIEW_REQUIRED
    if existing:
        return CANDIDATE_GPR_OVERLAY_REVIEW
    if missing or evidence.target_reaction_ids:
        return MODEL_EXPANSION_REQUIRED
    return MANUAL_REVIEW_REQUIRED


def _risk_text(
    evidence: GeneRuleEvidence,
    review_status: str,
    existing: tuple[str, ...],
    missing: tuple[str, ...],
) -> str:
    if review_status == CANDIDATE_GPR_OVERLAY_REVIEW:
        return (
            "High-confidence external locus evidence maps to existing model reaction(s), "
            "but overlay remains review-only and must not be written to formal GPR automatically."
        )
    if review_status == MODEL_EXPANSION_REQUIRED:
        return (
            "External evidence does not map to a current model reaction; model expansion would be required "
            f"before simulation use. Missing reactions: {', '.join(missing) or 'none'}."
        )
    return (
        "Evidence is below the executable overlay threshold or incomplete; keep this candidate in manual review. "
        f"rule_status={evidence.rule_status or 'unknown'}; existing reactions={', '.join(existing) or 'none'}."
    )


def _warnings(
    candidate: dict[str, object],
    evidence: GeneRuleEvidence,
    review_status: str,
    existing: tuple[str, ...],
    missing: tuple[str, ...],
    gene_id: str,
) -> tuple[str, ...]:
    warnings = [
        "Review row only: do not auto-write formal GPR or enable simulation without explicit review.",
        "External annotation, EC/KEGG/UniProt, and homology evidence are not phenotype evidence.",
    ]
    source_gene_id = str(candidate.get("gene_id") or "").strip()
    if source_gene_id and gene_id and source_gene_id != gene_id:
        warnings.append(f"Candidate gene_id {source_gene_id} differs from external locus {gene_id}.")
    if review_status == CANDIDATE_GPR_OVERLAY_REVIEW:
        warnings.append("Candidate may enter candidate_gpr_overlay_review, not executable KO/OE inputs.")
    elif review_status == MODEL_EXPANSION_REQUIRED:
        warnings.append("No current model reaction is available for this evidence row.")
    else:
        warnings.append("Keep as manual_review_required until stronger evidence is available.")
    if not existing and not missing:
        warnings.append("No target reaction was supplied by the external evidence cache.")
    return tuple(dict.fromkeys(warnings))


def _evidence_for_candidate(
    candidate: dict[str, object],
    evidence_by_name: dict[str, GeneRuleEvidence],
) -> GeneRuleEvidence | None:
    for alias in _candidate_aliases(candidate):
        evidence = evidence_by_name.get(alias)
        if evidence is not None:
            return evidence
    return None


def _candidate_aliases(candidate: dict[str, object]) -> tuple[str, ...]:
    values = (
        candidate.get("source_common_name"),
        candidate.get("homology_query_symbol"),
        candidate.get("standard_symbol"),
        candidate.get("display_name"),
    )
    aliases: list[str] = []
    for value in values:
        text = str(value or "").strip()
        if not text:
            continue
        aliases.append(_normalise_name(text))
    for value in values:
        text = str(value or "").strip()
        if not text:
            continue
        for part in _split_name(text):
            aliases.append(_normalise_name(part))
    return tuple(dict.fromkeys(alias for alias in aliases if alias))


def _normalise_evidence(evidence_by_name: dict[str, GeneRuleEvidence]) -> dict[str, GeneRuleEvidence]:
    result: dict[str, GeneRuleEvidence] = {}
    for key, record in evidence_by_name.items():
        names = [key, record.common_name]
        for name in names:
            normalised = _normalise_name(name)
            if normalised:
                result[normalised] = record
    return result


def _split_name(value: str) -> tuple[str, ...]:
    pieces = [value]
    seen = {value}
    for separator in ("/", "（", "(", "+", " "):
        next_pieces: list[str] = []
        for piece in pieces:
            for part in (item.strip(" ）)") for item in piece.split(separator) if item.strip(" ）)")):
                if part in seen:
                    continue
                seen.add(part)
                next_pieces.append(part)
        pieces.extend(next_pieces)
    return tuple(pieces)


def _normalise_name(value: object) -> str:
    return str(value or "").strip().upper()


def _normalise_target_context(value: str | None) -> str | None:
    if value is None or str(value).strip() == "":
        return None
    text = str(value).strip()
    lowered = text.lower()
    if lowered == "hlf":
        return "hLF"
    if lowered == "opn":
        return "OPN"
    if lowered == "shared":
        return "shared"
    raise ValueError(f"Unknown target_context: {value}")


def _row_from_payload(payload: dict[str, object]) -> HlfOpnGprOverlayReview:
    return HlfOpnGprOverlayReview(
        target_context=str(payload.get("target_context") or ""),
        gene_id=str(payload.get("gene_id") or ""),
        source_common_name=str(payload.get("source_common_name") or ""),
        source_candidate_gene_id=str(payload.get("source_candidate_gene_id") or ""),
        candidate_role=str(payload.get("candidate_role") or ""),
        recommended_intervention=str(payload.get("recommended_intervention") or ""),
        review_status=str(payload.get("review_status") or ""),
        evidence_confidence=str(payload.get("evidence_confidence") or ""),
        evidence_sources=_tuple_strings(payload.get("evidence_sources")),
        external_ids=_dict_strings(payload.get("external_ids")),
        protein_name=str(payload.get("protein_name") or ""),
        target_reaction_ids=_tuple_strings(payload.get("target_reaction_ids")),
        existing_model_reaction_ids=_tuple_strings(payload.get("existing_model_reaction_ids")),
        missing_model_reaction_ids=_tuple_strings(payload.get("missing_model_reaction_ids")),
        rule_status=str(payload.get("rule_status") or ""),
        recommended_action=str(payload.get("recommended_action") or ""),
        risk=str(payload.get("risk") or ""),
        warnings=_tuple_strings(payload.get("warnings")),
        homology_review_status=str(payload.get("homology_review_status") or ""),
        rule_transfer_status=str(payload.get("rule_transfer_status") or ""),
        source_cache=str(payload.get("source_cache") or "gene_rule_evidence_cache"),
    )


def _sort_key(row: HlfOpnGprOverlayReview) -> tuple[int, int, str, str]:
    status_rank = {
        CANDIDATE_GPR_OVERLAY_REVIEW: 0,
        MANUAL_REVIEW_REQUIRED: 1,
        MODEL_EXPANSION_REQUIRED: 2,
    }.get(row.review_status, 3)
    target_rank = {"hLF": 0, "OPN": 1, "shared": 2}.get(row.target_context, 3)
    return (status_rank, target_rank, row.source_common_name, row.gene_id)


def _count_by(rows: tuple[HlfOpnGprOverlayReview, ...], field_name: str) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in rows:
        value = str(getattr(row, field_name) or "")
        counts[value] = counts.get(value, 0) + 1
    return dict(sorted(counts.items()))


def _tuple_strings(value: object) -> tuple[str, ...]:
    if value is None:
        return ()
    if isinstance(value, str):
        return (value,) if value.strip() else ()
    if isinstance(value, (list, tuple, set)):
        return tuple(str(item).strip() for item in value if str(item).strip())
    return ()


def _dict_strings(value: object) -> dict[str, str]:
    if not isinstance(value, dict):
        return {}
    return {str(key): str(item) for key, item in value.items() if str(key).strip() and str(item).strip()}


__all__ = [
    "CANDIDATE_GPR_OVERLAY_REVIEW",
    "HLF_OPN_OVERLAY_REVIEW_SCHEMA_VERSION",
    "HlfOpnGprOverlayReview",
    "MANUAL_REVIEW_REQUIRED",
    "MODEL_EXPANSION_REQUIRED",
    "build_hlf_opn_gpr_overlay_review_rows",
    "filter_hlf_opn_gpr_overlay_review_rows",
    "load_hlf_opn_gpr_overlay_review_cache",
    "summarize_hlf_opn_gpr_overlay_review_rows",
    "write_hlf_opn_gpr_overlay_review_cache",
]
