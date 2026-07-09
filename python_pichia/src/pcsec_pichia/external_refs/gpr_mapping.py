from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
from typing import Iterable, Mapping

from pcsec_pichia.external_refs.cache_io import write_external_reference_cache_bundle
from pcsec_pichia.external_refs.gpr_candidates import (
    CONFLICTING_GPR_SOURCES,
    GENE_MAPPING_REQUIRED,
    MODEL_GENE_MAPPED,
    MODEL_GPR_CONFIRMED,
    MODEL_REACTION_MAPPED,
    REACTION_MAPPING_REQUIRED,
    SOURCE_RULE_MISSING,
)
from pcsec_pichia.external_refs.schema import ExternalGprCandidateEvidence


EXTERNAL_GPR_MAPPING_CANDIDATES_FILENAME = "external_gpr_candidate_evidence.jsonl"
EXTERNAL_GPR_MAPPING_MANIFEST_FILENAME = "external_gpr_mapping_manifest.json"
EXTERNAL_GPR_MAPPING_REPORT_FILENAME = "external_gpr_mapping_report.md"
MANUAL_REVIEW_REQUIRED = "manual_review_required"
NOT_IN_CURRENT_MODEL = "not_in_current_model"
MAPPED_EXTERNAL_GPR_CONFIDENCE = "mapped_external_gpr"


@dataclass(frozen=True)
class ExternalGprMappingOutputs:
    candidates_path: Path
    manifest_path: Path
    report_path: Path
    candidate_count: int
    status_counts: Mapping[str, int]


def map_external_gpr_candidates_to_model(
    candidates: Iterable[ExternalGprCandidateEvidence],
    *,
    current_model_reaction_ids: Iterable[str],
    current_model_gene_ids: Iterable[str],
    reaction_crosswalk: Mapping[str, str] | None = None,
    gene_crosswalk: Mapping[str, str] | None = None,
    homology_gene_crosswalk: Mapping[str, str] | None = None,
    external_name_gene_crosswalk: Mapping[str, str] | None = None,
) -> tuple[ExternalGprCandidateEvidence, ...]:
    """Map external GPR candidates onto current Pichia reaction/gene IDs.

    The result remains evidence-only. A confirmed mapping means the external
    rule can be reviewed against current model identifiers; it is not written
    into the current GEM and is not treated as ``model_gpr_executable``.
    """

    current_reactions = tuple(current_model_reaction_ids)
    current_genes = tuple(current_model_gene_ids)
    reaction_lookup = _normalized_lookup(current_reactions)
    gene_lookup = _normalized_lookup(current_genes)
    reaction_map = _normalized_crosswalk(reaction_crosswalk or {})
    gene_map = _normalized_crosswalk(
        {
            **dict(homology_gene_crosswalk or {}),
            **dict(external_name_gene_crosswalk or {}),
            **dict(gene_crosswalk or {}),
        }
    )
    mapped = tuple(
        _map_one_candidate(
            candidate,
            reaction_lookup=reaction_lookup,
            gene_lookup=gene_lookup,
            reaction_crosswalk=reaction_map,
            gene_crosswalk=gene_map,
        )
        for candidate in candidates
    )
    return tuple(sorted(_flag_conflicting_mapped_rules(mapped), key=lambda item: item.cache_key))


def write_external_gpr_mapping_outputs(
    candidates: Iterable[ExternalGprCandidateEvidence],
    output_dir: Path,
) -> ExternalGprMappingOutputs:
    resolved = tuple(candidates)
    manifest = write_external_reference_cache_bundle(
        resolved,
        output_dir,
        query_count=len(resolved),
        records_filename=EXTERNAL_GPR_MAPPING_CANDIDATES_FILENAME,
        manifest_filename=EXTERNAL_GPR_MAPPING_MANIFEST_FILENAME,
    )
    candidates_path = output_dir / EXTERNAL_GPR_MAPPING_CANDIDATES_FILENAME
    manifest_path = output_dir / EXTERNAL_GPR_MAPPING_MANIFEST_FILENAME
    report_path = output_dir / EXTERNAL_GPR_MAPPING_REPORT_FILENAME
    report_path.write_text(render_external_gpr_mapping_report(resolved), encoding="utf-8")
    return ExternalGprMappingOutputs(
        candidates_path=candidates_path,
        manifest_path=manifest_path,
        report_path=report_path,
        candidate_count=manifest.record_count,
        status_counts=_status_counts(resolved),
    )


def render_external_gpr_mapping_report(
    candidates: tuple[ExternalGprCandidateEvidence, ...],
) -> str:
    status_counts = _status_counts(candidates)
    lines = [
        "# External GPR Mapping Report",
        "",
        "External GPR candidates are mapped evidence only; rules are not written into the current Pichia GEM.",
        "",
        "## Status Counts",
        "",
        "| status | count |",
        "|---|---:|",
    ]
    for status, count in sorted(status_counts.items()):
        lines.append(f"| {status} | {count} |")
    lines.extend(
        [
            "",
            "## Candidates",
            "",
            "| status | source | model | external reaction | mapped reaction | mapped genes | reasons |",
            "|---|---|---|---|---|---|---|",
        ]
    )
    for candidate in candidates:
        lines.append(
            f"| {candidate.candidate_status} | {candidate.provenance.source_database} | "
            f"{candidate.external_model_id} | {candidate.external_reaction_id} | "
            f"{candidate.mapped_pichia_reaction_id or ''} | "
            f"{'; '.join(candidate.mapped_pichia_gene_ids)} | "
            f"{'; '.join(candidate.blocking_reasons)} |"
        )
    return "\n".join(lines) + "\n"


def _map_one_candidate(
    candidate: ExternalGprCandidateEvidence,
    *,
    reaction_lookup: Mapping[str, str],
    gene_lookup: Mapping[str, str],
    reaction_crosswalk: Mapping[str, str],
    gene_crosswalk: Mapping[str, str],
) -> ExternalGprCandidateEvidence:
    mapped_reaction, reaction_reasons = _mapped_reaction(candidate, reaction_lookup, reaction_crosswalk)
    mapped_genes, gene_reasons = _mapped_genes(candidate, gene_lookup, gene_crosswalk)
    status, reasons = _mapping_status(candidate, mapped_reaction, mapped_genes, reaction_reasons + gene_reasons)
    warnings = tuple(dict.fromkeys((*candidate.mapping_warnings, *reasons)))
    return replace(
        candidate,
        candidate_status=status,
        mapped_pichia_reaction_id=mapped_reaction,
        mapped_pichia_gene_ids=mapped_genes,
        gene_mapping_status=MODEL_GENE_MAPPED if mapped_genes else GENE_MAPPING_REQUIRED,
        reaction_mapping_status=MODEL_REACTION_MAPPED if mapped_reaction else REACTION_MAPPING_REQUIRED,
        gpr_transfer_status=status,
        confidence=MAPPED_EXTERNAL_GPR_CONFIDENCE if status == MODEL_GPR_CONFIRMED else MANUAL_REVIEW_REQUIRED,
        blocking_reasons=reasons,
        mapping_warnings=warnings,
    )


def _mapped_reaction(
    candidate: ExternalGprCandidateEvidence,
    reaction_lookup: Mapping[str, str],
    reaction_crosswalk: Mapping[str, str],
) -> tuple[str | None, tuple[str, ...]]:
    explicit_candidates = (
        candidate.mapped_pichia_reaction_id,
        reaction_crosswalk.get(_normalize_token(candidate.external_reaction_id)),
    )
    reasons: list[str] = []
    for value in explicit_candidates:
        token = _normalize_token(value)
        if token in reaction_lookup:
            return reaction_lookup[token], ()
        if token:
            reasons.append(f"mapped reaction is not present in current Pichia model: {value}")
    token = _normalize_token(candidate.external_reaction_id)
    if token in reaction_lookup:
        return reaction_lookup[token], ()
    return None, tuple(dict.fromkeys(reasons))


def _mapped_genes(
    candidate: ExternalGprCandidateEvidence,
    gene_lookup: Mapping[str, str],
    gene_crosswalk: Mapping[str, str],
) -> tuple[tuple[str, ...], tuple[str, ...]]:
    mapped: list[str] = []
    reasons: list[str] = []
    for value in candidate.mapped_pichia_gene_ids:
        token = _normalize_token(value)
        if token in gene_lookup:
            mapped.append(gene_lookup[token])
        elif token:
            reasons.append(f"mapped gene is not present in current Pichia model: {value}")
    for token in _candidate_gene_tokens(candidate):
        mapped_value = gene_crosswalk.get(_normalize_token(token))
        mapped_token = _normalize_token(mapped_value)
        if mapped_token in gene_lookup:
            mapped.append(gene_lookup[mapped_token])
        elif mapped_token:
            reasons.append(f"mapped gene is not present in current Pichia model: {mapped_value}")
        elif _normalize_token(token) in gene_lookup:
            mapped.append(gene_lookup[_normalize_token(token)])
    return tuple(dict.fromkeys(mapped)), tuple(dict.fromkeys(reasons))


def _mapping_status(
    candidate: ExternalGprCandidateEvidence,
    mapped_reaction: str | None,
    mapped_genes: tuple[str, ...],
    not_in_model_reasons: tuple[str, ...],
) -> tuple[str, tuple[str, ...]]:
    if candidate.candidate_status == CONFLICTING_GPR_SOURCES or candidate.gpr_transfer_status == CONFLICTING_GPR_SOURCES:
        return CONFLICTING_GPR_SOURCES, ("conflicting external GPR rules require manual review",)
    if not str(candidate.external_gene_rule or "").strip():
        return SOURCE_RULE_MISSING, ("external source reaction lacks a gene rule",)
    if candidate.candidate_status == MANUAL_REVIEW_REQUIRED or candidate.gpr_transfer_status == MANUAL_REVIEW_REQUIRED:
        return MANUAL_REVIEW_REQUIRED, candidate.blocking_reasons or ("external GPR candidate requires manual review",)
    if not_in_model_reasons:
        return NOT_IN_CURRENT_MODEL, not_in_model_reasons
    if not mapped_reaction:
        return REACTION_MAPPING_REQUIRED, ("external reaction is not mapped to a current Pichia model reaction",)
    if not mapped_genes:
        return GENE_MAPPING_REQUIRED, ("external gene rule is not mapped to a current Pichia model gene",)
    return MODEL_GPR_CONFIRMED, ()


def _flag_conflicting_mapped_rules(
    candidates: tuple[ExternalGprCandidateEvidence, ...],
) -> tuple[ExternalGprCandidateEvidence, ...]:
    rules_by_reaction: dict[str, set[str]] = {}
    for candidate in candidates:
        if candidate.mapped_pichia_reaction_id:
            rules_by_reaction.setdefault(_normalize_token(candidate.mapped_pichia_reaction_id), set()).add(
                _normalize_rule(candidate.external_gene_rule)
            )
    conflict_keys = {
        key
        for key, rules in rules_by_reaction.items()
        if len({rule for rule in rules if rule}) > 1
    }
    if not conflict_keys:
        return candidates
    flagged: list[ExternalGprCandidateEvidence] = []
    for candidate in candidates:
        key = _normalize_token(candidate.mapped_pichia_reaction_id)
        if key in conflict_keys:
            reason = "conflicting external GPR rules require manual review"
            flagged.append(
                replace(
                    candidate,
                    candidate_status=CONFLICTING_GPR_SOURCES,
                    gpr_transfer_status=CONFLICTING_GPR_SOURCES,
                    confidence=MANUAL_REVIEW_REQUIRED,
                    blocking_reasons=tuple(dict.fromkeys((*candidate.blocking_reasons, reason))),
                    mapping_warnings=tuple(dict.fromkeys((*candidate.mapping_warnings, reason))),
                )
            )
        else:
            flagged.append(candidate)
    return tuple(flagged)


def _candidate_gene_tokens(candidate: ExternalGprCandidateEvidence) -> tuple[str, ...]:
    return tuple(
        dict.fromkeys(
            token
            for token in (
                *_split_rule_tokens(candidate.external_gene_rule),
                candidate.query_gene_id or "",
                candidate.pichia_gene_id or "",
            )
            if str(token).strip()
        )
    )


def _status_counts(candidates: Iterable[ExternalGprCandidateEvidence]) -> Mapping[str, int]:
    counts: dict[str, int] = {}
    for candidate in candidates:
        counts[candidate.candidate_status] = counts.get(candidate.candidate_status, 0) + 1
    return counts


def _normalized_lookup(values: Iterable[str]) -> Mapping[str, str]:
    return {_normalize_token(value): value for value in values if _normalize_token(value)}


def _normalized_crosswalk(crosswalk: Mapping[str, str]) -> Mapping[str, str]:
    return {
        _normalize_token(key): value
        for key, value in crosswalk.items()
        if _normalize_token(key) and str(value or "").strip()
    }


def _split_rule_tokens(value: object) -> tuple[str, ...]:
    text = str(value or "")
    for sep in ("(", ")", "/", ",", ";", "|", "+"):
        text = text.replace(sep, " ")
    return tuple(part for part in text.split() if part.strip().lower() not in {"and", "or"})


def _normalize_rule(value: object) -> str:
    return " ".join(_normalize_token(part) for part in _split_rule_tokens(value))


def _normalize_token(value: object) -> str:
    return str(value or "").strip().upper().replace(" ", "").replace("_", "").replace("-", "")


__all__ = [
    "EXTERNAL_GPR_MAPPING_CANDIDATES_FILENAME",
    "EXTERNAL_GPR_MAPPING_MANIFEST_FILENAME",
    "EXTERNAL_GPR_MAPPING_REPORT_FILENAME",
    "ExternalGprMappingOutputs",
    "map_external_gpr_candidates_to_model",
    "render_external_gpr_mapping_report",
    "write_external_gpr_mapping_outputs",
]
