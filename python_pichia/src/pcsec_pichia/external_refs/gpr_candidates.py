from __future__ import annotations

from dataclasses import replace
from typing import Iterable, Mapping

from pcsec_pichia.external_refs.schema import (
    ExternalGeneFunctionEvidence,
    ExternalGprCandidateEvidence,
    ExternalReactionAssociation,
)


EXTERNAL_GPR_CANDIDATE = "external_gpr_candidate"
MODEL_GPR_CONFIRMED = "model_gpr_confirmed"
REACTION_MAPPING_REQUIRED = "reaction_mapping_required"
GENE_MAPPING_REQUIRED = "gene_mapping_required"
CONFLICTING_GPR_SOURCES = "conflicting_gpr_sources"
SOURCE_RULE_MISSING = "source_rule_missing"
MODEL_REACTION_MAPPED = "model_reaction_mapped"
MODEL_GENE_MAPPED = "model_gene_mapped"
EXTERNAL_REACTION_ONLY = "external_reaction_only"
EXTERNAL_GENE_RULE_ONLY = "external_gene_rule_only"


def build_external_gpr_candidates(
    *,
    pichia_gene_id: str,
    query_gene_id: str,
    gene_function_evidence: Iterable[ExternalGeneFunctionEvidence],
    reaction_associations: Iterable[ExternalReactionAssociation],
    current_model_reaction_ids: Iterable[str],
    reaction_crosswalk: Mapping[str, str] | None = None,
    current_model_gene_ids: Iterable[str] = (),
) -> tuple[ExternalGprCandidateEvidence, ...]:
    """Build candidate-only external GPR evidence without mutating the Pichia model."""

    current_reactions = tuple(current_model_reaction_ids)
    current_genes = tuple(current_model_gene_ids)
    crosswalk = dict(reaction_crosswalk or {})
    supporting_gene_evidence = tuple(
        dict.fromkeys(item.cache_key for item in gene_function_evidence)
    )
    candidates = [
        _candidate_from_association(
            pichia_gene_id=pichia_gene_id,
            query_gene_id=query_gene_id,
            association=association,
            current_model_reaction_ids=current_reactions,
            reaction_crosswalk=crosswalk,
            current_model_gene_ids=current_genes,
            supporting_gene_evidence=supporting_gene_evidence,
        )
        for association in reaction_associations
    ]
    return tuple(sorted(_flag_conflicting_rules(candidates), key=lambda item: item.cache_key))


def classify_gpr_transfer_status(
    *,
    gene_mapping_status: str,
    reaction_mapping_status: str,
    source_gene_rule: str | None,
    mapped_model_reaction_id: str | None,
    in_current_model_gene_index: bool,
) -> tuple[str, tuple[str, ...]]:
    """Classify external GPR transfer without treating it as current model truth."""

    if not str(source_gene_rule or "").strip():
        return SOURCE_RULE_MISSING, ("external source reaction lacks a gene rule",)
    if reaction_mapping_status != MODEL_REACTION_MAPPED or not mapped_model_reaction_id:
        return REACTION_MAPPING_REQUIRED, (
            "external reaction is not mapped to a current Pichia model reaction",
        )
    if gene_mapping_status != MODEL_GENE_MAPPED or not in_current_model_gene_index:
        return GENE_MAPPING_REQUIRED, (
            "external gene rule is not mapped to a current Pichia model gene",
        )
    return MODEL_GPR_CONFIRMED, ()


def _candidate_from_association(
    *,
    pichia_gene_id: str,
    query_gene_id: str,
    association: ExternalReactionAssociation,
    current_model_reaction_ids: tuple[str, ...],
    reaction_crosswalk: Mapping[str, str],
    current_model_gene_ids: tuple[str, ...],
    supporting_gene_evidence: tuple[str, ...],
) -> ExternalGprCandidateEvidence:
    mapped_reaction = _mapped_reaction_id(
        association,
        current_model_reaction_ids=current_model_reaction_ids,
        reaction_crosswalk=reaction_crosswalk,
    )
    reaction_status = MODEL_REACTION_MAPPED if mapped_reaction else EXTERNAL_REACTION_ONLY
    mapped_genes = _mapped_gene_ids(
        pichia_gene_id=pichia_gene_id,
        query_gene_id=query_gene_id,
        association=association,
        current_model_gene_ids=current_model_gene_ids,
    )
    gene_status = MODEL_GENE_MAPPED if mapped_genes else EXTERNAL_GENE_RULE_ONLY
    in_current_model_gene_index = bool(
        mapped_genes
        or (current_model_gene_ids and _contains_token(current_model_gene_ids, pichia_gene_id))
    )
    status, blocking_reasons = classify_gpr_transfer_status(
        gene_mapping_status=gene_status,
        reaction_mapping_status=reaction_status,
        source_gene_rule=association.gene_rule,
        mapped_model_reaction_id=mapped_reaction,
        in_current_model_gene_index=in_current_model_gene_index,
    )
    return ExternalGprCandidateEvidence(
        provenance=association.provenance,
        external_model_id=association.external_model_id,
        external_reaction_id=association.external_reaction_id,
        external_gene_rule=association.gene_rule,
        candidate_status=status,
        pichia_gene_id=pichia_gene_id,
        query_gene_id=query_gene_id,
        mapped_pichia_reaction_id=mapped_reaction,
        mapped_pichia_gene_ids=mapped_genes,
        gene_mapping_status=gene_status,
        reaction_mapping_status=reaction_status,
        gpr_transfer_status=status,
        confidence=_confidence_for_status(status),
        supporting_gene_evidence=supporting_gene_evidence,
        blocking_reasons=blocking_reasons,
        mapping_warnings=tuple(dict.fromkeys((*association.provenance.warnings, *blocking_reasons))),
    )


def _mapped_reaction_id(
    association: ExternalReactionAssociation,
    *,
    current_model_reaction_ids: tuple[str, ...],
    reaction_crosswalk: Mapping[str, str],
) -> str | None:
    current_by_token = {_normalize_token(reaction_id): reaction_id for reaction_id in current_model_reaction_ids}
    for candidate in (
        association.mapped_pichia_reaction_id,
        reaction_crosswalk.get(association.external_reaction_id),
        reaction_crosswalk.get(_normalize_token(association.external_reaction_id)),
        association.external_reaction_id,
    ):
        token = _normalize_token(candidate)
        if token and token in current_by_token:
            return current_by_token[token]
    return None


def _mapped_gene_ids(
    *,
    pichia_gene_id: str,
    query_gene_id: str,
    association: ExternalReactionAssociation,
    current_model_gene_ids: tuple[str, ...],
) -> tuple[str, ...]:
    mapped = tuple(
        gene_id
        for gene_id in association.mapped_pichia_gene_ids
        if _normalize_token(gene_id) == _normalize_token(pichia_gene_id)
    )
    if mapped:
        return mapped
    if (
        current_model_gene_ids
        and _contains_token(current_model_gene_ids, pichia_gene_id)
        and _source_rule_contains_query(association, query_gene_id)
    ):
        return (pichia_gene_id,)
    return ()


def _flag_conflicting_rules(
    candidates: Iterable[ExternalGprCandidateEvidence],
) -> tuple[ExternalGprCandidateEvidence, ...]:
    resolved = tuple(candidates)
    rules_by_reaction: dict[str, set[str]] = {}
    for candidate in resolved:
        key = candidate.mapped_pichia_reaction_id or candidate.external_reaction_id
        rules_by_reaction.setdefault(_normalize_token(key), set()).add(
            _normalize_rule(candidate.external_gene_rule)
        )
    conflict_keys = {
        key
        for key, rules in rules_by_reaction.items()
        if len({rule for rule in rules if rule}) > 1
    }
    if not conflict_keys:
        return resolved
    flagged: list[ExternalGprCandidateEvidence] = []
    for candidate in resolved:
        key = _normalize_token(candidate.mapped_pichia_reaction_id or candidate.external_reaction_id)
        if key in conflict_keys:
            reason = "conflicting external GPR rules"
            reasons = tuple(dict.fromkeys((*candidate.blocking_reasons, reason)))
            flagged.append(
                replace(
                    candidate,
                    candidate_status=CONFLICTING_GPR_SOURCES,
                    gpr_transfer_status=CONFLICTING_GPR_SOURCES,
                    confidence="manual_review_required",
                    blocking_reasons=reasons,
                    mapping_warnings=tuple(dict.fromkeys((*candidate.mapping_warnings, reason))),
                )
            )
        else:
            flagged.append(candidate)
    return tuple(flagged)


def _source_rule_contains_query(
    association: ExternalReactionAssociation,
    query_gene_id: str,
) -> bool:
    query = _normalize_token(query_gene_id)
    source_tokens = {
        _normalize_token(value)
        for value in (*association.external_gene_ids, *_split_rule_tokens(association.gene_rule))
    }
    return bool(query and query in source_tokens)


def _confidence_for_status(status: str) -> str:
    if status == MODEL_GPR_CONFIRMED:
        return "mapped_external_gpr"
    if status == EXTERNAL_GPR_CANDIDATE:
        return "external_candidate"
    return "manual_review_required"


def _contains_token(values: Iterable[str], query: str) -> bool:
    token = _normalize_token(query)
    return bool(token and token in {_normalize_token(value) for value in values})


def _split_rule_tokens(value: object) -> tuple[str, ...]:
    text = str(value or "")
    for sep in ("(", ")", "/", ",", ";", "|", "+"):
        text = text.replace(sep, " ")
    return tuple(
        part
        for part in text.split()
        if part.strip().lower() not in {"and", "or"}
    )


def _normalize_rule(value: object) -> str:
    return " ".join(_normalize_token(part) for part in _split_rule_tokens(value))


def _normalize_token(value: object) -> str:
    return str(value or "").strip().upper().replace(" ", "").replace("_", "").replace("-", "")


__all__ = [
    "CONFLICTING_GPR_SOURCES",
    "EXTERNAL_GPR_CANDIDATE",
    "GENE_MAPPING_REQUIRED",
    "MODEL_GENE_MAPPED",
    "MODEL_GPR_CONFIRMED",
    "MODEL_REACTION_MAPPED",
    "REACTION_MAPPING_REQUIRED",
    "SOURCE_RULE_MISSING",
    "build_external_gpr_candidates",
    "classify_gpr_transfer_status",
]
