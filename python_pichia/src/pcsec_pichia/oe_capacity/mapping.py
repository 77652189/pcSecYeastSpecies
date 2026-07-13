from __future__ import annotations

import hashlib
import json
from collections import Counter
from pathlib import Path
from typing import Any, Mapping, Sequence

from pcsec_pichia.oe_capacity.schema import (
    ConfidenceLevel,
    EvidenceSourceType,
    GeneCapacityCatalog,
    GeneCapacityCoverage,
    GeneCapacityValidationIssue,
    GeneCapacityValidationResult,
    GeneEnzymeReactionMapping,
    GPRRole,
    OECapacityValidationError,
    OEExecutionStatus,
)
from pcsec_pichia.screens.gene_interventions import plan_gene_overexpression


def build_gene_enzyme_reaction_catalog(
    model: Any,
    metabolic: Any,
    combined: Any,
    external_evidence: Sequence[Mapping[str, Any]] | None = None,
) -> GeneCapacityCatalog:
    fingerprint = _model_fingerprint(model)
    reaction_ids = {str(reaction_id) for reaction_id in getattr(model, "rxns", ())}
    metabolic_enzymes = tuple(str(item) for item in getattr(metabolic, "enzymes", ()))
    combined_enzymes = tuple(str(item) for item in getattr(combined, "enzymes", ()))
    enzyme_ids = tuple(dict.fromkeys((*metabolic_enzymes, *combined_enzymes)))
    source_refs = _source_refs(metabolic, combined)
    evidence = tuple(external_evidence or ())
    mappings: list[GeneEnzymeReactionMapping] = []

    gene_ids = tuple(
        str(item[0])
        for item in sorted(
            getattr(model, "gene_index", {}).items(),
            key=lambda item: int(item[1]),
        )
    )
    for gene_id in gene_ids:
        plan = plan_gene_overexpression(model, gene_id)
        if not plan.gpr_rules:
            mappings.append(
                _mapping(
                    fingerprint=fingerprint,
                    gene_id=gene_id,
                    enzyme_id="",
                    reaction_id="",
                    gpr_rule="",
                    role=GPRRole.UNRESOLVED,
                    status=OEExecutionStatus.UNRESOLVED,
                    confidence=ConfidenceLevel.LOW,
                    formation_id="",
                    source_refs=source_refs,
                    external_refs=(),
                    missing=("model_gpr_rule", "enzyme_id", "reaction_id"),
                )
            )
            continue
        for rule_payload in plan.gpr_rules:
            reaction_id = str(rule_payload.get("reaction_id") or "")
            role = GPRRole(str(rule_payload.get("gpr_role") or GPRRole.UNRESOLVED.value))
            candidates = tuple(
                enzyme_id
                for enzyme_id in enzyme_ids
                if _reaction_id_for_enzyme(metabolic, enzyme_id) == reaction_id
            )
            external_refs = _matching_external_refs(evidence, gene_id, reaction_id)
            if not candidates:
                mappings.append(
                    _mapping(
                        fingerprint=fingerprint,
                        gene_id=gene_id,
                        enzyme_id="",
                        reaction_id=reaction_id,
                        gpr_rule=str(
                            rule_payload.get("gr_rule")
                            or rule_payload.get("rule")
                            or ""
                        ),
                        role=role,
                        status=_unresolved_status(role),
                        confidence=ConfidenceLevel.LOW,
                        formation_id="",
                        source_refs=source_refs,
                        external_refs=external_refs,
                        missing=("enzyme_id", "formation_or_dilution_reaction_id"),
                    )
                )
                continue
            for enzyme_id in candidates:
                formation_id = _formation_id_for_enzyme(metabolic, enzyme_id)
                status = _execution_status(
                    role=role,
                    candidate_count=len(candidates),
                    formation_available=formation_id in reaction_ids,
                )
                missing = ()
                if formation_id not in reaction_ids:
                    missing = ("formation_or_dilution_reaction_id",)
                mappings.append(
                    _mapping(
                        fingerprint=fingerprint,
                        gene_id=gene_id,
                        enzyme_id=enzyme_id,
                        reaction_id=reaction_id,
                        gpr_rule=str(
                            rule_payload.get("gr_rule")
                            or rule_payload.get("rule")
                            or ""
                        ),
                        role=role,
                        status=status,
                        confidence=(
                            ConfidenceLevel.HIGH
                            if status is OEExecutionStatus.GENE_LEVEL_EXECUTABLE
                            else ConfidenceLevel.MEDIUM
                        ),
                        formation_id=formation_id if formation_id in reaction_ids else "",
                        source_refs=source_refs,
                        external_refs=external_refs,
                        missing=missing,
                    )
                )

    mapped_external_keys = {
        (mapping.gene_id, mapping.reaction_id) for mapping in mappings
    }
    for item in evidence:
        gene_id = str(item.get("gene_id") or "")
        reaction_id = str(item.get("reaction_id") or "")
        if not gene_id or (gene_id, reaction_id) in mapped_external_keys:
            continue
        source_type = _external_source_type(item.get("source_type"))
        mappings.append(
            GeneEnzymeReactionMapping(
                mapping_id=_mapping_id(
                    fingerprint, gene_id, str(item.get("enzyme_id") or ""), reaction_id, "external"
                ),
                model_fingerprint="",
                gene_id=gene_id,
                enzyme_id=str(item.get("enzyme_id") or ""),
                reaction_id=reaction_id,
                gpr_rule=str(item.get("gpr_rule") or ""),
                gpr_role=GPRRole(str(item.get("gpr_role") or GPRRole.UNRESOLVED.value)),
                mapping_source=source_type,
                mapping_confidence=ConfidenceLevel.LOW,
                execution_status=OEExecutionStatus.EXTERNAL_EVIDENCE_ONLY,
                source_ref=str(item.get("source_ref") or "external evidence"),
                source_version=str(item.get("source_version") or ""),
                missing_information=("current_model_gene_enzyme_reaction_mapping",),
                warnings=("External evidence does not create a current-model executable mapping.",),
            )
        )

    catalog = GeneCapacityCatalog(
        model_fingerprint=fingerprint,
        mappings=tuple(sorted(mappings, key=lambda item: item.mapping_id)),
        source_refs=source_refs,
    )
    catalog.validate()
    return catalog


def validate_gene_capacity_catalog(
    catalog: GeneCapacityCatalog,
) -> GeneCapacityValidationResult:
    issues: list[GeneCapacityValidationIssue] = []
    try:
        catalog.validate()
    except OECapacityValidationError as exc:
        issues.append(
            GeneCapacityValidationIssue(
                code="catalog_invalid",
                message=str(exc),
                severity="error",
            )
        )
    for mapping in catalog.mappings:
        for missing in mapping.missing_information:
            issues.append(
                GeneCapacityValidationIssue(
                    code="mapping_missing_information",
                    message=f"{mapping.mapping_id} is missing {missing}.",
                    severity="warning",
                    mapping_id=mapping.mapping_id,
                    field_name=missing,
                )
            )
    result = GeneCapacityValidationResult(issues=tuple(issues))
    result.validate()
    return result


def summarize_gene_capacity_catalog(
    catalog: GeneCapacityCatalog,
) -> GeneCapacityCoverage:
    catalog.validate()
    role_counts = Counter(mapping.gpr_role.value for mapping in catalog.mappings)
    status_counts = Counter(
        mapping.execution_status.value for mapping in catalog.mappings
    )
    coverage = GeneCapacityCoverage(
        total_mappings=len(catalog.mappings),
        gene_count=len({mapping.gene_id for mapping in catalog.mappings}),
        reaction_count=len(
            {mapping.reaction_id for mapping in catalog.mappings if mapping.reaction_id}
        ),
        enzyme_count=len(
            {mapping.enzyme_id for mapping in catalog.mappings if mapping.enzyme_id}
        ),
        by_role=tuple(sorted(role_counts.items())),
        by_status=tuple(sorted(status_counts.items())),
    )
    coverage.validate()
    return coverage


def _mapping(
    *,
    fingerprint: str,
    gene_id: str,
    enzyme_id: str,
    reaction_id: str,
    gpr_rule: str,
    role: GPRRole,
    status: OEExecutionStatus,
    confidence: ConfidenceLevel,
    formation_id: str,
    source_refs: tuple[str, ...],
    external_refs: tuple[str, ...],
    missing: tuple[str, ...],
) -> GeneEnzymeReactionMapping:
    warnings = tuple(
        f"External traceability only: {source_ref}" for source_ref in external_refs
    )
    return GeneEnzymeReactionMapping(
        mapping_id=_mapping_id(fingerprint, gene_id, enzyme_id, reaction_id, role.value),
        model_fingerprint=fingerprint,
        gene_id=gene_id,
        enzyme_id=enzyme_id,
        reaction_id=reaction_id,
        gpr_rule=gpr_rule,
        gpr_role=role,
        enzyme_variable_id=formation_id,
        formation_or_dilution_reaction_id=formation_id,
        mapping_source=EvidenceSourceType.CURRENT_MODEL,
        mapping_confidence=confidence,
        execution_status=status,
        source_ref="; ".join(source_refs),
        warnings=warnings,
        missing_information=missing,
    )


def _execution_status(
    *,
    role: GPRRole,
    candidate_count: int,
    formation_available: bool,
) -> OEExecutionStatus:
    if role is GPRRole.ISOENZYME or candidate_count > 1:
        return OEExecutionStatus.ISOENZYME_AMBIGUOUS
    if role is GPRRole.COMPLEX_SUBUNIT:
        return OEExecutionStatus.COMPLEX_LIMITED
    if role is GPRRole.MIXED:
        return OEExecutionStatus.PARTIAL_MAPPING
    if role is GPRRole.UNRESOLVED or not formation_available:
        return OEExecutionStatus.PARTIAL_MAPPING
    return OEExecutionStatus.GENE_LEVEL_EXECUTABLE


def _unresolved_status(role: GPRRole) -> OEExecutionStatus:
    if role is GPRRole.ISOENZYME:
        return OEExecutionStatus.ISOENZYME_AMBIGUOUS
    if role is GPRRole.COMPLEX_SUBUNIT:
        return OEExecutionStatus.COMPLEX_LIMITED
    return OEExecutionStatus.PARTIAL_MAPPING


def _reaction_id_for_enzyme(metabolic: Any, enzyme_id: str) -> str:
    if hasattr(metabolic, "reaction_id_for_enzyme"):
        return str(metabolic.reaction_id_for_enzyme(enzyme_id))
    return enzyme_id.replace("_complex", "")


def _formation_id_for_enzyme(metabolic: Any, enzyme_id: str) -> str:
    if hasattr(metabolic, "formation_reaction_id_for_enzyme"):
        return str(metabolic.formation_reaction_id_for_enzyme(enzyme_id))
    return f"{enzyme_id}_formation"


def _model_fingerprint(model: Any) -> str:
    payload = {
        "rxns": [str(item) for item in getattr(model, "rxns", ())],
        "genes": [str(item) for item in getattr(model, "genes", ())],
        "gene_index": sorted(
            (str(key), int(value))
            for key, value in getattr(model, "gene_index", {}).items()
        ),
        "rules": [str(item) for item in getattr(model, "rules", ())],
        "gr_rules": [str(item) for item in getattr(model, "gr_rules", ())],
    }
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _mapping_id(
    fingerprint: str,
    gene_id: str,
    enzyme_id: str,
    reaction_id: str,
    role: str,
) -> str:
    raw = "::".join((fingerprint, gene_id, enzyme_id, reaction_id, role))
    return "oe-map-" + hashlib.sha256(raw.encode("utf-8")).hexdigest()[:16]


def _source_refs(metabolic: Any, combined: Any) -> tuple[str, ...]:
    refs: list[str] = []
    metabolic_source = getattr(metabolic, "source_file", None)
    if metabolic_source:
        refs.append(str(Path(metabolic_source)))
    refs.extend(str(Path(item)) for item in getattr(combined, "source_files", ()) if item)
    return tuple(dict.fromkeys(refs))


def _matching_external_refs(
    evidence: Sequence[Mapping[str, Any]],
    gene_id: str,
    reaction_id: str,
) -> tuple[str, ...]:
    return tuple(
        str(item.get("source_ref") or "external evidence")
        for item in evidence
        if str(item.get("gene_id") or "") == gene_id
        and str(item.get("reaction_id") or "") == reaction_id
    )


def _external_source_type(value: Any) -> EvidenceSourceType:
    try:
        source_type = EvidenceSourceType(str(value))
    except ValueError:
        return EvidenceSourceType.EXTERNAL_PICHIA_MODEL
    if source_type in {
        EvidenceSourceType.CURRENT_MODEL,
        EvidenceSourceType.LOCAL_ENZYME_DATA,
        EvidenceSourceType.REVIEWED_PICHIA_MAPPING,
    }:
        return EvidenceSourceType.EXTERNAL_PICHIA_MODEL
    return source_type


__all__ = [
    "build_gene_enzyme_reaction_catalog",
    "summarize_gene_capacity_catalog",
    "validate_gene_capacity_catalog",
]
