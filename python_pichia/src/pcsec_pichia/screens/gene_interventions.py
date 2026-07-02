from __future__ import annotations

import re
from dataclasses import asdict, dataclass
from typing import Any

from pcsec_pichia.screens.candidate_resolution import reactions_for_gene
from pcsec_pichia.screens.gene_perturbation_map import build_reaction_perturbation_mapping
from pcsec_pichia.services.gene_evidence import (
    GeneExternalEvidence,
    evidence_for_gene,
    recommendation_tier_for_candidate,
    resolve_gene_identifier,
)


RULE_TOKEN_PATTERN = re.compile(r"x\((\d+)\)")


@dataclass(frozen=True)
class GeneInterventionPlan:
    gene_id: str
    intervention_type: str
    resolved: bool
    affected_reactions: tuple[str, ...]
    inactive_reactions: tuple[str, ...]
    executable_reactions: tuple[str, ...]
    explain_only_reactions: tuple[str, ...]
    gpr_rules: tuple[dict[str, str], ...]
    gpr_role: str
    capacity_effect: str
    simulation_basis: str
    mapping_confidence: str
    warnings: tuple[str, ...]
    ko_support_status: str = ""
    oe_support_status: str = ""
    support_reason: str = ""
    missing_information: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, object]:
        return asdict(self)

    def candidate_fields(self) -> dict[str, object]:
        return {
            "affected_reactions": list(self.affected_reactions),
            "inactive_reactions": list(self.inactive_reactions),
            "inactive_reactions_preview": list(self.inactive_reactions[:10]),
            "inactive_reaction_count": len(self.inactive_reactions),
            "gpr_rules": list(self.gpr_rules),
            "gpr_role": self.gpr_role,
            "capacity_effect": self.capacity_effect,
            "simulation_basis": self.simulation_basis,
            "mapping_confidence": self.mapping_confidence,
            "ko_support_status": self.ko_support_status,
            "oe_support_status": self.oe_support_status,
            "support_reason": self.support_reason,
            "missing_information": list(self.missing_information),
            "warnings": list(self.warnings),
        }


@dataclass(frozen=True)
class GeneCapabilityProfile:
    gene_id: str
    canonical_gene_id: str
    aliases: tuple[str, ...]
    resolved: bool
    affected_reactions: tuple[str, ...]
    inactive_reactions_if_ko: tuple[str, ...]
    oe_executable_reactions: tuple[str, ...]
    oe_explain_only_reactions: tuple[str, ...]
    gpr_rules: tuple[dict[str, str], ...]
    gpr_role: str
    ko_support_status: str
    oe_support_status: str
    support_reason: str
    missing_information: tuple[str, ...]
    confidence: str
    database_annotation_sources: tuple[str, ...]
    database_annotation_confidence: str
    model_gpr_executable: bool
    oe_reaction_proxy: bool
    phenotype_evidence: dict[str, dict[str, object]]
    recommendation_tier: dict[str, str]
    recommendation_tier_reason: dict[str, str]

    def to_dict(self) -> dict[str, object]:
        payload = asdict(self)
        for key in (
            "aliases",
            "affected_reactions",
            "inactive_reactions_if_ko",
            "oe_executable_reactions",
            "oe_explain_only_reactions",
            "gpr_rules",
            "missing_information",
            "database_annotation_sources",
        ):
            payload[key] = list(payload[key])
        return payload


def plan_gene_knockout(model: Any, gene_id: str) -> GeneInterventionPlan:
    entries = _reaction_entries_for_gene(model, gene_id)
    if not _gene_exists(model, gene_id):
        return _unresolved_plan(gene_id, "KO", "gene_not_found_or_no_gpr_reactions")
    if not entries:
        return GeneInterventionPlan(
            gene_id=gene_id,
            intervention_type="KO",
            resolved=True,
            affected_reactions=(),
            inactive_reactions=(),
            executable_reactions=(),
            explain_only_reactions=(),
            gpr_rules=(),
            gpr_role="unresolved",
            capacity_effect="no_reaction_disabled",
            simulation_basis="gpr_gene_deletion",
            mapping_confidence="low",
            ko_support_status="ko_no_gpr_effect",
            support_reason="Gene exists in the model, but no reaction GPR currently references it.",
            missing_information=("model_gpr_rule",),
            warnings=("Gene exists in the model, but no reaction GPR currently references it.",),
        )

    inactive_entries: list[str] = []
    unevaluable_entries: list[str] = []
    gene_number = _gene_number(model, gene_id)
    for entry in entries:
        evaluation = _evaluate_rule_after_ko(entry["rule"], gene_number)
        if evaluation is False:
            inactive_entries.append(entry["reaction_id"])
        elif evaluation is None:
            unevaluable_entries.append(entry["reaction_id"])
    inactive = tuple(inactive_entries)
    role = _combined_gpr_role(tuple(_rule_role(entry["rule"], entry["gr_rule"]) for entry in entries))
    warnings: list[str] = []
    if not inactive:
        warnings.append("KO gene is present in GPR rules, but no reaction becomes inactive after AND/OR evaluation.")
    if unevaluable_entries:
        warnings.append("Some gene-associated reactions lack evaluable model rule tokens; they are not disabled by KO planning.")
    missing_information = ("model_rule_token_mapping",) if unevaluable_entries else ()
    return GeneInterventionPlan(
        gene_id=gene_id,
        intervention_type="KO",
        resolved=True,
        affected_reactions=tuple(entry["reaction_id"] for entry in entries),
        inactive_reactions=inactive,
        executable_reactions=inactive,
        explain_only_reactions=(),
        gpr_rules=_gpr_rule_payload(entries),
        gpr_role=role,
        capacity_effect="disables_reactions" if inactive else "no_reaction_disabled",
        simulation_basis="gpr_gene_deletion",
        mapping_confidence="high" if inactive else "medium",
        ko_support_status="ko_runnable_gpr_gene_deletion" if inactive else "ko_no_reaction_disabled",
        support_reason=(
            "GPR deletion disables one or more model reactions."
            if inactive
            else (
                "Gene is referenced by textual GPR, but numeric rule tokens are missing for one or more reactions."
                if unevaluable_entries
                else "Gene is in GPR rules, but AND/OR evaluation leaves all associated reactions active."
            )
        ),
        missing_information=missing_information,
        warnings=tuple(warnings),
    )


def plan_gene_overexpression(
    model: Any,
    gene_id: str,
    complex_subunits: dict[str, list[dict[str, object]]] | None = None,
) -> GeneInterventionPlan:
    entries = _reaction_entries_for_gene(model, gene_id)
    if not _gene_exists(model, gene_id):
        return _unresolved_plan(gene_id, "OE_gene_proxy", "gene_not_found_or_no_gpr_reactions")
    if not entries:
        return GeneInterventionPlan(
            gene_id=gene_id,
            intervention_type="OE_gene_proxy",
            resolved=True,
            affected_reactions=(),
            inactive_reactions=(),
            executable_reactions=(),
            explain_only_reactions=(),
            gpr_rules=(),
            gpr_role="unresolved",
            capacity_effect="no_gpr_effect",
            simulation_basis="explain_only",
            mapping_confidence="low",
            oe_support_status="oe_no_gpr_effect",
            support_reason="Gene exists in the model, but no reaction GPR currently references it.",
            missing_information=("model_gpr_rule",),
            warnings=("Gene exists in the model, but no reaction GPR currently references it.",),
        )

    executable: list[str] = []
    explain_only: list[str] = []
    roles: list[str] = []
    warnings: list[str] = []
    for entry in entries:
        reaction_id = entry["reaction_id"]
        role = _rule_role(entry["rule"], entry["gr_rule"])
        mapping = build_reaction_perturbation_mapping(reaction_id, complex_subunits or {})
        if mapping.mapping_level == "complex_subunit" or role in {"complex_subunit", "mixed"}:
            explain_only.append(reaction_id)
            roles.append("complex_subunit" if mapping.mapping_level == "complex_subunit" else role)
            continue
        executable.append(reaction_id)
        roles.append(role)

    combined_role = _combined_gpr_role(tuple(roles))
    if explain_only:
        warnings.append(
            "Some OE gene targets are complex-subunit or mixed GPR reactions; single-gene OE is explain-only for those rows."
        )
    if executable and explain_only:
        capacity_effect = "partial_reaction_capacity_proxy"
        confidence = "medium"
    elif executable:
        capacity_effect = "reaction_capacity_proxy"
        confidence = "high" if combined_role == "single_gene" else "medium"
    else:
        capacity_effect = "complex_subunit_limited" if combined_role == "complex_subunit" else "manual_review_required"
        confidence = "low"
    return GeneInterventionPlan(
        gene_id=gene_id,
        intervention_type="OE_gene_proxy",
        resolved=True,
        affected_reactions=tuple(entry["reaction_id"] for entry in entries),
        inactive_reactions=(),
        executable_reactions=tuple(executable),
        explain_only_reactions=tuple(explain_only),
        gpr_rules=_gpr_rule_payload(entries),
        gpr_role=combined_role,
        capacity_effect=capacity_effect,
        simulation_basis="reaction_level_capacity_proxy" if executable else "explain_only",
        mapping_confidence=confidence,
        oe_support_status=_oe_support_status(executable, explain_only, capacity_effect),
        support_reason=_oe_support_reason(executable, explain_only, capacity_effect),
        missing_information=_oe_missing_information(executable, explain_only, capacity_effect),
        warnings=tuple(warnings),
    )


def build_gene_capability_profile(
    model: Any,
    gene_id: str,
    complex_subunits: dict[str, list[dict[str, object]]] | None = None,
    aliases: tuple[str, ...] = (),
    evidence_by_gene: dict[str, GeneExternalEvidence] | None = None,
    target_protein_context: str | None = None,
) -> GeneCapabilityProfile:
    input_gene_id = str(gene_id or "").strip()
    evidence_records = evidence_by_gene or {}
    external_evidence = evidence_for_gene(input_gene_id, evidence_records)
    canonical_gene_id = resolve_gene_identifier(input_gene_id, evidence_records, model=model)
    profile_aliases = tuple(aliases or (external_evidence.aliases if external_evidence else ()))
    if not canonical_gene_id or not _gene_exists(model, canonical_gene_id):
        evidence_fields = _capability_evidence_fields(
            canonical_gene_id or input_gene_id,
            profile_aliases,
            external_evidence,
            target_protein_context,
            resolved=False,
        )
        return GeneCapabilityProfile(
            gene_id=input_gene_id,
            canonical_gene_id=canonical_gene_id,
            aliases=profile_aliases,
            resolved=False,
            affected_reactions=(),
            inactive_reactions_if_ko=(),
            oe_executable_reactions=(),
            oe_explain_only_reactions=(),
            gpr_rules=(),
            gpr_role="unresolved",
            ko_support_status="unresolved_gene",
            oe_support_status="unresolved_gene",
            support_reason="Gene ID is not present in the current model gene_index.",
            missing_information=("model_gene_id_or_alias_mapping",),
            confidence="unresolved",
            **evidence_fields,
        )

    entries = _reaction_entries_for_gene(model, canonical_gene_id)
    if not entries:
        evidence_fields = _capability_evidence_fields(
            canonical_gene_id,
            profile_aliases,
            external_evidence,
            target_protein_context,
            resolved=True,
        )
        return GeneCapabilityProfile(
            gene_id=input_gene_id,
            canonical_gene_id=canonical_gene_id,
            aliases=profile_aliases,
            resolved=True,
            affected_reactions=(),
            inactive_reactions_if_ko=(),
            oe_executable_reactions=(),
            oe_explain_only_reactions=(),
            gpr_rules=(),
            gpr_role="unresolved",
            ko_support_status="ko_no_gpr_effect",
            oe_support_status="oe_no_gpr_effect",
            support_reason="Gene exists in the model, but no reaction GPR currently references it.",
            missing_information=("model_gpr_rule",),
            confidence="low",
            **evidence_fields,
        )

    ko_plan = plan_gene_knockout(model, canonical_gene_id)
    oe_plan = plan_gene_overexpression(model, canonical_gene_id, complex_subunits=complex_subunits)
    missing = tuple(dict.fromkeys((*ko_plan.missing_information, *oe_plan.missing_information)))
    evidence_fields = _capability_evidence_fields(
        canonical_gene_id,
        profile_aliases,
        external_evidence,
        target_protein_context,
        resolved=True,
        model_gpr_executable=ko_plan.ko_support_status == "ko_runnable_gpr_gene_deletion",
        oe_reaction_proxy=bool(oe_plan.executable_reactions),
    )
    return GeneCapabilityProfile(
        gene_id=input_gene_id,
        canonical_gene_id=canonical_gene_id,
        aliases=profile_aliases,
        resolved=True,
        affected_reactions=ko_plan.affected_reactions,
        inactive_reactions_if_ko=ko_plan.inactive_reactions,
        oe_executable_reactions=oe_plan.executable_reactions,
        oe_explain_only_reactions=oe_plan.explain_only_reactions,
        gpr_rules=ko_plan.gpr_rules,
        gpr_role=_combined_gpr_role((ko_plan.gpr_role, oe_plan.gpr_role)),
        ko_support_status=ko_plan.ko_support_status or "ko_no_gpr_effect",
        oe_support_status=oe_plan.oe_support_status or "oe_explain_only_no_capacity_model",
        support_reason=_combine_support_reason(ko_plan.support_reason, oe_plan.support_reason),
        missing_information=missing,
        confidence=_combine_confidence(ko_plan.mapping_confidence, oe_plan.mapping_confidence),
        **evidence_fields,
    )


def build_all_gene_capability_catalog(
    model: Any,
    complex_subunits: dict[str, list[dict[str, object]]] | None = None,
) -> tuple[GeneCapabilityProfile, ...]:
    return tuple(
        build_gene_capability_profile(model, gene_id, complex_subunits=complex_subunits)
        for gene_id in _model_gene_ids(model)
    )


def _capability_evidence_fields(
    gene_id: str,
    aliases: tuple[str, ...],
    external_evidence: GeneExternalEvidence | None,
    target_protein_context: str | None,
    *,
    resolved: bool,
    model_gpr_executable: bool = False,
    oe_reaction_proxy: bool = False,
) -> dict[str, object]:
    database_sources = tuple(external_evidence.evidence_sources) if external_evidence else ()
    database_confidence = external_evidence.evidence_confidence if external_evidence else ""
    ko_tier, ko_reason, ko_evidence = recommendation_tier_for_candidate(
        gene_id=gene_id,
        intervention_type="KO",
        target_protein_context=target_protein_context,
        model_gpr_executable=model_gpr_executable,
        resolved=resolved,
        database_annotation_available=bool(external_evidence),
        aliases=aliases,
    )
    oe_tier, oe_reason, oe_evidence = recommendation_tier_for_candidate(
        gene_id=gene_id,
        intervention_type="OE",
        target_protein_context=target_protein_context,
        oe_reaction_proxy=oe_reaction_proxy,
        resolved=resolved,
        database_annotation_available=bool(external_evidence),
        aliases=aliases,
    )
    phenotype_payload = {
        key: evidence.to_dict()
        for key, evidence in (("KO", ko_evidence), ("OE", oe_evidence))
        if evidence is not None
    }
    return {
        "database_annotation_sources": database_sources,
        "database_annotation_confidence": database_confidence,
        "model_gpr_executable": model_gpr_executable,
        "oe_reaction_proxy": oe_reaction_proxy,
        "phenotype_evidence": phenotype_payload,
        "recommendation_tier": {"KO": ko_tier, "OE": oe_tier},
        "recommendation_tier_reason": {"KO": ko_reason, "OE": oe_reason},
    }


def _reaction_entries_for_gene(model: Any, gene_id: str) -> tuple[dict[str, str], ...]:
    reaction_ids = set(reactions_for_gene(model, gene_id))
    entries: list[dict[str, str]] = []
    for reaction_id, rule, gr_rule in zip(
        getattr(model, "rxns", []),
        getattr(model, "rules", []),
        getattr(model, "gr_rules", []),
    ):
        if str(reaction_id) not in reaction_ids:
            continue
        entries.append(
            {
                "reaction_id": str(reaction_id),
                "rule": str(rule or ""),
                "gr_rule": str(gr_rule or ""),
            }
        )
    return tuple(entries)


def _model_gene_ids(model: Any) -> tuple[str, ...]:
    genes = tuple(str(gene_id) for gene_id in getattr(model, "genes", ()) if str(gene_id))
    if genes:
        return genes
    gene_index = getattr(model, "gene_index", {})
    return tuple(str(item[0]) for item in sorted(gene_index.items(), key=lambda item: int(item[1])))


def _gene_exists(model: Any, gene_id: str) -> bool:
    return gene_id in getattr(model, "gene_index", {})


def _gene_number(model: Any, gene_id: str) -> int:
    return int(getattr(model, "gene_index", {})[gene_id]) + 1


def _gpr_rule_payload(entries: tuple[dict[str, str], ...]) -> tuple[dict[str, str], ...]:
    return tuple(
        {
            "reaction_id": entry["reaction_id"],
            "rule": entry["rule"],
            "gr_rule": entry["gr_rule"],
            "gpr_role": _rule_role(entry["rule"], entry["gr_rule"]),
        }
        for entry in entries
    )


def _rule_role(rule: str, gr_rule: str = "") -> str:
    text = f"{rule} {gr_rule}".strip().lower()
    if not text or text == "[]":
        return "unresolved"
    has_and = "&" in text or " and " in f" {text} "
    has_or = "|" in text or " or " in f" {text} "
    if has_and and has_or:
        return "mixed"
    if has_and:
        return "complex_subunit"
    if has_or:
        return "isoenzyme"
    return "single_gene"


def _combined_gpr_role(roles: tuple[str, ...]) -> str:
    role_set = {role for role in roles if role and role != "unresolved"}
    if not role_set:
        return "unresolved"
    if "mixed" in role_set or len(role_set) > 1:
        return "mixed"
    return next(iter(role_set))


def _evaluate_rule_after_ko(rule: str, knocked_gene_number: int) -> bool | None:
    expression = str(rule or "").replace("&", " and ").replace("|", " or ")
    if not RULE_TOKEN_PATTERN.search(expression):
        return None

    def repl(match: re.Match[str]) -> str:
        gene_number = int(match.group(1))
        return "False" if gene_number == knocked_gene_number else "True"

    expression = RULE_TOKEN_PATTERN.sub(repl, expression)
    if "x(" in expression:
        return None
    try:
        return bool(eval(expression, {"__builtins__": {}}, {}))
    except Exception:
        return None


def _unresolved_plan(gene_id: str, intervention_type: str, warning: str) -> GeneInterventionPlan:
    ko_support_status = "unresolved_gene" if intervention_type == "KO" else ""
    oe_support_status = "unresolved_gene" if intervention_type == "OE_gene_proxy" else ""
    return GeneInterventionPlan(
        gene_id=gene_id,
        intervention_type=intervention_type,
        resolved=False,
        affected_reactions=(),
        inactive_reactions=(),
        executable_reactions=(),
        explain_only_reactions=(),
        gpr_rules=(),
        gpr_role="unresolved",
        capacity_effect="unresolved",
        simulation_basis="unresolved",
        mapping_confidence="unresolved",
        ko_support_status=ko_support_status,
        oe_support_status=oe_support_status,
        support_reason="Gene ID is not present in the current model gene_index or has no model GPR reactions.",
        missing_information=("model_gene_id_or_alias_mapping",),
        warnings=(warning,),
    )


def _oe_support_status(
    executable: list[str],
    explain_only: list[str],
    capacity_effect: str,
) -> str:
    if executable:
        return "oe_runnable_reaction_proxy"
    if capacity_effect == "no_gpr_effect":
        return "oe_no_gpr_effect"
    if explain_only and capacity_effect == "complex_subunit_limited":
        return "oe_explain_only_complex_subunit"
    if explain_only:
        return "oe_explain_only_no_capacity_model"
    return "oe_no_gpr_effect"


def _oe_support_reason(
    executable: list[str],
    explain_only: list[str],
    capacity_effect: str,
) -> str:
    if executable and explain_only:
        return "Some reactions can run as reaction-level OE proxy; complex-subunit or mixed GPR reactions remain explain-only."
    if executable:
        return "Gene-associated reactions can run as reaction-level capacity proxy."
    if explain_only and capacity_effect == "complex_subunit_limited":
        return "Single-gene OE of a complex subunit is not treated as a reliable complex capacity increase."
    if explain_only:
        return "Gene is model-associated, but current evidence is insufficient for a quantitative OE proxy."
    return "Gene has no model GPR effect for OE."


def _oe_missing_information(
    executable: list[str],
    explain_only: list[str],
    capacity_effect: str,
) -> tuple[str, ...]:
    missing: list[str] = []
    if explain_only:
        missing.append("gene_expression_to_complex_capacity_model")
    if capacity_effect in {"complex_subunit_limited", "manual_review_required"}:
        missing.append("complex_assembly_or_limiting_subunit_evidence")
    if not executable and not explain_only:
        missing.append("model_gpr_rule")
    return tuple(dict.fromkeys(missing))


def _combine_support_reason(*items: str) -> str:
    return " ".join(item.strip() for item in items if item and item.strip())


def _combine_confidence(*items: str) -> str:
    values = tuple(item for item in items if item)
    if "high" in values:
        return "high"
    if "medium" in values:
        return "medium"
    if "low" in values:
        return "low"
    return "unresolved"


__all__ = [
    "GeneCapabilityProfile",
    "GeneInterventionPlan",
    "build_all_gene_capability_catalog",
    "build_gene_capability_profile",
    "plan_gene_knockout",
    "plan_gene_overexpression",
]
