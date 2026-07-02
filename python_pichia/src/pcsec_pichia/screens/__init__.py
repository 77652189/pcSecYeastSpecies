from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from pcsec_pichia.screens._prototype_adapter import (
    AminoAcidStoichiometry,
    CobraModel,
    CombinedEnzymeData,
    MetabolicEnzymeData,
    SecretoryEnzymeData,
    TargetSpec,
    build_supported_target_model,
    build_target_enzymedata,
    classify_candidate_effect,
    classify_secretory_process,
    default_ko_genes,
    default_oe_reactions,
    run_ko_screen,
    run_oe_screen,
    run_pcsec_growth_tradeoff,
    run_pcsec_ko_screen,
    run_pcsec_oe_screen,
    run_pcsec_reaction_ko_screen,
    solve_pcsec_maximize,
)
from pcsec_pichia.screens.candidate_resolution import (
    reactions_for_gene,
    resolve_oe_gene_reactions,
    split_existing_genes,
    split_existing_reactions,
)
from pcsec_pichia.screens.gene_perturbation_map import (
    GenePerturbationMapping,
    GenePerturbationMapResult,
    GeneReactionMapping,
    build_gene_perturbation_map,
    build_reaction_perturbation_mapping,
)
from pcsec_pichia.screens.gene_interventions import (
    GeneCapabilityProfile,
    GeneInterventionPlan,
    build_all_gene_capability_catalog,
    build_gene_capability_profile,
    plan_gene_knockout,
    plan_gene_overexpression,
)
from pcsec_pichia.screens.planning import ScreenPlanResult, build_screen_plan


@dataclass(frozen=True)
class ScreenResult:
    target_id: str
    screen_type: str
    success: bool
    candidate_count: int
    rows: tuple[dict[str, Any], ...]
    constraint_counts: dict[str, int]
    baseline_objective_value: float | None
    result_status: str
    matlab_alignment_status: str


def run_knockout_screen(
    model: CobraModel,
    target: TargetSpec,
    amino_acids: AminoAcidStoichiometry,
    metabolic: MetabolicEnzymeData,
    secretory: SecretoryEnzymeData,
    combined: CombinedEnzymeData,
    genes: list[str],
    growth_rate: float = 0.10,
    write_ribosome_translation_constraint: bool = False,
    write_misfolding_constraints: bool = False,
) -> ScreenResult:
    if not genes:
        return _empty_unsolved_screen_result(target.target_id, "knockout")

    prepared = _prepare_screen_inputs(
        model,
        target,
        amino_acids,
        metabolic,
        secretory,
        combined,
        growth_rate,
        write_ribosome_translation_constraint=write_ribosome_translation_constraint,
        write_misfolding_constraints=write_misfolding_constraints,
    )
    if not prepared["baseline_success"]:
        return _empty_screen_result(target.target_id, "knockout", prepared)

    plans = {gene_id: plan_gene_knockout(prepared["fixed_model"], gene_id) for gene_id in genes}
    raw_by_gene = {
        gene_id: _solve_gene_knockout_plan(
            prepared,
            plans[gene_id],
            metabolic=metabolic,
            growth_rate=growth_rate,
            write_ribosome_translation_constraint=write_ribosome_translation_constraint,
            write_misfolding_constraints=write_misfolding_constraints,
        )
        for gene_id in genes
        if plans[gene_id].inactive_reactions
    }
    rows = tuple(
        _normalize_screen_row(
            {
                **(
                    raw_by_gene.get(gene_id)
                    or (
                        _unresolved_gene_knockout_row(plans[gene_id], prepared["baseline"].objective_value)
                        if not plans[gene_id].resolved
                        else _no_effect_gene_knockout_row(plans[gene_id], prepared["baseline"].objective_value)
                    )
                ),
                **_gene_plan_fields(plans[gene_id]),
            },
            target_id=target.target_id,
            screen_type="knockout",
            intervention_type="KO",
            baseline_objective_value=prepared["baseline"].objective_value,
            complex_subunits=prepared["secretory"].complex_subunits,
            input_gene_id=gene_id,
        )
        for gene_id in genes
    )
    return _screen_result(target.target_id, "knockout", rows, prepared)


def _solve_gene_knockout_plan(
    prepared: dict[str, Any],
    plan: GeneInterventionPlan,
    metabolic: MetabolicEnzymeData,
    growth_rate: float,
    write_ribosome_translation_constraint: bool,
    write_misfolding_constraints: bool,
) -> dict[str, Any]:
    changes = {reaction_id: (0.0, 0.0) for reaction_id in plan.inactive_reactions}
    solved, counts = solve_pcsec_maximize(
        prepared["fixed_model"].with_bounds(changes),
        prepared["exchange_reaction_id"],
        metabolic=metabolic,
        secretory=prepared["secretory"],
        combined=prepared["combined"],
        mu=growth_rate,
        key_reactions=("BIOMASS", "Ex_glc_D", prepared["exchange_reaction_id"]),
        write_ribosome_translation_constraint=write_ribosome_translation_constraint,
        write_misfolding_constraints=write_misfolding_constraints,
    )
    baseline = prepared["baseline"]
    return {
        "gene": plan.gene_id,
        "inactive_reaction_count": len(plan.inactive_reactions),
        "inactive_reactions_preview": list(plan.inactive_reactions[:10]),
        "inactive_reactions": list(plan.inactive_reactions),
        "status": solved.status,
        "success": solved.success,
        "objective_value": solved.objective_value,
        "delta_vs_baseline": (
            solved.objective_value - baseline.objective_value
            if solved.success and baseline.objective_value is not None and solved.objective_value is not None
            else None
        ),
        "constraint_counts": counts,
    }


def run_overexpression_screen(
    model: CobraModel,
    target: TargetSpec,
    amino_acids: AminoAcidStoichiometry,
    metabolic: MetabolicEnzymeData,
    secretory: SecretoryEnzymeData,
    combined: CombinedEnzymeData,
    reactions: list[str],
    growth_rate: float = 0.10,
    factor: float = 2.0,
    intervention_type: str = "OE_reaction",
    input_gene_ids_by_reaction: dict[str, str] | None = None,
    gene_intervention_plans_by_gene: dict[str, GeneInterventionPlan] | None = None,
    write_ribosome_translation_constraint: bool = False,
    write_misfolding_constraints: bool = False,
) -> ScreenResult:
    if not reactions:
        return _empty_unsolved_screen_result(target.target_id, "overexpression")

    prepared = _prepare_screen_inputs(
        model,
        target,
        amino_acids,
        metabolic,
        secretory,
        combined,
        growth_rate,
        write_ribosome_translation_constraint=write_ribosome_translation_constraint,
        write_misfolding_constraints=write_misfolding_constraints,
    )
    if not prepared["baseline_success"]:
        return _empty_screen_result(target.target_id, "overexpression", prepared)

    raw_rows = run_pcsec_oe_screen(
        prepared["fixed_model"],
        prepared["baseline"],
        reactions,
        prepared["exchange_reaction_id"],
        metabolic=metabolic,
        secretory=prepared["secretory"],
        combined=prepared["combined"],
        mu=growth_rate,
        factor=factor,
        write_ribosome_translation_constraint=write_ribosome_translation_constraint,
        write_misfolding_constraints=write_misfolding_constraints,
    )
    rows = tuple(
        _normalize_screen_row(
            {**row, **_gene_plan_fields(_plan_for_reaction(row, input_gene_ids_by_reaction, gene_intervention_plans_by_gene))},
            target_id=target.target_id,
            screen_type="overexpression",
            intervention_type=intervention_type,
            baseline_objective_value=prepared["baseline"].objective_value,
            complex_subunits=prepared["secretory"].complex_subunits,
            input_gene_id=(input_gene_ids_by_reaction or {}).get(str(row.get("reaction"))),
        )
        for row in raw_rows
    )
    return _screen_result(target.target_id, "overexpression", rows, prepared)


def explain_only_gene_overexpression_rows(
    target_id: str,
    plans: tuple[GeneInterventionPlan, ...],
    baseline_objective_value: float | None,
    complex_subunits: dict[str, list[dict[str, object]]] | None = None,
) -> tuple[dict[str, Any], ...]:
    rows: list[dict[str, Any]] = []
    for plan in plans:
        reaction_id = (plan.explain_only_reactions or plan.affected_reactions or (None,))[0]
        raw = {
            "gene": plan.gene_id,
            "reaction": reaction_id,
            "success": False,
            "status": _explain_only_oe_status(plan),
            "objective_value": None,
            "delta_vs_baseline": None,
            **_gene_plan_fields(plan),
        }
        rows.append(
            _normalize_screen_row(
                raw,
                target_id=target_id,
                screen_type="overexpression",
                intervention_type="OE_gene_proxy",
                baseline_objective_value=baseline_objective_value,
                complex_subunits=complex_subunits,
                input_gene_id=plan.gene_id,
            )
        )
    return tuple(rows)


def run_reaction_knockout_screen(
    model: CobraModel,
    target: TargetSpec,
    amino_acids: AminoAcidStoichiometry,
    metabolic: MetabolicEnzymeData,
    secretory: SecretoryEnzymeData,
    combined: CombinedEnzymeData,
    reactions: list[str],
    growth_rate: float = 0.10,
    write_ribosome_translation_constraint: bool = False,
    write_misfolding_constraints: bool = False,
) -> ScreenResult:
    if not reactions:
        return _empty_unsolved_screen_result(target.target_id, "knockout")

    prepared = _prepare_screen_inputs(
        model,
        target,
        amino_acids,
        metabolic,
        secretory,
        combined,
        growth_rate,
        write_ribosome_translation_constraint=write_ribosome_translation_constraint,
        write_misfolding_constraints=write_misfolding_constraints,
    )
    if not prepared["baseline_success"]:
        return _empty_screen_result(target.target_id, "knockout", prepared)

    raw_rows = run_pcsec_reaction_ko_screen(
        prepared["fixed_model"],
        prepared["baseline"],
        reactions,
        prepared["exchange_reaction_id"],
        metabolic=metabolic,
        secretory=prepared["secretory"],
        combined=prepared["combined"],
        mu=growth_rate,
        write_ribosome_translation_constraint=write_ribosome_translation_constraint,
        write_misfolding_constraints=write_misfolding_constraints,
    )
    rows = tuple(
        _normalize_screen_row(
            row,
            target_id=target.target_id,
            screen_type="knockout",
            intervention_type="KO_reaction",
            baseline_objective_value=prepared["baseline"].objective_value,
            complex_subunits=prepared["secretory"].complex_subunits,
        )
        for row in raw_rows
    )
    return _screen_result(target.target_id, "knockout", rows, prepared)


def summarize_screen_result(result: ScreenResult) -> dict[str, Any]:
    return {
        "target_id": result.target_id,
        "screen_type": result.screen_type,
        "success": result.success,
        "candidate_count": result.candidate_count,
        "rows": result.rows,
        "constraint_counts": result.constraint_counts,
        "baseline_objective_value": result.baseline_objective_value,
        "result_status": result.result_status,
        "matlab_alignment_status": result.matlab_alignment_status,
    }


def _prepare_screen_inputs(
    model: CobraModel,
    target: TargetSpec,
    amino_acids: AminoAcidStoichiometry,
    metabolic: MetabolicEnzymeData,
    secretory: SecretoryEnzymeData,
    combined: CombinedEnzymeData,
    growth_rate: float,
    write_ribosome_translation_constraint: bool = False,
    write_misfolding_constraints: bool = False,
) -> dict[str, Any]:
    build = build_supported_target_model(model, target, amino_acids)
    if not build.supported or build.model is None or build.exchange_reaction_id is None:
        return {
            "baseline_success": False,
            "baseline": None,
            "constraint_counts": {},
            "fixed_model": None,
            "secretory": None,
            "combined": None,
            "exchange_reaction_id": build.exchange_reaction_id,
        }

    target_enzymedata = build_target_enzymedata(target, build.model, secretory)
    target_secretory = secretory.with_reaction_coefficients(target_enzymedata.reaction_coefficients)
    target_combined = combined.with_target(target_enzymedata)
    fixed_model = build.model.with_bounds({"BIOMASS": (growth_rate, growth_rate)})
    baseline, counts = solve_pcsec_maximize(
        fixed_model,
        build.exchange_reaction_id,
        metabolic=metabolic,
        secretory=target_secretory,
        combined=target_combined,
        mu=growth_rate,
        write_ribosome_translation_constraint=write_ribosome_translation_constraint,
        write_misfolding_constraints=write_misfolding_constraints,
    )
    return {
        "baseline_success": baseline.success,
        "baseline": baseline,
        "constraint_counts": {str(key): int(value) for key, value in counts.items()},
        "fixed_model": fixed_model,
        "secretory": target_secretory,
        "combined": target_combined,
        "exchange_reaction_id": build.exchange_reaction_id,
    }


def _normalize_screen_row(
    row: dict[str, Any],
    target_id: str,
    screen_type: str,
    intervention_type: str,
    baseline_objective_value: float | None,
    complex_subunits: dict[str, list[dict[str, object]]] | None,
    input_gene_id: str | None = None,
    resolved_reaction_id: str | None = None,
) -> dict[str, Any]:
    gene_id = row.get("gene")
    reaction_id = row.get("reaction")
    if reaction_id is None:
        preview = row.get("inactive_reactions_preview") or []
        reaction_id = preview[0] if preview else None
    complex_id = str(reaction_id).replace("_formation", "") if reaction_id else ""
    subunits = (complex_subunits or {}).get(complex_id, [])
    objective_value = row.get("objective_value")
    delta = row.get("delta_vs_baseline")
    normalized_delta = float(delta) if delta is not None else None
    normalized_baseline = float(baseline_objective_value) if baseline_objective_value is not None else None
    relative_delta = _relative_delta(normalized_delta, normalized_baseline)
    resolved_reaction = resolved_reaction_id or (str(reaction_id) if reaction_id is not None else None)
    process_code = classify_secretory_process(resolved_reaction)
    effect_code = classify_candidate_effect(bool(row.get("success")), relative_delta)
    solver_status_label = _solver_status_label(row.get("status"), bool(row.get("success")))
    mapping = build_reaction_perturbation_mapping(resolved_reaction, complex_subunits)
    default_basis = _default_simulation_basis(intervention_type)
    return {
        **row,
        "target_id": target_id,
        "screen_type": screen_type,
        "intervention_type": intervention_type,
        "candidate_id": str(gene_id or reaction_id or ""),
        "gene_id": str(gene_id) if gene_id is not None else None,
        "canonical_gene_id": str(gene_id) if gene_id is not None else None,
        "reaction_id": str(reaction_id) if reaction_id is not None else None,
        "input_gene_id": input_gene_id,
        "resolved_reaction_id": resolved_reaction,
        "objective_value": float(objective_value) if objective_value is not None else None,
        "baseline_objective_value": baseline_objective_value,
        "delta_objective": normalized_delta,
        "effect_label": _effect_label(effect_code, row.get("status")),
        "solver_status_label": solver_status_label,
        "failure_reason": None if row.get("success") else solver_status_label,
        "secretory_process": mapping.secretory_process or _secretory_process_label(process_code),
        "mapping_level": mapping.mapping_level,
        "mapping_confidence": mapping.mapping_confidence,
        "mapping_interpretation": mapping.interpretation,
        "complex_id": mapping.complex_id or (complex_id or None),
        "complex_subunit_ids": list(mapping.complex_subunit_ids) or [str(item["subunit_id"]) for item in subunits],
        "complex_subunit_stoichiometry": list(mapping.complex_subunit_stoichiometry) or [float(item["stoichiometry"]) for item in subunits],
        "affected_reactions": row.get("affected_reactions") or ([resolved_reaction] if resolved_reaction else []),
        "inactive_reactions": row.get("inactive_reactions") or row.get("inactive_reactions_preview") or [],
        "inactive_reaction_count": int(row.get("inactive_reaction_count") or 0),
        "gpr_rules": row.get("gpr_rules") or [],
        "gpr_role": row.get("gpr_role") or _default_gpr_role(intervention_type),
        "capacity_effect": row.get("capacity_effect") or _default_capacity_effect(intervention_type),
        "simulation_basis": row.get("simulation_basis") or default_basis,
        "ko_support_status": row.get("ko_support_status") or _default_ko_support_status(intervention_type, row.get("status")),
        "oe_support_status": row.get("oe_support_status") or _default_oe_support_status(intervention_type, row.get("status")),
        "support_reason": row.get("support_reason") or _default_support_reason(intervention_type, row.get("status")),
        "missing_information": row.get("missing_information") or [],
        "warnings": row.get("warnings") or [],
    }


def _plan_for_reaction(
    row: dict[str, Any],
    input_gene_ids_by_reaction: dict[str, str] | None,
    gene_intervention_plans_by_gene: dict[str, GeneInterventionPlan] | None,
) -> GeneInterventionPlan | None:
    reaction_id = str(row.get("reaction")) if row.get("reaction") is not None else ""
    gene_text = (input_gene_ids_by_reaction or {}).get(reaction_id)
    if not gene_text:
        return None
    first_gene = gene_text.split(",")[0]
    return (gene_intervention_plans_by_gene or {}).get(first_gene)


def _gene_plan_fields(plan: GeneInterventionPlan | None) -> dict[str, Any]:
    if plan is None:
        return {}
    return dict(plan.candidate_fields())


def _no_effect_gene_knockout_row(
    plan: GeneInterventionPlan,
    baseline_objective_value: float | None,
) -> dict[str, Any]:
    return {
        "gene": plan.gene_id,
        "reaction": None,
        "success": True,
        "status": "no_reaction_disabled",
        "objective_value": baseline_objective_value,
        "delta_vs_baseline": 0.0 if baseline_objective_value is not None else None,
        "inactive_reaction_count": 0,
        "inactive_reactions_preview": [],
    }


def _unresolved_gene_knockout_row(
    plan: GeneInterventionPlan,
    baseline_objective_value: float | None,
) -> dict[str, Any]:
    return {
        "gene": plan.gene_id,
        "reaction": None,
        "success": False,
        "status": "unresolved_gene",
        "objective_value": None,
        "baseline_objective_value": baseline_objective_value,
        "delta_vs_baseline": None,
        "inactive_reaction_count": 0,
        "inactive_reactions_preview": [],
    }


def _explain_only_oe_status(plan: GeneInterventionPlan) -> str:
    if plan.oe_support_status == "oe_no_gpr_effect":
        return "not_run_no_gpr_effect"
    if plan.oe_support_status == "oe_explain_only_no_capacity_model":
        return "not_run_gene_oe_proxy"
    return "not_run_complex_subunit_limited"


def _default_gpr_role(intervention_type: str) -> str:
    if intervention_type in {"KO", "OE_gene_proxy"}:
        return "unresolved"
    return "reaction_level"


def _default_capacity_effect(intervention_type: str) -> str:
    if intervention_type == "KO_reaction":
        return "reaction_disabled"
    if intervention_type == "OE_reaction":
        return "reaction_capacity_proxy"
    if intervention_type == "OE_gene_proxy":
        return "reaction_capacity_proxy"
    return "unknown"


def _default_simulation_basis(intervention_type: str) -> str:
    if intervention_type == "KO":
        return "gpr_gene_deletion"
    if intervention_type == "KO_reaction":
        return "reaction_deletion"
    if intervention_type in {"OE_gene_proxy", "OE_reaction"}:
        return "reaction_level_capacity_proxy"
    return "unknown"


def _default_ko_support_status(intervention_type: str, status: object | None) -> str:
    if intervention_type == "KO":
        if str(status) == "no_reaction_disabled":
            return "ko_no_reaction_disabled"
        return "ko_runnable_gpr_gene_deletion"
    if intervention_type == "KO_reaction":
        return "reaction_level_diagnostic"
    return ""


def _default_oe_support_status(intervention_type: str, status: object | None) -> str:
    if intervention_type == "OE_gene_proxy":
        if str(status) == "not_run_complex_subunit_limited":
            return "oe_explain_only_complex_subunit"
        if str(status) == "not_run_no_gpr_effect":
            return "oe_no_gpr_effect"
        if str(status) == "not_run_gene_oe_proxy":
            return "oe_explain_only_no_capacity_model"
        return "oe_runnable_reaction_proxy"
    if intervention_type == "OE_reaction":
        return "reaction_level_diagnostic"
    return ""


def _default_support_reason(intervention_type: str, status: object | None) -> str:
    if str(status) == "no_reaction_disabled":
        return "Gene deletion leaves all associated model reactions active under GPR AND/OR evaluation."
    if str(status) in {"not_run_complex_subunit_limited", "not_run_gene_oe_proxy"}:
        return "Single-gene OE of a complex subunit is explain-only; it is not a reliable capacity increase."
    if str(status) == "not_run_no_gpr_effect":
        return "Gene exists in the model, but no reaction GPR currently references it."
    if intervention_type == "KO":
        return "Gene KO is simulated by disabling reactions whose GPR rule becomes false."
    if intervention_type in {"OE_gene_proxy", "OE_reaction"}:
        return "OE is represented as a reaction-level capacity proxy."
    return ""


def _relative_delta(delta: float | None, baseline_value: float | None) -> float | None:
    if delta is None or baseline_value is None or abs(baseline_value) < 1e-15:
        return None
    return float(delta) / abs(float(baseline_value))


def _effect_label(effect_code: str, raw_status: object | None = None) -> str:
    if effect_code in {"strong_improvement", "weak_improvement"}:
        return "提升分泌"
    if effect_code in {"strong_decrease", "weak_decrease"}:
        return "降低分泌"
    if effect_code == "neutral":
        return "无明显变化"
    if effect_code == "infeasible_at_fixed_mu":
        if str(raw_status) == "2":
            return "约束不可行"
        if str(raw_status) in {"not_run_complex_subunit_limited", "not_run_gene_oe_proxy", "not_run_no_gpr_effect"}:
            return "未运行"
        if str(raw_status) in {"missing_reaction", "unresolved_gene", "unresolved_reaction"}:
            return "未解析"
        return "求解失败"
    return "未解析"


def _solver_status_label(raw_status: object | None, success: bool) -> str:
    status = "" if raw_status is None else str(raw_status)
    if status == "no_reaction_disabled":
        return "未运行：GPR 未失活任何反应"
    if success:
        return "求解成功"
    labels = {
        "2": "约束不可行",
        "3": "目标无界",
        "4": "求解器数值错误",
        "missing_reaction": "反应未找到",
        "unresolved_gene": "基因未解析",
        "unresolved_reaction": "反应未解析",
        "not_run_complex_subunit_limited": "仅解释，未求解",
        "not_run_gene_oe_proxy": "仅解释，未求解",
        "not_run_no_gpr_effect": "仅解释，模型无 GPR 影响",
        "no_reaction_disabled": "未运行：GPR 未失活任何反应",
        "missing_objective": "目标反应未找到",
    }
    return labels.get(status, "求解失败")


def _secretory_process_label(process_code: str) -> str:
    labels = {
        "ribosome": "翻译",
        "proteasome_degradation": "蛋白降解",
        "disulfide_folding": "ER 折叠 / DSB",
        "n_glycan_processing": "N-糖基化 NG",
        "o_glycan_processing": "O-糖基化 OG",
        "chaperone_folding": "ER 折叠 / 分子伴侣",
        "erad_misfolding": "错误折叠 / ERAD",
        "er_translocation": "ER 转运",
        "er_to_golgi_transport": "ER 到 Golgi 转运",
        "golgi_surface_transport": "Golgi 到胞外运输",
        "secretory_capacity": "分泌容量",
        "metabolic_or_other": "代谢或其它反应",
        "unknown": "未解析",
    }
    return labels.get(process_code, process_code)


def _screen_result(
    target_id: str,
    screen_type: str,
    rows: tuple[dict[str, Any], ...],
    prepared: dict[str, Any],
) -> ScreenResult:
    return ScreenResult(
        target_id=target_id,
        screen_type=screen_type,
        success=bool(rows) and all(bool(row.get("success")) for row in rows),
        candidate_count=len(rows),
        rows=rows,
        constraint_counts=prepared["constraint_counts"],
        baseline_objective_value=prepared["baseline"].objective_value,
        result_status="draft",
        matlab_alignment_status="pending",
    )


def _empty_screen_result(target_id: str, screen_type: str, prepared: dict[str, Any]) -> ScreenResult:
    baseline = prepared.get("baseline")
    return ScreenResult(
        target_id=target_id,
        screen_type=screen_type,
        success=False,
        candidate_count=0,
        rows=(),
        constraint_counts=prepared.get("constraint_counts", {}),
        baseline_objective_value=getattr(baseline, "objective_value", None),
        result_status="draft",
        matlab_alignment_status="pending",
    )


def _empty_unsolved_screen_result(target_id: str, screen_type: str) -> ScreenResult:
    return ScreenResult(
        target_id=target_id,
        screen_type=screen_type,
        success=False,
        candidate_count=0,
        rows=(),
        constraint_counts={},
        baseline_objective_value=None,
        result_status="draft_skipped_empty_screen",
        matlab_alignment_status="pending",
    )


__all__ = [
    "ScreenResult",
    "ScreenPlanResult",
    "GeneCapabilityProfile",
    "GenePerturbationMapping",
    "GenePerturbationMapResult",
    "GeneReactionMapping",
    "GeneInterventionPlan",
    "build_all_gene_capability_catalog",
    "build_gene_capability_profile",
    "build_gene_perturbation_map",
    "build_reaction_perturbation_mapping",
    "build_screen_plan",
    "default_ko_genes",
    "default_oe_reactions",
    "explain_only_gene_overexpression_rows",
    "plan_gene_knockout",
    "plan_gene_overexpression",
    "reactions_for_gene",
    "resolve_oe_gene_reactions",
    "split_existing_genes",
    "run_knockout_screen",
    "run_ko_screen",
    "run_overexpression_screen",
    "run_oe_screen",
    "run_reaction_knockout_screen",
    "run_pcsec_growth_tradeoff",
    "run_pcsec_ko_screen",
    "run_pcsec_oe_screen",
    "run_pcsec_reaction_ko_screen",
    "split_existing_reactions",
    "summarize_screen_result",
]
