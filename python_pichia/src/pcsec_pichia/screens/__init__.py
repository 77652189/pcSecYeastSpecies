from __future__ import annotations

from collections.abc import Iterable
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
from pcsec_pichia.strain_modifications import StrainModifications, apply_strain_modifications


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


# R2 (ADR-004): default OE capacity multipliers to sweep for a dose-response shape. Spans
# modest to aggressive over-expression and includes the legacy single 2.0x so it stays
# visible as one point on the curve. 1.0 is the no-OE baseline and is added separately.
DEFAULT_OE_DOSE_RESPONSE_FACTORS: tuple[float, ...] = (1.25, 1.5, 2.0, 3.0, 5.0, 8.0)


@dataclass(frozen=True)
class OeDoseResponseSweepResult:
    """R2 (ADR-004) raw OE factor sweep for a set of reactions (screens layer, no shape logic).

    Re-solves the same target LP while scaling each reaction's OE capacity multiplier across a
    factor grid, so the caller (analysis.classify_oe_dose_response_sweep) can classify the
    dose-response *shape*. It carries only relative objective values/deltas; it never produces
    an absolute capacity. Shape classification is deliberately kept out of this layer.
    """

    target_id: str
    enabled: bool
    success: bool
    tested_factors: tuple[float, ...]
    baseline_objective: float | None
    reactions: tuple[str, ...]
    reaction_points: dict[str, tuple[tuple[float, float | None], ...]]
    sweep_rows: tuple[dict[str, Any], ...]
    result_status: str
    warnings: tuple[str, ...]
    matlab_alignment_status: str = "pending"


def _apply_modifications_and_rebaseline(
    prepared: dict[str, Any],
    modifications: StrainModifications,
    metabolic: MetabolicEnzymeData,
    growth_rate: float,
    write_ribosome_translation_constraint: bool,
    write_misfolding_constraints: bool,
) -> dict[str, Any]:
    """把改造（stacked KO/OE）应用到已 prepare 的野生型基线上，并重解**改造后基线**。

    返回覆盖了 fixed_model/secretory/combined/baseline/baseline_success 的 prepared 副本，供 KO
    筛查、OE 剂量响应等在**改造后菌株**上叠加单扰动（ADR-004 #1 迭代2 分层复用）。改造后不可行则
    baseline_success=False。空改造时调用方不应进入此函数（默认 None → 野生型路径 byte-identical）。
    """
    fixed_model, secretory_eff, combined_eff, _applied, _warnings = apply_strain_modifications(
        prepared["fixed_model"], prepared["secretory"], prepared["combined"], modifications
    )
    baseline, _counts = solve_pcsec_maximize(
        fixed_model,
        prepared["exchange_reaction_id"],
        metabolic=metabolic,
        secretory=secretory_eff,
        combined=combined_eff,
        mu=growth_rate,
        key_reactions=("BIOMASS", "Ex_glc_D", prepared["exchange_reaction_id"]),
        write_ribosome_translation_constraint=write_ribosome_translation_constraint,
        write_misfolding_constraints=write_misfolding_constraints,
    )
    updated = dict(prepared)
    updated["fixed_model"] = fixed_model
    updated["secretory"] = secretory_eff
    updated["combined"] = combined_eff
    updated["baseline"] = baseline
    updated["baseline_success"] = bool(getattr(baseline, "success", False))
    return updated


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
    strain_modifications: StrainModifications | None = None,
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

    # ADR-004 #1 迭代2：改造后 KO 候选——把 stacked KO/OE 叠进基线再逐个测 KO，delta 相对**改造后**
    # 菌株（不是野生型）。默认 None → 野生型路径 byte-identical。改造后不可行则优雅返回空。
    if strain_modifications is not None and not strain_modifications.is_empty():
        prepared = _apply_modifications_and_rebaseline(
            prepared, strain_modifications, metabolic, growth_rate,
            write_ribosome_translation_constraint, write_misfolding_constraints,
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


def run_oe_dose_response_sweep(
    model: CobraModel,
    target: TargetSpec,
    amino_acids: AminoAcidStoichiometry,
    metabolic: MetabolicEnzymeData,
    secretory: SecretoryEnzymeData,
    combined: CombinedEnzymeData,
    reactions: list[str],
    factors: Iterable[float] = DEFAULT_OE_DOSE_RESPONSE_FACTORS,
    growth_rate: float = 0.10,
    write_ribosome_translation_constraint: bool = False,
    write_misfolding_constraints: bool = False,
    strain_modifications: StrainModifications | None = None,
) -> OeDoseResponseSweepResult:
    """R2 (ADR-004): sweep the OE capacity multiplier over a factor grid for each reaction.

    Prepares the target LP once (shared baseline) and re-solves each reaction at every factor
    via run_pcsec_oe_screen. Returns the raw (factor, objective) points per reaction anchored at
    the no-OE baseline (factor 1.0); it does not classify the shape (see analysis layer).

    `strain_modifications` (opt-in, ADR-004 #1 迭代候选): when set, the stacked KO/OE are applied to
    the prepared strain and the factor-1.0 anchor is re-solved on the *modified* strain, so each
    swept reaction's dose-response is measured on top of the already-modified strain (谁松开能再涨、
    涨到哪饱和). `None`/empty leaves the wildtype-baseline sweep byte-identical.
    """

    warnings = (
        "OE dose-response sweep re-solves the target LP at several capacity multipliers; it is an "
        "opt-in relative probe and does not change the default single-run objective.",
        "Objective values are relative model secretion, not absolute titers; a factor is a capacity "
        "multiplier, not a measured expression level.",
    )
    sweep_factors = tuple(sorted({float(f) for f in factors if float(f) > 0.0 and abs(float(f) - 1.0) > 1e-9}))
    unique_reactions = tuple(dict.fromkeys(str(r) for r in reactions))
    if not unique_reactions or not sweep_factors:
        return OeDoseResponseSweepResult(
            target_id=target.target_id,
            enabled=True,
            success=False,
            tested_factors=sweep_factors,
            baseline_objective=None,
            reactions=unique_reactions,
            reaction_points={},
            sweep_rows=(),
            result_status="draft_oe_dose_response_no_reactions",
            warnings=(*warnings, "No reactions or no OE factors > 1.0 to sweep."),
        )

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
        return OeDoseResponseSweepResult(
            target_id=target.target_id,
            enabled=True,
            success=False,
            tested_factors=sweep_factors,
            baseline_objective=None,
            reactions=unique_reactions,
            reaction_points={},
            sweep_rows=(),
            result_status="draft_oe_dose_response_unavailable",
            warnings=(*warnings, "Baseline secretion solve did not succeed; cannot build a dose-response."),
        )

    if strain_modifications is not None and not strain_modifications.is_empty():
        # 在改造后菌株上重锚 factor 1.0，使每个候选的剂量响应都测在已改造菌株之上（与 KO 筛查同款 helper）。
        prepared = _apply_modifications_and_rebaseline(
            prepared, strain_modifications, metabolic, growth_rate,
            write_ribosome_translation_constraint, write_misfolding_constraints,
        )
    fixed_model = prepared["fixed_model"]
    secretory_effective = prepared["secretory"]
    combined_effective = prepared["combined"]
    baseline = prepared["baseline"]
    baseline_objective = baseline.objective_value
    points_by_reaction: dict[str, list[tuple[float, float | None]]] = {
        reaction_id: [(1.0, baseline_objective)] for reaction_id in unique_reactions
    }
    sweep_rows: list[dict[str, Any]] = []
    for factor in sweep_factors:
        raw_rows = run_pcsec_oe_screen(
            fixed_model,
            baseline,
            list(unique_reactions),
            prepared["exchange_reaction_id"],
            metabolic=metabolic,
            secretory=secretory_effective,
            combined=combined_effective,
            mu=growth_rate,
            factor=factor,
            write_ribosome_translation_constraint=write_ribosome_translation_constraint,
            write_misfolding_constraints=write_misfolding_constraints,
        )
        for row in raw_rows:
            reaction_id = str(row.get("reaction"))
            objective_value = row.get("objective_value") if row.get("success") else None
            points_by_reaction.setdefault(reaction_id, []).append((factor, objective_value))
            sweep_rows.append(
                {
                    "reaction": reaction_id,
                    "factor": factor,
                    "objective_value": objective_value,
                    "delta_vs_baseline": row.get("delta_vs_baseline"),
                    "success": bool(row.get("success")),
                    "status": row.get("status"),
                    "capacity_basis": row.get("capacity_basis"),
                }
            )

    reaction_points = {
        reaction_id: tuple(points) for reaction_id, points in points_by_reaction.items()
    }
    return OeDoseResponseSweepResult(
        target_id=target.target_id,
        enabled=True,
        success=baseline_objective is not None,
        tested_factors=sweep_factors,
        baseline_objective=baseline_objective,
        reactions=unique_reactions,
        reaction_points=reaction_points,
        sweep_rows=tuple(sweep_rows),
        result_status="draft_oe_dose_response",
        warnings=warnings,
    )


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


def prepare_screen_inputs(
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
    """Prepare one target with the canonical screen baseline path.

    Service facades that need a prepared model should use this public adapter
    instead of depending on the private implementation helper below.
    """

    return _prepare_screen_inputs(
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
        "external_model_sources": row.get("external_model_sources") or [],
        "gpr_source_priority": row.get("gpr_source_priority") or {},
        "external_gpr_candidate_count": int(row.get("external_gpr_candidate_count") or 0),
        "best_external_gpr_source": row.get("best_external_gpr_source") or "",
        "external_gpr_mapping_status": row.get("external_gpr_mapping_status") or {},
        "external_gpr_conflict_warnings": row.get("external_gpr_conflict_warnings") or [],
        "manual_review_reasons": row.get("manual_review_reasons") or [],
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
    "DEFAULT_OE_DOSE_RESPONSE_FACTORS",
    "OeDoseResponseSweepResult",
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
    "run_oe_dose_response_sweep",
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
