from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable

import numpy as np
from scipy import sparse
from scipy.optimize import linprog

from pcsec_pichia.probe import (
    AminoAcidStoichiometry,
    CobraModel,
    CombinedEnzymeData,
    MetabolicEnzymeData,
    SecretoryEnzymeData,
    SolveResult,
    TargetSpec,
    build_pcsec_constraint_matrices,
    build_supported_target_model,
    build_target_enzymedata,
    run_pcsec_growth_tradeoff,
    solve_maximize,
    solve_pcsec_maximize,
    target_secretion_smoke,
)
from pcsec_pichia.media import CarbonSourceFormulation, MATLAB_LEGACY_COST_MEDIUM_EXCHANGES, list_carbon_source_formulations


@dataclass(frozen=True)
class SecretionSimulationResult:
    success: bool
    target_id: str
    objective_value: float | None
    growth_rate: float
    secretion_flux: float | None
    status: str
    message: str
    constraint_counts: dict[str, int]
    result_status: str
    target_parameter_status: str
    matlab_alignment_status: str
    exchange_reaction_id: str | None
    build_status: str
    lp_sensitivity: dict[str, tuple[float, ...]] | None = None
    key_fluxes: dict[str, float] | None = None
    growth_reaction_id: str = "BIOMASS"
    open_growth_reaction_ids: tuple[str, ...] = ("BIOMASS",)
    growth_reaction_status: str = "single_growth_reaction"


@dataclass(frozen=True)
class GrowthTradeoffResult:
    target_id: str
    success: bool
    growth_points: tuple[float, ...]
    tradeoff_rows: tuple[dict[str, Any], ...]
    result_status: str
    target_parameter_status: str
    matlab_alignment_status: str
    growth_reaction_id: str = "BIOMASS"
    open_growth_reaction_ids: tuple[str, ...] = ("BIOMASS",)
    growth_reaction_status: str = "single_growth_reaction"


@dataclass(frozen=True)
class ProteinCostSlopeCompatibilityResult:
    target_id: str
    enabled: bool
    success: bool
    growth_rates: tuple[float, ...]
    secretion_ratios: tuple[float, ...]
    rows: tuple[dict[str, Any], ...]
    glucose_cost_slopes: tuple[dict[str, Any], ...]
    ribosome_cost_slopes: tuple[dict[str, Any], ...]
    result_status: str
    warnings: tuple[str, ...]
    medium_compatibility_mode: str = "corrected"
    medium_bound_overrides: tuple[dict[str, Any], ...] = ()
    secretion_ratio_policy: str = "explicit_absolute_ratios"
    capacity_reference: float | None = None
    capacity_fractions: tuple[float, ...] = ()


@dataclass(frozen=True)
class MixedCarbonObjectiveResult:
    target_id: str
    success: bool
    objective_mode: str
    carbon_weights: dict[str, float]
    growth_rate: float
    target_exchange_ratio: float
    objective_value: float | None
    target_exchange_reaction: str | None
    carbon_fluxes: dict[str, float]
    carbon_costs: dict[str, float]
    total_weighted_carbon_cost: float | None
    status: str
    message: str
    constraint_counts: dict[str, int]
    result_status: str
    warnings: tuple[str, ...]
    growth_reaction_id: str = "BIOMASS"
    open_growth_reaction_ids: tuple[str, ...] = ("BIOMASS",)
    growth_reaction_status: str = "single_growth_reaction"


def solve_secretion_capacity(
    model: CobraModel,
    target: TargetSpec,
    amino_acids: AminoAcidStoichiometry,
    metabolic: MetabolicEnzymeData,
    secretory: SecretoryEnzymeData,
    combined: CombinedEnzymeData,
    growth_rate: float = 0.10,
    write_ribosome_translation_constraint: bool = False,
    write_misfolding_constraints: bool = False,
) -> SecretionSimulationResult:
    """Solve fixed-growth pcSec secretion capacity for one target."""

    prepared = _prepare_target_pcsec_inputs(model, target, amino_acids, secretory, combined)
    if prepared["model"] is None or prepared["exchange_reaction_id"] is None:
        return SecretionSimulationResult(
            success=False,
            target_id=target.target_id,
            objective_value=None,
            growth_rate=float(growth_rate),
            secretion_flux=None,
            status=str(prepared["build_status"]),
            message=str(prepared["reason"]),
            constraint_counts={},
            result_status="draft",
            target_parameter_status=_target_parameter_status(target),
            matlab_alignment_status="pending",
            exchange_reaction_id=None,
            build_status=str(prepared["build_status"]),
            lp_sensitivity=None,
            key_fluxes={},
        )

    growth_context = _growth_reaction_context(prepared["model"])
    fixed_model = prepared["model"].with_bounds(_fixed_growth_bounds(growth_context, growth_rate))
    solved, counts = solve_pcsec_maximize(
        fixed_model,
        prepared["exchange_reaction_id"],
        metabolic=metabolic,
        secretory=prepared["secretory"],
        combined=prepared["combined"],
        mu=growth_rate,
        key_reactions=(
            growth_context["growth_reaction_id"],
            *growth_context["open_growth_reaction_ids"],
            "Ex_glc_D",
            "Ex_glyc",
            "Ex_meoh",
            "Ex_o2",
            prepared["exchange_reaction_id"],
        ),
        write_ribosome_translation_constraint=write_ribosome_translation_constraint,
        write_misfolding_constraints=write_misfolding_constraints,
    )
    return SecretionSimulationResult(
        success=solved.success,
        target_id=target.target_id,
        objective_value=solved.objective_value,
        growth_rate=float(growth_rate),
        secretion_flux=solved.objective_value if solved.success else None,
        status=solved.status,
        message=solved.message,
        constraint_counts={str(key): int(value) for key, value in counts.items()},
        result_status="draft",
        target_parameter_status=_target_parameter_status(target),
        matlab_alignment_status="pending",
        exchange_reaction_id=prepared["exchange_reaction_id"],
        build_status=str(prepared["build_status"]),
        lp_sensitivity=solved.sensitivity,
        key_fluxes=dict(solved.fluxes),
        growth_reaction_id=str(growth_context["growth_reaction_id"]),
        open_growth_reaction_ids=tuple(str(item) for item in growth_context["open_growth_reaction_ids"]),
        growth_reaction_status=str(growth_context["growth_reaction_status"]),
    )


def run_growth_tradeoff(
    model: CobraModel,
    target: TargetSpec,
    amino_acids: AminoAcidStoichiometry,
    metabolic: MetabolicEnzymeData,
    secretory: SecretoryEnzymeData,
    combined: CombinedEnzymeData,
    growth_points: Iterable[float],
    write_ribosome_translation_constraint: bool = False,
    write_misfolding_constraints: bool = False,
) -> GrowthTradeoffResult:
    """Run a small fixed-growth pcSec secretion tradeoff grid for one target."""

    prepared = _prepare_target_pcsec_inputs(model, target, amino_acids, secretory, combined)
    normalized_points = tuple(sorted({float(value) for value in growth_points if float(value) > 0}))
    if prepared["model"] is None or prepared["exchange_reaction_id"] is None:
        return GrowthTradeoffResult(
            target_id=target.target_id,
            success=False,
            growth_points=normalized_points,
            tradeoff_rows=(),
            result_status="draft",
            target_parameter_status=_target_parameter_status(target),
            matlab_alignment_status="pending",
        )

    growth_context = _growth_reaction_context(prepared["model"])
    rows: list[dict[str, object]] = []
    for mu in normalized_points:
        fixed_model = prepared["model"].with_bounds(_fixed_growth_bounds(growth_context, mu))
        solved, counts = solve_pcsec_maximize(
            fixed_model,
            prepared["exchange_reaction_id"],
            metabolic=metabolic,
            secretory=prepared["secretory"],
            combined=prepared["combined"],
            mu=mu,
            key_reactions=(
                growth_context["growth_reaction_id"],
                *growth_context["open_growth_reaction_ids"],
                "Ex_glc_D",
                "Ex_glyc",
                "Ex_meoh",
                "Ex_o2",
                prepared["exchange_reaction_id"],
            ),
            write_ribosome_translation_constraint=write_ribosome_translation_constraint,
            write_misfolding_constraints=write_misfolding_constraints,
        )
        secretion = solved.objective_value if solved.success else None
        rows.append(
            {
                "mu": mu,
                "success": solved.success,
                "status": solved.status,
                "secretion_flux": secretion,
                "secretion_per_biomass": secretion / mu if secretion is not None and mu > 0 else None,
                "message": solved.message,
                "constraint_counts": counts,
                "growth_reaction_id": growth_context["growth_reaction_id"],
                "open_growth_reaction_ids": growth_context["open_growth_reaction_ids"],
                "growth_reaction_status": growth_context["growth_reaction_status"],
            }
        )
    return GrowthTradeoffResult(
        target_id=target.target_id,
        success=bool(rows) and all(bool(row.get("success")) for row in rows),
        growth_points=normalized_points,
        tradeoff_rows=tuple(rows),
        result_status="draft",
        target_parameter_status=_target_parameter_status(target),
        matlab_alignment_status="pending",
        growth_reaction_id=str(growth_context["growth_reaction_id"]),
        open_growth_reaction_ids=tuple(str(item) for item in growth_context["open_growth_reaction_ids"]),
        growth_reaction_status=str(growth_context["growth_reaction_status"]),
    )


def run_protein_cost_slope_compatibility(
    model: CobraModel,
    target: TargetSpec,
    amino_acids: AminoAcidStoichiometry,
    metabolic: MetabolicEnzymeData,
    secretory: SecretoryEnzymeData,
    combined: CombinedEnzymeData,
    growth_rates: Iterable[float] = (0.05, 0.10),
    secretion_ratios: Iterable[float] = (5e-7, 1e-6, 5e-6, 1e-5, 2e-5),
    write_ribosome_translation_constraint: bool = False,
    write_misfolding_constraints: bool = False,
    medium_compatibility_mode: str = "corrected",
) -> ProteinCostSlopeCompatibilityResult:
    """Run optional MATLAB-style protein cost slope probes.

    This fixes target secretion and growth, then maximizes Ex_glc_D with the
    existing pcSec LP. It is an opt-in compatibility probe, not the default
    corrected secretion-capacity objective.
    """

    normalized_mu = tuple(sorted({float(value) for value in growth_rates if float(value) > 0}))
    normalized_ratios = tuple(sorted({float(value) for value in secretion_ratios if float(value) > 0}))
    normalized_medium_mode = _normalize_cost_slope_medium_mode(medium_compatibility_mode)
    prepared = _prepare_target_pcsec_inputs(model, target, amino_acids, secretory, combined)
    warnings = [
        "MATLAB-compatible cost slope mode is an opt-in Python draft probe; it does not replace the default corrected pipeline.",
        "Definition: fix target exchange ratios and growth, then optimize Ex_glc_D to estimate glucose cost slopes.",
        "Ribosome slope uses Mach_Ribosome_complex_formation flux when that reaction is available; otherwise it is reported as unavailable.",
    ]
    if prepared["model"] is None or prepared["exchange_reaction_id"] is None:
        return ProteinCostSlopeCompatibilityResult(
            target_id=target.target_id,
            enabled=True,
            success=False,
            growth_rates=normalized_mu,
            secretion_ratios=normalized_ratios,
            rows=(),
            glucose_cost_slopes=(),
            ribosome_cost_slopes=(),
            result_status="draft_cost_slope_unavailable",
            warnings=tuple([*warnings, str(prepared["reason"])]),
            medium_compatibility_mode=normalized_medium_mode,
            medium_bound_overrides=(),
        )

    exchange_reaction_id = str(prepared["exchange_reaction_id"])
    ribosome_reaction_id = "Mach_Ribosome_complex_formation"
    medium_bounds, medium_overrides, medium_warnings = _cost_slope_medium_bound_overrides(
        prepared["model"],
        normalized_medium_mode,
    )
    warnings.extend(medium_warnings)
    growth_context = _growth_reaction_context(prepared["model"])
    warnings.extend(_growth_context_warnings(growth_context))
    rows: list[dict[str, Any]] = []
    for mu in normalized_mu:
        for ratio in normalized_ratios:
            fixed_model = prepared["model"].with_bounds(
                {
                    **_fixed_growth_bounds(growth_context, mu),
                    exchange_reaction_id: (ratio, ratio),
                    **medium_bounds,
                }
            )
            key_reactions = (
                growth_context["growth_reaction_id"],
                *growth_context["open_growth_reaction_ids"],
                "Ex_glc_D",
                "Ex_glyc",
                "Ex_meoh",
                "Ex_o2",
                exchange_reaction_id,
                ribosome_reaction_id,
            )
            solved, counts = solve_pcsec_maximize(
                fixed_model,
                "Ex_glc_D",
                metabolic=metabolic,
                secretory=prepared["secretory"],
                combined=prepared["combined"],
                mu=mu,
                key_reactions=key_reactions,
                write_ribosome_translation_constraint=write_ribosome_translation_constraint,
                write_misfolding_constraints=write_misfolding_constraints,
            )
            glucose_flux = _optional_flux(solved.fluxes.get("Ex_glc_D"))
            glucose_cost = _uptake_cost(glucose_flux)
            ribosome_flux = _optional_flux(solved.fluxes.get(ribosome_reaction_id))
            rows.append(
                {
                    "mu": mu,
                    "target_exchange_ratio": ratio,
                    "success": solved.success,
                    "status": solved.status,
                    "message": solved.message,
                    "objective_reaction": "Ex_glc_D",
                    "objective_value": solved.objective_value,
                    "target_exchange_reaction": exchange_reaction_id,
                    "glucose_flux": glucose_flux,
                    "glucose_cost": glucose_cost,
                    "glucose_cost_status": _uptake_cost_status(glucose_flux),
                    "ribosome_reaction": ribosome_reaction_id,
                    "ribosome_flux": ribosome_flux,
                    "ribosome_cost": abs(ribosome_flux) if ribosome_flux is not None else None,
                    "constraint_counts": counts,
                    "medium_compatibility_mode": normalized_medium_mode,
                    "growth_reaction_id": growth_context["growth_reaction_id"],
                    "open_growth_reaction_ids": growth_context["open_growth_reaction_ids"],
                    "growth_reaction_status": growth_context["growth_reaction_status"],
                }
            )

    if any(row.get("glucose_cost_status") == "non_uptake_flux" for row in rows):
        warnings.append(
            "At least one cost-slope row produced positive Ex_glc_D flux; that row is not treated as glucose uptake cost."
        )
    glucose_slopes = _cost_slopes(rows, "glucose_cost")
    ribosome_slopes = _cost_slopes(rows, "ribosome_cost")
    comparable_rows = [row for row in rows if row.get("success")]
    return ProteinCostSlopeCompatibilityResult(
        target_id=target.target_id,
        enabled=True,
        success=bool(comparable_rows) and bool(glucose_slopes),
        growth_rates=normalized_mu,
        secretion_ratios=normalized_ratios,
        rows=tuple(rows),
        glucose_cost_slopes=tuple(glucose_slopes),
        ribosome_cost_slopes=tuple(ribosome_slopes),
        result_status="draft_matlab_compatible_cost_slope",
        warnings=tuple(warnings),
        medium_compatibility_mode=normalized_medium_mode,
        medium_bound_overrides=tuple(medium_overrides),
    )


def run_mixed_carbon_objective_probe(
    model: CobraModel,
    target: TargetSpec,
    amino_acids: AminoAcidStoichiometry,
    metabolic: MetabolicEnzymeData,
    secretory: SecretoryEnzymeData,
    combined: CombinedEnzymeData,
    *,
    growth_rate: float = 0.10,
    target_exchange_ratio: float,
    carbon_weights: dict[str, float],
    objective_mode: str = "weighted_carbon_uptake",
    write_ribosome_translation_constraint: bool = False,
    write_misfolding_constraints: bool = False,
) -> MixedCarbonObjectiveResult:
    """Opt-in mixed-carbon probe; not used by the corrected default pipeline."""

    if objective_mode != "weighted_carbon_uptake":
        raise ValueError(f"Unsupported mixed-carbon objective mode: {objective_mode}")
    normalized_weights = {
        str(reaction_id): float(weight)
        for reaction_id, weight in carbon_weights.items()
        if float(weight) > 0
    }
    warnings = [
        "Mixed-carbon objective probe is opt-in and does not change the corrected default pipeline.",
        "It minimizes a weighted uptake-cost proxy after fixing growth and target exchange.",
    ]
    if not normalized_weights:
        return MixedCarbonObjectiveResult(
            target_id=target.target_id,
            success=False,
            objective_mode=objective_mode,
            carbon_weights={},
            growth_rate=float(growth_rate),
            target_exchange_ratio=float(target_exchange_ratio),
            objective_value=None,
            target_exchange_reaction=None,
            carbon_fluxes={},
            carbon_costs={},
            total_weighted_carbon_cost=None,
            status="missing_carbon_weights",
            message="At least one positive carbon objective weight is required.",
            constraint_counts={},
            result_status="draft_mixed_carbon_objective_unavailable",
            warnings=tuple(warnings),
        )

    prepared = _prepare_target_pcsec_inputs(model, target, amino_acids, secretory, combined)
    if prepared["model"] is None or prepared["exchange_reaction_id"] is None:
        return MixedCarbonObjectiveResult(
            target_id=target.target_id,
            success=False,
            objective_mode=objective_mode,
            carbon_weights=normalized_weights,
            growth_rate=float(growth_rate),
            target_exchange_ratio=float(target_exchange_ratio),
            objective_value=None,
            target_exchange_reaction=None,
            carbon_fluxes={},
            carbon_costs={},
            total_weighted_carbon_cost=None,
            status=str(prepared["build_status"]),
            message=str(prepared["reason"]),
            constraint_counts={},
            result_status="draft_mixed_carbon_objective_unavailable",
            warnings=tuple(warnings),
        )

    missing = tuple(reaction_id for reaction_id in normalized_weights if reaction_id not in prepared["model"].reaction_index)
    if missing:
        return MixedCarbonObjectiveResult(
            target_id=target.target_id,
            success=False,
            objective_mode=objective_mode,
            carbon_weights=normalized_weights,
            growth_rate=float(growth_rate),
            target_exchange_ratio=float(target_exchange_ratio),
            objective_value=None,
            target_exchange_reaction=str(prepared["exchange_reaction_id"]),
            carbon_fluxes={},
            carbon_costs={},
            total_weighted_carbon_cost=None,
            status="missing_carbon_reactions",
            message=f"Carbon objective reactions not found: {', '.join(missing)}",
            constraint_counts={},
            result_status="draft_mixed_carbon_objective_unavailable",
            warnings=tuple(warnings),
        )

    growth_context = _growth_reaction_context(prepared["model"])
    warnings.extend(_growth_context_warnings(growth_context))
    fixed_model = prepared["model"].with_bounds(
        {
            **_fixed_growth_bounds(growth_context, float(growth_rate)),
            str(prepared["exchange_reaction_id"]): (float(target_exchange_ratio), float(target_exchange_ratio)),
        }
    )
    solved, counts = _solve_weighted_carbon_uptake_minimize(
        fixed_model,
        carbon_weights=normalized_weights,
        metabolic=metabolic,
        secretory=prepared["secretory"],
        combined=prepared["combined"],
        mu=float(growth_rate),
        key_reactions=(
            *normalized_weights.keys(),
            str(prepared["exchange_reaction_id"]),
            growth_context["growth_reaction_id"],
            *growth_context["open_growth_reaction_ids"],
            "Ex_o2",
        ),
        write_ribosome_translation_constraint=write_ribosome_translation_constraint,
        write_misfolding_constraints=write_misfolding_constraints,
    )
    carbon_fluxes = {
        reaction_id: float(solved.fluxes[reaction_id])
        for reaction_id in normalized_weights
        if reaction_id in solved.fluxes
    }
    carbon_costs = {
        reaction_id: float(solved.fluxes.get(f"uptake_cost::{reaction_id}", 0.0))
        for reaction_id in normalized_weights
    }
    total_cost = (
        sum(normalized_weights[reaction_id] * carbon_costs.get(reaction_id, 0.0) for reaction_id in normalized_weights)
        if solved.success
        else None
    )
    return MixedCarbonObjectiveResult(
        target_id=target.target_id,
        success=solved.success,
        objective_mode=objective_mode,
        carbon_weights=normalized_weights,
        growth_rate=float(growth_rate),
        target_exchange_ratio=float(target_exchange_ratio),
        objective_value=solved.objective_value,
        target_exchange_reaction=str(prepared["exchange_reaction_id"]),
        carbon_fluxes=carbon_fluxes,
        carbon_costs=carbon_costs,
        total_weighted_carbon_cost=total_cost,
        status=solved.status,
        message=solved.message,
        constraint_counts={str(key): int(value) for key, value in counts.items()},
        result_status="draft_mixed_carbon_objective",
        warnings=tuple(warnings),
        growth_reaction_id=str(growth_context["growth_reaction_id"]),
        open_growth_reaction_ids=tuple(str(item) for item in growth_context["open_growth_reaction_ids"]),
        growth_reaction_status=str(growth_context["growth_reaction_status"]),
    )


def summarize_mixed_carbon_objective_result(result: MixedCarbonObjectiveResult) -> dict[str, Any]:
    return {
        "target_id": result.target_id,
        "success": result.success,
        "objective_mode": result.objective_mode,
        "carbon_weights": result.carbon_weights,
        "growth_rate": result.growth_rate,
        "target_exchange_ratio": result.target_exchange_ratio,
        "objective_value": result.objective_value,
        "target_exchange_reaction": result.target_exchange_reaction,
        "carbon_fluxes": result.carbon_fluxes,
        "carbon_costs": result.carbon_costs,
        "total_weighted_carbon_cost": result.total_weighted_carbon_cost,
        "status": result.status,
        "message": result.message,
        "constraint_counts": result.constraint_counts,
        "result_status": result.result_status,
        "warnings": result.warnings,
        "growth_reaction_id": result.growth_reaction_id,
        "open_growth_reaction_ids": result.open_growth_reaction_ids,
        "growth_reaction_status": result.growth_reaction_status,
    }


def summarize_simulation_result(result: SecretionSimulationResult | GrowthTradeoffResult) -> dict[str, Any]:
    if isinstance(result, SecretionSimulationResult):
        return {
            "success": result.success,
            "target_id": result.target_id,
            "objective_value": result.objective_value,
            "growth_rate": result.growth_rate,
            "secretion_flux": result.secretion_flux,
            "status": result.status,
            "message": result.message,
            "constraint_counts": result.constraint_counts,
            "result_status": result.result_status,
            "target_parameter_status": result.target_parameter_status,
            "matlab_alignment_status": result.matlab_alignment_status,
            "exchange_reaction_id": result.exchange_reaction_id,
            "build_status": result.build_status,
            "growth_reaction_id": result.growth_reaction_id,
            "open_growth_reaction_ids": result.open_growth_reaction_ids,
            "growth_reaction_status": result.growth_reaction_status,
        }
    return {
        "success": result.success,
        "target_id": result.target_id,
        "growth_points": result.growth_points,
        "tradeoff_rows": result.tradeoff_rows,
        "result_status": result.result_status,
        "target_parameter_status": result.target_parameter_status,
        "matlab_alignment_status": result.matlab_alignment_status,
        "growth_reaction_id": result.growth_reaction_id,
        "open_growth_reaction_ids": result.open_growth_reaction_ids,
        "growth_reaction_status": result.growth_reaction_status,
    }


def _optional_flux(value: object) -> float | None:
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


_BIOMASS_REACTION_IDS: tuple[str, ...] = ("BIOMASS", "BIOMASS_glyc", "BIOMASS_meoh")


def _growth_reaction_context(model: CobraModel) -> dict[str, object]:
    open_reactions = _open_growth_reactions(model)
    formulation = _matching_carbon_source_formulation(model, open_reactions)
    selected = _select_growth_reaction(model, open_reactions, formulation)
    status = "single_growth_reaction" if len(open_reactions) <= 1 else "multiple_growth_reactions_selected"
    if not open_reactions:
        status = "fallback_growth_reaction"
    return {
        "growth_reaction_id": selected,
        "open_growth_reaction_ids": open_reactions or (selected,),
        "growth_reaction_status": status,
        "carbon_source_id": None if formulation is None else formulation.carbon_source_id,
        "carbon_source_formulation_status": None if formulation is None else formulation.formulation_status,
        "carbon_objective_weights": {} if formulation is None else formulation.carbon_objective_weights,
        "formulation_warnings": () if formulation is None else formulation.warnings,
    }


def _fixed_growth_bounds(growth_context: dict[str, object], growth_rate: float) -> dict[str, tuple[float, float]]:
    selected = str(growth_context["growth_reaction_id"])
    open_reactions = tuple(str(item) for item in growth_context["open_growth_reaction_ids"])
    bounds = {selected: (float(growth_rate), float(growth_rate))}
    for reaction_id in open_reactions:
        if reaction_id != selected:
            bounds[reaction_id] = (0.0, 0.0)
    return bounds


def _growth_context_warnings(growth_context: dict[str, object]) -> list[str]:
    warnings = [str(item) for item in growth_context.get("formulation_warnings", ())]
    if growth_context["growth_reaction_status"] != "multiple_growth_reactions_selected":
        return warnings
    selected = str(growth_context["growth_reaction_id"])
    open_reactions = ", ".join(str(item) for item in growth_context["open_growth_reaction_ids"])
    warnings.append(
        (
            "Multiple biomass reactions are open for this carbon-source boundary "
            f"({open_reactions}); fixed-growth probes select {selected} and close the others. "
            "This is a draft mixed-carbon convention, not a calibrated total-growth constraint."
        )
    )
    return warnings


def _open_growth_reactions(model: CobraModel) -> tuple[str, ...]:
    open_reactions: list[str] = []
    for reaction_id in _BIOMASS_REACTION_IDS:
        index = model.reaction_index.get(reaction_id)
        if index is None:
            continue
        if float(model.ub[index]) > 0.0:
            open_reactions.append(reaction_id)
    return tuple(open_reactions)


def _matching_carbon_source_formulation(
    model: CobraModel,
    open_reactions: tuple[str, ...],
) -> CarbonSourceFormulation | None:
    for formulation in list_carbon_source_formulations():
        if not all(_can_uptake(model, reaction_id) for reaction_id in formulation.active_uptake_reaction_ids):
            continue
        if any(_can_uptake(model, reaction_id) for reaction_id in formulation.blocked_uptake_reaction_ids):
            continue
        if not all(reaction_id in open_reactions for reaction_id in formulation.candidate_growth_reaction_ids):
            continue
        return formulation
    return None


def _select_growth_reaction(
    model: CobraModel,
    open_reactions: tuple[str, ...],
    formulation: CarbonSourceFormulation | None = None,
) -> str:
    if not open_reactions:
        return "BIOMASS"
    if formulation is not None and formulation.selected_growth_reaction_id in open_reactions:
        return formulation.selected_growth_reaction_id
    uptake = {
        "glucose": _can_uptake(model, "Ex_glc_D"),
        "glycerol": _can_uptake(model, "Ex_glyc"),
        "methanol": _can_uptake(model, "Ex_meoh"),
    }
    if uptake["glucose"] and "BIOMASS" in open_reactions:
        return "BIOMASS"
    if uptake["glycerol"] and "BIOMASS_glyc" in open_reactions:
        return "BIOMASS_glyc"
    if uptake["methanol"] and "BIOMASS_meoh" in open_reactions:
        return "BIOMASS_meoh"
    return open_reactions[0]


def _can_uptake(model: CobraModel, reaction_id: str) -> bool:
    index = model.reaction_index.get(reaction_id)
    return index is not None and float(model.lb[index]) < 0.0


def _uptake_cost(flux: float | None) -> float | None:
    if flux is None:
        return None
    if flux > 0:
        return None
    return -flux


def _uptake_cost_status(flux: float | None) -> str:
    if flux is None:
        return "unavailable"
    if flux > 0:
        return "non_uptake_flux"
    if flux < 0:
        return "uptake_flux"
    return "zero_flux"


def _normalize_cost_slope_medium_mode(value: str | None) -> str:
    mode = str(value or "corrected").strip().lower()
    aliases = {
        "": "corrected",
        "none": "corrected",
        "python_corrected": "corrected",
        "matlab": "matlab_legacy_cost",
        "matlab_legacy": "matlab_legacy_cost",
        "matlab_medium_compatibility": "matlab_legacy_cost",
    }
    mode = aliases.get(mode, mode)
    if mode not in {"corrected", "matlab_legacy_cost"}:
        raise ValueError(f"Unsupported cost slope medium compatibility mode: {value}")
    return mode


def _cost_slope_medium_bound_overrides(
    model: CobraModel,
    mode: str,
) -> tuple[dict[str, tuple[float, float]], list[dict[str, Any]], list[str]]:
    if mode == "corrected":
        return {}, [], [
            "Cost slope medium compatibility mode: corrected; no MATLAB legacy medium bounds are applied.",
        ]
    bounds: dict[str, tuple[float, float]] = {}
    overrides: list[dict[str, Any]] = []
    warnings = [
        "Cost slope medium compatibility mode: matlab_legacy_cost; applying 9 MATLAB legacy exchange lower bounds for historical protein-cost comparison only.",
    ]
    for reaction_id in MATLAB_LEGACY_COST_MEDIUM_EXCHANGES:
        index = model.reaction_index.get(reaction_id)
        if index is None:
            warnings.append(f"MATLAB legacy medium bound reaction not found: {reaction_id}.")
            continue
        current_lower = float(model.lb[index])
        current_upper = float(model.ub[index])
        bounds[reaction_id] = (0.0, current_upper)
        overrides.append(
            {
                "reaction_id": reaction_id,
                "legacy_lower_bound": 0.0,
                "corrected_lower_bound": current_lower,
                "upper_bound": current_upper,
                "reason": "MATLAB SimulateProteinCost legacy medium compatibility",
            }
        )
    return bounds, overrides, warnings


def _solve_weighted_carbon_uptake_minimize(
    model: CobraModel,
    *,
    carbon_weights: dict[str, float],
    metabolic: MetabolicEnzymeData,
    secretory: SecretoryEnzymeData,
    combined: CombinedEnzymeData,
    mu: float,
    key_reactions: Iterable[str] = (),
    total_protein_content: float = 0.37,
    unmodeled_er_protein_fraction: float = 0.040,
    mitochondrial_protein_fraction: float = 0.05,
    write_ribosome_translation_constraint: bool = False,
    write_misfolding_constraints: bool = False,
) -> tuple[SolveResult, dict[str, int]]:
    reaction_index = model.reaction_index
    carbon_reactions = tuple(carbon_weights)
    reaction_count = len(model.rxns)
    aux_count = len(carbon_reactions)
    objective = np.zeros(reaction_count + aux_count, dtype=float)
    for aux_offset, reaction_id in enumerate(carbon_reactions):
        objective[reaction_count + aux_offset] = float(carbon_weights[reaction_id])
    A_eq, b_eq, A_ub, b_ub, counts = build_pcsec_constraint_matrices(
        model,
        metabolic=metabolic,
        secretory=secretory,
        combined=combined,
        mu=mu,
        total_protein_content=total_protein_content,
        unmodeled_er_protein_fraction=unmodeled_er_protein_fraction,
        mitochondrial_protein_fraction=mitochondrial_protein_fraction,
        write_ribosome_translation_constraint=write_ribosome_translation_constraint,
        write_misfolding_constraints=write_misfolding_constraints,
    )
    A_eq_ext = sparse.hstack(
        [A_eq, sparse.csr_matrix((A_eq.shape[0], aux_count), dtype=float)],
        format="csr",
    )
    A_ub_ext = sparse.hstack(
        [A_ub, sparse.csr_matrix((A_ub.shape[0], aux_count), dtype=float)],
        format="csr",
    )
    uptake_rows = sparse.lil_matrix((aux_count, reaction_count + aux_count), dtype=float)
    for aux_offset, reaction_id in enumerate(carbon_reactions):
        uptake_rows[aux_offset, reaction_index[reaction_id]] = -1.0
        uptake_rows[aux_offset, reaction_count + aux_offset] = -1.0
    A_ub_ext = sparse.vstack([A_ub_ext, uptake_rows.tocsr()], format="csr")
    b_ub_ext = np.concatenate([np.asarray(b_ub, dtype=float), np.zeros(aux_count, dtype=float)])
    bounds = list(zip(model.lb.tolist(), model.ub.tolist()))
    bounds.extend((0.0, 1000.0) for _ in carbon_reactions)
    result = linprog(
        c=objective,
        A_eq=A_eq_ext,
        b_eq=b_eq,
        A_ub=A_ub_ext,
        b_ub=b_ub_ext,
        bounds=bounds,
        method="highs",
        options={"presolve": True, "disp": False},
    )
    fluxes: dict[str, float] = {}
    if result.success and result.x is not None:
        for reaction_id in {*carbon_weights.keys(), *key_reactions}:
            index = reaction_index.get(reaction_id)
            if index is not None:
                fluxes[reaction_id] = float(result.x[index])
        for aux_offset, reaction_id in enumerate(carbon_reactions):
            fluxes[f"uptake_cost::{reaction_id}"] = float(result.x[reaction_count + aux_offset])
    return (
        SolveResult(
            objective="weighted_carbon_uptake",
            status=str(result.status),
            success=bool(result.success),
            objective_value=float(result.fun) if result.success else None,
            message=str(result.message),
            fluxes=fluxes,
            sensitivity=None,
        ),
        counts,
    )


def _cost_slopes(rows: list[dict[str, Any]], cost_key: str) -> list[dict[str, Any]]:
    slopes: list[dict[str, Any]] = []
    mu_values = sorted({float(row["mu"]) for row in rows})
    for mu in mu_values:
        points = [
            (float(row["target_exchange_ratio"]), float(row[cost_key]))
            for row in rows
            if row.get("success") and row.get("mu") == mu and row.get(cost_key) is not None
        ]
        if len(points) < 2:
            slopes.append(
                {
                    "mu": mu,
                    "cost_key": cost_key,
                    "success": False,
                    "slope": None,
                    "point_count": len(points),
                    "status": "insufficient_points",
                }
            )
            continue
        slope = _linear_slope(points)
        slopes.append(
            {
                "mu": mu,
                "cost_key": cost_key,
                "success": slope is not None,
                "slope": slope,
                "point_count": len(points),
                "status": "slope_estimated" if slope is not None else "zero_variance",
            }
        )
    return slopes


def _linear_slope(points: list[tuple[float, float]]) -> float | None:
    x_values = [point[0] for point in points]
    y_values = [point[1] for point in points]
    x_mean = sum(x_values) / len(x_values)
    y_mean = sum(y_values) / len(y_values)
    denominator = sum((x_value - x_mean) ** 2 for x_value in x_values)
    if denominator <= 0:
        return None
    numerator = sum((x_value - x_mean) * (y_value - y_mean) for x_value, y_value in points)
    return float(numerator / denominator)


def _prepare_target_pcsec_inputs(
    model: CobraModel,
    target: TargetSpec,
    amino_acids: AminoAcidStoichiometry,
    secretory: SecretoryEnzymeData,
    combined: CombinedEnzymeData,
) -> dict[str, Any]:
    build = build_supported_target_model(model, target, amino_acids)
    if not build.supported or build.model is None or build.exchange_reaction_id is None:
        return {
            "model": None,
            "exchange_reaction_id": build.exchange_reaction_id,
            "secretory": None,
            "combined": None,
            "build_status": build.status,
            "reason": build.reason,
        }
    target_enzymedata = build_target_enzymedata(target, build.model, secretory)
    return {
        "model": build.model,
        "exchange_reaction_id": build.exchange_reaction_id,
        "secretory": secretory.with_reaction_coefficients(target_enzymedata.reaction_coefficients),
        "combined": combined.with_target(target_enzymedata),
        "target_enzymedata": target_enzymedata,
        "build_status": build.status,
        "reason": build.reason,
    }


def _target_parameter_status(target: TargetSpec) -> str:
    pending_targets = {"hlf"}
    if target.target_id.lower() in pending_targets or target.protein_id.lower() in pending_targets:
        return "draft_matlab_alignment_pending"
    return "draft"


__all__ = [
    "GrowthTradeoffResult",
    "MixedCarbonObjectiveResult",
    "ProteinCostSlopeCompatibilityResult",
    "MATLAB_LEGACY_COST_MEDIUM_EXCHANGES",
    "SecretionSimulationResult",
    "build_supported_target_model",
    "build_target_enzymedata",
    "run_growth_tradeoff",
    "run_mixed_carbon_objective_probe",
    "run_pcsec_growth_tradeoff",
    "run_protein_cost_slope_compatibility",
    "solve_maximize",
    "solve_pcsec_maximize",
    "solve_secretion_capacity",
    "summarize_mixed_carbon_objective_result",
    "summarize_simulation_result",
    "target_secretion_smoke",
]
