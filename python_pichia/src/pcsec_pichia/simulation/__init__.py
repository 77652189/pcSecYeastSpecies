from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable

from pcsec_pichia.probe import (
    AminoAcidStoichiometry,
    CobraModel,
    CombinedEnzymeData,
    MetabolicEnzymeData,
    SecretoryEnzymeData,
    TargetSpec,
    build_supported_target_model,
    build_target_enzymedata,
    run_pcsec_growth_tradeoff,
    solve_maximize,
    solve_pcsec_maximize,
    target_secretion_smoke,
)
from pcsec_pichia.media import MATLAB_LEGACY_COST_MEDIUM_EXCHANGES


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


@dataclass(frozen=True)
class GrowthTradeoffResult:
    target_id: str
    success: bool
    growth_points: tuple[float, ...]
    tradeoff_rows: tuple[dict[str, Any], ...]
    result_status: str
    target_parameter_status: str
    matlab_alignment_status: str

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

    fixed_model = prepared["model"].with_bounds({"BIOMASS": (growth_rate, growth_rate)})
    solved, counts = solve_pcsec_maximize(
        fixed_model,
        prepared["exchange_reaction_id"],
        metabolic=metabolic,
        secretory=prepared["secretory"],
        combined=prepared["combined"],
        mu=growth_rate,
        key_reactions=("BIOMASS", "Ex_glc_D", "Ex_o2", prepared["exchange_reaction_id"]),
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

    rows = run_pcsec_growth_tradeoff(
        prepared["model"],
        prepared["exchange_reaction_id"],
        metabolic=metabolic,
        secretory=prepared["secretory"],
        combined=prepared["combined"],
        mu_points=normalized_points,
        write_ribosome_translation_constraint=write_ribosome_translation_constraint,
        write_misfolding_constraints=write_misfolding_constraints,
    )
    return GrowthTradeoffResult(
        target_id=target.target_id,
        success=bool(rows) and all(bool(row.get("success")) for row in rows),
        growth_points=normalized_points,
        tradeoff_rows=tuple(rows),
        result_status="draft",
        target_parameter_status=_target_parameter_status(target),
        matlab_alignment_status="pending",
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
    """Run optional MATLAB-style protein cost slope probes."""

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
    rows: list[dict[str, Any]] = []
    for mu in normalized_mu:
        for ratio in normalized_ratios:
            fixed_model = prepared["model"].with_bounds(
                {
                    "BIOMASS": (mu, mu),
                    exchange_reaction_id: (ratio, ratio),
                    **medium_bounds,
                }
            )
            solved, counts = solve_pcsec_maximize(
                fixed_model,
                "Ex_glc_D",
                metabolic=metabolic,
                secretory=prepared["secretory"],
                combined=prepared["combined"],
                mu=mu,
                key_reactions=("BIOMASS", "Ex_glc_D", "Ex_o2", exchange_reaction_id, ribosome_reaction_id),
                write_ribosome_translation_constraint=write_ribosome_translation_constraint,
                write_misfolding_constraints=write_misfolding_constraints,
            )
            glucose_flux = _optional_flux(solved.fluxes.get("Ex_glc_D"))
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
                    "glucose_cost": _uptake_cost(glucose_flux),
                    "glucose_cost_status": _uptake_cost_status(glucose_flux),
                    "ribosome_reaction": ribosome_reaction_id,
                    "ribosome_flux": ribosome_flux,
                    "ribosome_cost": abs(ribosome_flux) if ribosome_flux is not None else None,
                    "constraint_counts": counts,
                    "medium_compatibility_mode": normalized_medium_mode,
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
        }
    return {
        "success": result.success,
        "target_id": result.target_id,
        "growth_points": result.growth_points,
        "tradeoff_rows": result.tradeoff_rows,
        "result_status": result.result_status,
        "target_parameter_status": result.target_parameter_status,
        "matlab_alignment_status": result.matlab_alignment_status,
    }


def _optional_flux(value: object) -> float | None:
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


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
    "ProteinCostSlopeCompatibilityResult",
    "MATLAB_LEGACY_COST_MEDIUM_EXCHANGES",
    "SecretionSimulationResult",
    "build_supported_target_model",
    "build_target_enzymedata",
    "run_growth_tradeoff",
    "run_protein_cost_slope_compatibility",
    "run_pcsec_growth_tradeoff",
    "solve_maximize",
    "solve_pcsec_maximize",
    "solve_secretion_capacity",
    "summarize_simulation_result",
    "target_secretion_smoke",
]
