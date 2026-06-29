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


@dataclass(frozen=True)
class GrowthTradeoffResult:
    target_id: str
    success: bool
    growth_points: tuple[float, ...]
    tradeoff_rows: tuple[dict[str, Any], ...]
    result_status: str
    target_parameter_status: str
    matlab_alignment_status: str


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
    "SecretionSimulationResult",
    "build_supported_target_model",
    "build_target_enzymedata",
    "run_growth_tradeoff",
    "run_pcsec_growth_tradeoff",
    "solve_maximize",
    "solve_pcsec_maximize",
    "solve_secretion_capacity",
    "summarize_simulation_result",
    "target_secretion_smoke",
]
