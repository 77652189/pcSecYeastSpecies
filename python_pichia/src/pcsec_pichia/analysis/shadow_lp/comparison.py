from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Mapping

from pcsec_pichia.analysis.shadow_lp.backends import SolverBackend
from pcsec_pichia.analysis.shadow_lp.secretion_capacity import solve_shadow_secretion_capacity
from pcsec_pichia.probe import (
    AminoAcidStoichiometry,
    CobraModel,
    CombinedEnzymeData,
    MetabolicEnzymeData,
    SecretoryEnzymeData,
    TargetSpec,
)
from pcsec_pichia.simulation import SecretionSimulationResult, solve_secretion_capacity


@dataclass(frozen=True)
class SecretionCapacityComparisonResult:
    target_id: str
    growth_rate: float
    reference_result: SecretionSimulationResult
    shadow_result: SecretionSimulationResult
    objective_abs_diff: float | None
    objective_rel_diff: float | None
    constraint_count_diff: int | None
    key_flux_diffs: Mapping[str, Mapping[str, float | None]]
    reference_status_category: str
    shadow_status_category: str
    status_match: bool
    within_tolerance: bool
    tolerance: float
    warnings: tuple[str, ...]

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def compare_secretion_capacity(
    model: CobraModel,
    target: TargetSpec,
    amino_acids: AminoAcidStoichiometry,
    metabolic: MetabolicEnzymeData,
    secretory: SecretoryEnzymeData,
    combined: CombinedEnzymeData,
    growth_rate: float = 0.10,
    write_ribosome_translation_constraint: bool = False,
    write_misfolding_constraints: bool = False,
    backend: SolverBackend | None = None,
    solver_options: Mapping[str, Any] | None = None,
    relative_tolerance: float = 1e-4,
) -> SecretionCapacityComparisonResult:
    """Run reference and shadow capacity solvers on identical inputs and compare outputs."""

    reference = _run_reference(
        model,
        target,
        amino_acids,
        metabolic,
        secretory,
        combined,
        growth_rate=growth_rate,
        write_ribosome_translation_constraint=write_ribosome_translation_constraint,
        write_misfolding_constraints=write_misfolding_constraints,
    )
    shadow = solve_shadow_secretion_capacity(
        model,
        target,
        amino_acids,
        metabolic,
        secretory,
        combined,
        growth_rate=growth_rate,
        write_ribosome_translation_constraint=write_ribosome_translation_constraint,
        write_misfolding_constraints=write_misfolding_constraints,
        backend=backend,
        solver_options=solver_options,
    )
    reference_status = normalize_result_status(reference.success, reference.status, reference.message)
    shadow_status = normalize_result_status(shadow.success, shadow.status, shadow.message)
    abs_diff = _abs_diff(reference.objective_value, shadow.objective_value)
    rel_diff = _rel_diff(shadow.objective_value, reference.objective_value)
    constraint_diff = _constraint_total(shadow.constraint_counts) - _constraint_total(reference.constraint_counts)
    status_match = reference_status == shadow_status
    within = (
        status_match
        and reference.success
        and shadow.success
        and rel_diff is not None
        and rel_diff <= float(relative_tolerance)
        and constraint_diff == 0
    )
    return SecretionCapacityComparisonResult(
        target_id=target.target_id,
        growth_rate=float(growth_rate),
        reference_result=reference,
        shadow_result=shadow,
        objective_abs_diff=abs_diff,
        objective_rel_diff=rel_diff,
        constraint_count_diff=constraint_diff,
        key_flux_diffs=_key_flux_diffs(reference.key_fluxes or {}, shadow.key_fluxes or {}),
        reference_status_category=reference_status,
        shadow_status_category=shadow_status,
        status_match=status_match,
        within_tolerance=bool(within),
        tolerance=float(relative_tolerance),
        warnings=_comparison_warnings(reference, shadow, status_match, rel_diff, relative_tolerance),
    )


def normalize_result_status(success: bool, status: str | None, message: str | None = None) -> str:
    text = f"{status or ''} {message or ''}".strip().lower()
    if success or text in {"0", "optimal"} or "optimal" in text:
        return "optimal"
    if "infeasible" in text:
        return "infeasible"
    if "unbounded" in text:
        return "unbounded"
    if "time" in text or "iteration" in text or "limit" in text:
        return "timeout_iteration_limit"
    if "unavailable" in text or "too_large" in text or "not available" in text:
        return "unavailable_backend"
    if "numerical" in text or "precision" in text:
        return "numerical_failure"
    if "exception" in text:
        return "exception"
    return "exception"


def _run_reference(
    model: CobraModel,
    target: TargetSpec,
    amino_acids: AminoAcidStoichiometry,
    metabolic: MetabolicEnzymeData,
    secretory: SecretoryEnzymeData,
    combined: CombinedEnzymeData,
    *,
    growth_rate: float,
    write_ribosome_translation_constraint: bool,
    write_misfolding_constraints: bool,
) -> SecretionSimulationResult:
    try:
        return solve_secretion_capacity(
            model,
            target,
            amino_acids,
            metabolic,
            secretory,
            combined,
            growth_rate=growth_rate,
            write_ribosome_translation_constraint=write_ribosome_translation_constraint,
            write_misfolding_constraints=write_misfolding_constraints,
        )
    except Exception as exc:
        return SecretionSimulationResult(
            success=False,
            target_id=target.target_id,
            objective_value=None,
            growth_rate=float(growth_rate),
            secretion_flux=None,
            status="exception",
            message=str(exc),
            constraint_counts={},
            result_status="reference_exception",
            target_parameter_status="unknown",
            matlab_alignment_status="pending",
            exchange_reaction_id=None,
            build_status="exception",
            lp_sensitivity=None,
            key_fluxes={},
            warnings=(str(exc),),
            solver_mode="reference",
        )


def _constraint_total(counts: Mapping[str, int]) -> int:
    if "eq_total" in counts or "ub_total" in counts:
        return int(counts.get("eq_total", 0)) + int(counts.get("ub_total", 0))
    return int(sum(int(value) for value in counts.values()))


def _key_flux_diffs(
    reference_fluxes: Mapping[str, float],
    shadow_fluxes: Mapping[str, float],
) -> dict[str, dict[str, float | None]]:
    rows: dict[str, dict[str, float | None]] = {}
    for reaction_id in sorted(set(reference_fluxes) | set(shadow_fluxes)):
        reference_value = _optional_float(reference_fluxes.get(reaction_id))
        shadow_value = _optional_float(shadow_fluxes.get(reaction_id))
        rows[reaction_id] = {
            "reference": reference_value,
            "shadow": shadow_value,
            "abs_diff": _abs_diff(reference_value, shadow_value),
        }
    return rows


def _comparison_warnings(
    reference: SecretionSimulationResult,
    shadow: SecretionSimulationResult,
    status_match: bool,
    rel_diff: float | None,
    tolerance: float,
) -> tuple[str, ...]:
    warnings: list[str] = []
    warnings.extend(str(warning) for warning in reference.warnings)
    warnings.extend(str(warning) for warning in shadow.warnings)
    if not status_match:
        warnings.append("Reference and shadow solver status categories differ.")
    if rel_diff is None:
        warnings.append("Objective relative diff is unavailable.")
    elif rel_diff > float(tolerance):
        warnings.append("Objective relative diff exceeds tolerance.")
    return tuple(dict.fromkeys(warnings))


def _abs_diff(left: float | None, right: float | None) -> float | None:
    if left is None or right is None:
        return None
    return abs(float(left) - float(right))


def _rel_diff(left: float | None, right: float | None) -> float | None:
    if left is None or right in (None, 0):
        return None
    return abs(float(left) - float(right)) / abs(float(right))


def _optional_float(value: Any) -> float | None:
    if value is None:
        return None
    return float(value)


__all__ = [
    "SecretionCapacityComparisonResult",
    "compare_secretion_capacity",
    "normalize_result_status",
]
