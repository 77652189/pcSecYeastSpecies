from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Mapping

from pcsec_pichia.analysis.shadow_lp.comparison import SecretionCapacityComparisonResult, compare_secretion_capacity
from pcsec_pichia.analysis.shadow_lp.constraint_spec import ShadowConstraintConfig
from pcsec_pichia.analysis.shadow_lp.ladder import ShadowLadderResult
from pcsec_pichia.loading import load_pcsec_pichia_inputs
from pcsec_pichia.simulation import solve_secretion_capacity
from pcsec_pichia.targets import load_builtin_targets


@dataclass(frozen=True)
class ReferenceValidationResult:
    target_id: str
    reference_objective: float | None
    shadow_objective: float | None
    objective_abs_diff: float | None
    objective_rel_diff: float | None
    constraint_count_diff: int | None
    final_alignment_status: str
    reference_status: str
    reference_message: str
    reference_constraint_counts: Mapping[str, int]
    reference_exchange_reaction_id: str | None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class ShadowValidationMatrixCase:
    target_id: str
    growth_rate: float
    reference_status: str
    shadow_status: str
    objective_abs_diff: float | None
    objective_rel_diff: float | None
    constraint_count_diff: int | None
    key_flux_diffs: Mapping[str, Mapping[str, float | None]]
    within_tolerance: bool
    alignment_status: str
    warnings: tuple[str, ...]

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class ShadowValidationMatrixResult:
    target_ids: tuple[str, ...]
    growth_rates: tuple[float, ...]
    tolerance: float
    cases: tuple[ShadowValidationMatrixCase, ...]
    replacement_evidence_cases: tuple[ShadowValidationMatrixCase, ...]
    review_required_cases: tuple[ShadowValidationMatrixCase, ...]

    @property
    def all_required_defaults_aligned(self) -> bool:
        required = {
            (target_id, 0.10)
            for target_id in self.target_ids
            if target_id in {"hLF", "OPN_ALPHA_FULL_PROJECT"}
        }
        aligned = {
            (case.target_id, case.growth_rate)
            for case in self.cases
            if case.within_tolerance
        }
        return required.issubset(aligned)

    def to_dict(self) -> dict[str, Any]:
        payload = asdict(self)
        payload["all_required_defaults_aligned"] = self.all_required_defaults_aligned
        return payload


def solve_pcsec_reference_for_validation(
    target_id: str,
    root: Path | None = None,
    config: ShadowConstraintConfig | None = None,
) -> Any:
    """Run the pcSec reference path for validation-only comparison."""

    resolved_config = config or ShadowConstraintConfig()
    inputs = load_pcsec_pichia_inputs(root)
    targets = {target.target_id: target for target in load_builtin_targets(inputs.root)}
    try:
        target = targets[target_id]
    except KeyError as exc:
        raise KeyError(f"Unknown built-in target: {target_id}") from exc
    return solve_secretion_capacity(
        inputs.prepared_model,
        target,
        inputs.amino_acids,
        inputs.metabolic,
        inputs.secretory,
        inputs.combined,
        growth_rate=resolved_config.growth_rate,
        write_ribosome_translation_constraint=False,
        write_misfolding_constraints=False,
    )


def validate_shadow_ladder_against_reference(
    ladder: ShadowLadderResult,
    root: Path | None = None,
    config: ShadowConstraintConfig | None = None,
    relative_tolerance: float = 1e-4,
) -> ReferenceValidationResult:
    """Compare an already-solved shadow ladder with the validation-only pcSec reference."""

    reference = solve_pcsec_reference_for_validation(ladder.target_id, root=root, config=config)
    final = ladder.final_layer
    shadow_objective = final.objective
    reference_objective = reference.objective_value
    abs_diff = _abs_diff(shadow_objective, reference_objective)
    rel_diff = _rel_diff(shadow_objective, reference_objective)
    reference_constraint_total = _reference_constraint_total(reference.constraint_counts)
    constraint_count_diff = (
        final.constraint_count - reference_constraint_total
        if reference_constraint_total is not None
        else None
    )
    aligned = bool(reference.success) and rel_diff is not None and rel_diff <= relative_tolerance
    return ReferenceValidationResult(
        target_id=ladder.target_id,
        reference_objective=reference_objective,
        shadow_objective=shadow_objective,
        objective_abs_diff=abs_diff,
        objective_rel_diff=rel_diff,
        constraint_count_diff=constraint_count_diff,
        final_alignment_status="aligned" if aligned else "review_required",
        reference_status=str(reference.status),
        reference_message=str(reference.message),
        reference_constraint_counts={str(key): int(value) for key, value in reference.constraint_counts.items()},
        reference_exchange_reaction_id=reference.exchange_reaction_id,
    )


def attach_reference_validation(
    ladder: ShadowLadderResult,
    root: Path | None = None,
    config: ShadowConstraintConfig | None = None,
    relative_tolerance: float = 1e-4,
) -> ShadowLadderResult:
    validation = validate_shadow_ladder_against_reference(
        ladder,
        root=root,
        config=config,
        relative_tolerance=relative_tolerance,
    )
    return ladder.with_reference_validation(validation.to_dict())


def run_shadow_validation_matrix(
    root: Path | None = None,
    target_ids: tuple[str, ...] = ("hLF", "OPN_ALPHA_FULL_PROJECT"),
    growth_rates: tuple[float, ...] = (0.05, 0.10, 0.15),
    relative_tolerance: float = 1e-4,
) -> ShadowValidationMatrixResult:
    """Compare reference and shadow capacity across built-in targets and growth rates."""

    inputs = load_pcsec_pichia_inputs(root)
    targets = {target.target_id: target for target in load_builtin_targets(inputs.root)}
    cases: list[ShadowValidationMatrixCase] = []
    for target_id in target_ids:
        target = targets[target_id]
        for growth_rate in growth_rates:
            cases.append(
                _matrix_case_from_compare(
                    compare_secretion_capacity(
                        inputs.prepared_model,
                        target,
                        inputs.amino_acids,
                        inputs.metabolic,
                        inputs.secretory,
                        inputs.combined,
                        growth_rate=float(growth_rate),
                        relative_tolerance=relative_tolerance,
                    )
                )
            )
    replacement = tuple(case for case in cases if case.within_tolerance)
    review_required = tuple(case for case in cases if not case.within_tolerance)
    return ShadowValidationMatrixResult(
        target_ids=tuple(target_ids),
        growth_rates=tuple(float(value) for value in growth_rates),
        tolerance=float(relative_tolerance),
        cases=tuple(cases),
        replacement_evidence_cases=replacement,
        review_required_cases=review_required,
    )


def _reference_constraint_total(counts: Mapping[str, int]) -> int | None:
    if "eq_total" in counts or "ub_total" in counts:
        return int(counts.get("eq_total", 0)) + int(counts.get("ub_total", 0))
    return None


def _matrix_case_from_compare(
    comparison: SecretionCapacityComparisonResult,
) -> ShadowValidationMatrixCase:
    alignment_status = "aligned" if comparison.within_tolerance else "review_required"
    return ShadowValidationMatrixCase(
        target_id=comparison.target_id,
        growth_rate=comparison.growth_rate,
        reference_status=comparison.reference_status_category,
        shadow_status=comparison.shadow_status_category,
        objective_abs_diff=comparison.objective_abs_diff,
        objective_rel_diff=comparison.objective_rel_diff,
        constraint_count_diff=comparison.constraint_count_diff,
        key_flux_diffs=comparison.key_flux_diffs,
        within_tolerance=comparison.within_tolerance,
        alignment_status=alignment_status,
        warnings=comparison.warnings,
    )


def _abs_diff(left: float | None, right: float | None) -> float | None:
    if left is None or right is None:
        return None
    return abs(float(left) - float(right))


def _rel_diff(left: float | None, right: float | None) -> float | None:
    if left is None or right in (None, 0):
        return None
    return abs(float(left) - float(right)) / abs(float(right))
