from __future__ import annotations

from dataclasses import dataclass

from pcsec_pichia.experimental_feedback.schema import (
    ExperimentBundle,
    ExperimentalFeedbackError,
)


@dataclass(frozen=True)
class ExperimentValidationIssue:
    code: str
    message: str
    record_type: str = "bundle"
    record_id: str = ""


@dataclass(frozen=True)
class ExperimentValidationResult:
    bundle: ExperimentBundle
    is_valid: bool
    errors: tuple[ExperimentValidationIssue, ...] = ()
    warnings: tuple[ExperimentValidationIssue, ...] = ()


def validate_experiment_bundle(bundle: ExperimentBundle) -> ExperimentValidationResult:
    try:
        bundle.validate()
    except ExperimentalFeedbackError as exc:
        return ExperimentValidationResult(
            bundle=bundle,
            is_valid=False,
            errors=(
                ExperimentValidationIssue(
                    code="schema_validation_error",
                    message=str(exc),
                ),
            ),
        )
    return ExperimentValidationResult(bundle=bundle, is_valid=True)


__all__ = [
    "ExperimentValidationIssue",
    "ExperimentValidationResult",
    "validate_experiment_bundle",
]
