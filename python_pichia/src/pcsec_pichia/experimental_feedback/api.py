from __future__ import annotations

from typing import Any, Iterable

from pcsec_pichia.experimental_feedback.io import (
    load_experiment_bundle,
    write_experiment_feedback_cache,
)
from pcsec_pichia.experimental_feedback.quality import (
    ExperimentValidationResult,
    validate_experiment_bundle,
)
from pcsec_pichia.experimental_feedback.schema import ExperimentBundle, ExperimentalFeedbackError


class ExperimentalFeedbackPhaseError(ExperimentalFeedbackError):
    """Raised when a frozen public entrypoint belongs to a later Phase 1 round."""


def build_prediction_index(screen_runs: Iterable[Any]) -> Any:
    raise ExperimentalFeedbackPhaseError("build_prediction_index is implemented in Phase 1 Round 2.")


def link_experiments_to_predictions(bundle: ExperimentBundle, prediction_index: Any) -> Any:
    raise ExperimentalFeedbackPhaseError(
        "link_experiments_to_predictions is implemented in Phase 1 Round 2."
    )


def build_calibration_summary(
    validated_bundle: ExperimentBundle | ExperimentValidationResult,
    linkage_result: Any,
    config: Any,
) -> Any:
    raise ExperimentalFeedbackPhaseError("build_calibration_summary is implemented in Phase 1 Round 3.")


__all__ = [
    "ExperimentalFeedbackPhaseError",
    "build_calibration_summary",
    "build_prediction_index",
    "link_experiments_to_predictions",
    "load_experiment_bundle",
    "validate_experiment_bundle",
    "write_experiment_feedback_cache",
]
