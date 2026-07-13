from __future__ import annotations

from pcsec_pichia.experimental_feedback.calibration import build_calibration_summary
from pcsec_pichia.experimental_feedback.io import (
    load_experiment_bundle,
    write_experiment_feedback_cache,
)
from pcsec_pichia.experimental_feedback.linkage import (
    build_prediction_index,
    link_experiments_to_predictions,
)
from pcsec_pichia.experimental_feedback.quality import validate_experiment_bundle
from pcsec_pichia.experimental_feedback.schema import ExperimentalFeedbackError


class ExperimentalFeedbackPhaseError(ExperimentalFeedbackError):
    """Raised when a frozen public entrypoint belongs to a later Phase 1 round."""


__all__ = [
    "ExperimentalFeedbackPhaseError",
    "build_calibration_summary",
    "build_prediction_index",
    "link_experiments_to_predictions",
    "load_experiment_bundle",
    "validate_experiment_bundle",
    "write_experiment_feedback_cache",
]
