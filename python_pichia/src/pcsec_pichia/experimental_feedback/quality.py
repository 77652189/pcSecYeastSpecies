from __future__ import annotations

from dataclasses import dataclass

from pcsec_pichia.experimental_feedback.schema import (
    ExperimentBundle,
    ExperimentalFeedbackError,
    UnitValidationError,
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
    errors = [
        ExperimentValidationIssue(
            code=conflict.code,
            message=f"{conflict.record_type} {conflict.record_id} has conflicting payloads.",
            record_type=conflict.record_type,
            record_id=conflict.record_id,
        )
        for conflict in bundle.import_conflicts
    ]
    warnings = [
        ExperimentValidationIssue(code="import_warning", message=warning)
        for warning in bundle.warnings
    ]
    for experiment in bundle.experiments:
        condition = experiment.condition
        missing_fields = [
            name
            for name, value in (
                ("condition_description", condition.condition_description),
                ("sampling_time_h", condition.sampling_time_h),
            )
            if value is None or str(value).strip().lower() in {"", "missing", "unknown"}
        ]
        if missing_fields:
            warnings.append(
                ExperimentValidationIssue(
                    code="condition_missing",
                    message=f"missing condition fields: {', '.join(missing_fields)}",
                    record_type="experiment",
                    record_id=experiment.experiment_id,
                )
            )
    for record_type, records in (
        ("experiment", bundle.experiments),
        ("intervention", bundle.interventions),
        ("measurement", bundle.measurements),
        ("prediction_link", bundle.prediction_links),
    ):
        for record in records:
            try:
                record.validate()
            except ExperimentalFeedbackError as exc:
                errors.append(
                    ExperimentValidationIssue(
                        code=_issue_code(exc),
                        message=str(exc),
                        record_type=record_type,
                        record_id=_record_id(record_type, record),
                    )
                )
    experiment_ids = {record.experiment_id for record in bundle.experiments}
    intervention_keys = {
        (record.experiment_id, record.intervention_id) for record in bundle.interventions
    }
    for record_type, records in (
        ("intervention", bundle.interventions),
        ("measurement", bundle.measurements),
    ):
        for record in records:
            if record.experiment_id not in experiment_ids:
                errors.append(
                    ExperimentValidationIssue(
                        code="missing_experiment_reference",
                        message=f"{record_type} references missing experiment_id: {record.experiment_id}",
                        record_type=record_type,
                        record_id=_record_id(record_type, record),
                    )
                )
    for record in bundle.prediction_links:
        if (record.experiment_id, record.intervention_id) not in intervention_keys:
            errors.append(
                ExperimentValidationIssue(
                    code="missing_intervention_reference",
                    message=(
                        "prediction_link references missing intervention: "
                        f"{record.experiment_id}/{record.intervention_id}"
                    ),
                    record_type="prediction_link",
                    record_id=_record_id("prediction_link", record),
                )
            )
    try:
        bundle.validate()
    except ExperimentalFeedbackError as exc:
        if not any(issue.message == str(exc) for issue in errors):
            errors.append(
                ExperimentValidationIssue(
                    code=_issue_code(exc),
                    message=str(exc),
                )
            )
    errors = _dedupe_issues(errors)
    warnings = _dedupe_issues(warnings)
    return ExperimentValidationResult(
        bundle=bundle,
        is_valid=not errors,
        errors=tuple(errors),
        warnings=tuple(warnings),
    )


def _issue_code(exc: ExperimentalFeedbackError) -> str:
    return "unit_validation_error" if isinstance(exc, UnitValidationError) else "schema_validation_error"


def _record_id(record_type: str, record: object) -> str:
    if record_type == "experiment":
        return str(record.experiment_id)  # type: ignore[attr-defined]
    if record_type == "intervention":
        return f"{record.experiment_id}/{record.intervention_id}"  # type: ignore[attr-defined]
    if record_type == "measurement":
        return f"{record.experiment_id}/{record.measurement_id}"  # type: ignore[attr-defined]
    return f"{record.experiment_id}/{record.intervention_id}"  # type: ignore[attr-defined]


def _dedupe_issues(
    issues: list[ExperimentValidationIssue],
) -> list[ExperimentValidationIssue]:
    deduped: list[ExperimentValidationIssue] = []
    seen: set[tuple[str, str, str, str]] = set()
    for issue in issues:
        key = (issue.code, issue.message, issue.record_type, issue.record_id)
        if key not in seen:
            deduped.append(issue)
            seen.add(key)
    return deduped


__all__ = [
    "ExperimentValidationIssue",
    "ExperimentValidationResult",
    "validate_experiment_bundle",
]
