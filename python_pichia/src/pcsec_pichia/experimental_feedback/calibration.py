from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from enum import Enum
from pathlib import Path
from statistics import fmean
from typing import Mapping

from pcsec_pichia.experimental_feedback.linkage import PredictionLinkageResult
from pcsec_pichia.experimental_feedback.quality import (
    ExperimentValidationResult,
    validate_experiment_bundle,
)
from pcsec_pichia.experimental_feedback.schema import (
    ExperimentBundle,
    InterventionType,
    MeasurementRecord,
    MeasurementStatus,
    PredictionLinkRecord,
    PredictionLinkStatus,
    QualityStatus,
    SchemaValidationError,
)


@dataclass(frozen=True)
class CalibrationConfig:
    primary_assay_type: str = "titer"
    increase_threshold_ratio: float = 1.05
    decrease_threshold_ratio: float = 0.95
    top_k: tuple[int, ...] = (5, 10)
    baseline_hit_rate: float = 0.10
    minimum_rank_pairs: int = 3
    config_version: int = 1

    def validate(self) -> None:
        if not self.primary_assay_type:
            raise SchemaValidationError("primary_assay_type must be non-empty.")
        numeric_values = (
            self.increase_threshold_ratio,
            self.decrease_threshold_ratio,
            self.baseline_hit_rate,
        )
        if any(
            not isinstance(value, (int, float))
            or isinstance(value, bool)
            or not math.isfinite(float(value))
            for value in numeric_values
        ):
            raise SchemaValidationError("calibration thresholds and baseline must be finite.")
        if self.increase_threshold_ratio <= 1.0:
            raise SchemaValidationError("increase_threshold_ratio must be greater than 1.")
        if not 0 < self.decrease_threshold_ratio < 1.0:
            raise SchemaValidationError("decrease_threshold_ratio must be between 0 and 1.")
        if not self.top_k or any(
            not isinstance(value, int) or isinstance(value, bool) or value < 1
            for value in self.top_k
        ):
            raise SchemaValidationError("top_k must contain positive integers.")
        if len(set(self.top_k)) != len(self.top_k):
            raise SchemaValidationError("top_k values must be unique.")
        if not 0 < self.baseline_hit_rate <= 1.0:
            raise SchemaValidationError("baseline_hit_rate must be in (0, 1].")
        if not isinstance(self.minimum_rank_pairs, int) or isinstance(self.minimum_rank_pairs, bool):
            raise SchemaValidationError("minimum_rank_pairs must be an integer.")
        if self.minimum_rank_pairs < 2:
            raise SchemaValidationError("minimum_rank_pairs must be at least 2.")
        if self.config_version != 1:
            raise SchemaValidationError("config_version must be 1.")


@dataclass(frozen=True)
class CalibrationRecord:
    experiment_id: str
    intervention_id: str
    target_id: str
    evidence_id: str
    prediction_run_id: str
    prediction_rank: int | None
    evidence_tier: str
    recommendation_tier: str
    predicted_direction: str
    measurement_ids: tuple[str, ...]
    measurement_statuses: tuple[str, ...]
    control_experiment_ids: tuple[str, ...]
    control_measurement_ids: tuple[str, ...]
    eligibility_status: str
    ineligibility_reasons: tuple[str, ...]
    candidate_value: float | None = None
    control_value: float | None = None
    observed_ratio: float | None = None
    observed_direction: str = ""
    direction_consistent: bool | None = None
    hit: bool | None = None


@dataclass(frozen=True)
class TopKCalibrationMetric:
    k: int
    tested_count: int
    hit_count: int
    hit_rate: float | None
    relative_baseline_enrichment: float | None


@dataclass(frozen=True)
class EvidenceTierCalibrationMetric:
    evidence_tier: str
    tested_count: int
    hit_count: int
    hit_rate: float | None


@dataclass(frozen=True)
class TargetCalibrationSummary:
    target_id: str
    eligible_count: int
    ineligible_count: int
    direction_consistency_rate: float | None
    rank_correlation: float | None
    top_k_metrics: tuple[TopKCalibrationMetric, ...]
    evidence_tier_metrics: tuple[EvidenceTierCalibrationMetric, ...]


@dataclass(frozen=True)
class CalibrationSummary:
    config: CalibrationConfig
    records: tuple[CalibrationRecord, ...]
    targets: tuple[TargetCalibrationSummary, ...]


@dataclass(frozen=True)
class CalibrationOutputs:
    records_path: Path
    summary_path: Path
    manifest_path: Path


def build_calibration_summary(
    validated_bundle: ExperimentBundle | ExperimentValidationResult,
    linkage_result: PredictionLinkageResult,
    config: CalibrationConfig,
) -> CalibrationSummary:
    config.validate()
    if isinstance(validated_bundle, ExperimentValidationResult):
        if not validated_bundle.is_valid:
            raise SchemaValidationError("calibration requires a valid experiment bundle.")
        bundle = validated_bundle.bundle
    else:
        validation = validate_experiment_bundle(validated_bundle)
        if not validation.is_valid:
            raise SchemaValidationError("calibration requires a valid experiment bundle.")
        bundle = validation.bundle
    records = _build_records(bundle, linkage_result, config)
    return CalibrationSummary(
        config=config,
        records=records,
        targets=_build_target_summaries(records, config),
    )


def write_calibration_outputs(
    summary: CalibrationSummary,
    output_dir: str | Path,
) -> CalibrationOutputs:
    resolved = Path(output_dir)
    resolved.mkdir(parents=True, exist_ok=True)
    outputs = CalibrationOutputs(
        records_path=resolved / "calibration_records.jsonl",
        summary_path=resolved / "calibration_summary.json",
        manifest_path=resolved / "calibration_manifest.json",
    )
    with outputs.records_path.open("w", encoding="utf-8", newline="\n") as handle:
        for record in summary.records:
            handle.write(json.dumps(_json_ready(asdict(record)), ensure_ascii=False, sort_keys=True))
            handle.write("\n")
    outputs.summary_path.write_text(
        json.dumps(
            {
                "config": _json_ready(asdict(summary.config)),
                "targets": _json_ready([asdict(target) for target in summary.targets]),
                "record_count": len(summary.records),
                "eligible_count": sum(
                    record.eligibility_status == "eligible" for record in summary.records
                ),
                "ineligible_count": sum(
                    record.eligibility_status != "eligible" for record in summary.records
                ),
            },
            ensure_ascii=False,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    outputs.manifest_path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "generated_at": datetime.now(timezone.utc).isoformat(),
                "config": _json_ready(asdict(summary.config)),
                "record_count": len(summary.records),
                "target_count": len(summary.targets),
                "files": {
                    "records": outputs.records_path.name,
                    "summary": outputs.summary_path.name,
                },
                "mutates_recommendation_tier": False,
                "mutates_model_constraints": False,
                "top_k_denominator": "eligible tested predictions with prediction_rank <= k",
            },
            ensure_ascii=False,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    return outputs


def _build_records(
    bundle: ExperimentBundle,
    linkage_result: PredictionLinkageResult,
    config: CalibrationConfig,
) -> tuple[CalibrationRecord, ...]:
    experiments = {record.experiment_id: record for record in bundle.experiments}
    interventions_by_experiment: dict[str, list[object]] = {}
    for intervention in bundle.interventions:
        interventions_by_experiment.setdefault(intervention.experiment_id, []).append(intervention)
    measurements_by_experiment: dict[str, list[MeasurementRecord]] = {}
    for measurement in bundle.measurements:
        measurements_by_experiment.setdefault(measurement.experiment_id, []).append(measurement)
    controls = {
        experiment_id
        for experiment_id, interventions in interventions_by_experiment.items()
        if interventions and all(
            intervention.intervention_type is InterventionType.CONTROL for intervention in interventions
        )
    }
    links_by_intervention: dict[tuple[str, str], list[PredictionLinkRecord]] = {}
    for link in linkage_result.links:
        links_by_intervention.setdefault((link.experiment_id, link.intervention_id), []).append(link)
    records: list[CalibrationRecord] = []
    for experiment_id, interventions in interventions_by_experiment.items():
        non_controls = [
            intervention
            for intervention in interventions
            if intervention.intervention_type is not InterventionType.CONTROL
        ]
        for intervention in non_controls:
            experiment = experiments.get(experiment_id)
            if experiment is None:
                continue
            links = links_by_intervention.get((experiment_id, intervention.intervention_id), [])
            measurements = [
                measurement
                for measurement in measurements_by_experiment.get(experiment_id, [])
                if measurement.assay_type == config.primary_assay_type
            ]
            base = _record_base(experiment, intervention, links, measurements)
            if experiment.quality_status is not QualityStatus.VALID:
                records.append(
                    _ineligible(base, f"experiment_quality_status:{experiment.quality_status.value}")
                )
                continue
            if not _experiment_context_complete(experiment):
                records.append(_ineligible(base, "experiment_context_incomplete"))
                continue
            if len(non_controls) != 1:
                records.append(_ineligible(base, "combination_intervention_not_attributable"))
                continue
            if len(links) != 1 or links[0].status is not PredictionLinkStatus.MATCHED:
                reason = links[0].status.value if links else "prediction_link_missing"
                records.append(_ineligible(base, f"prediction_link:{reason}"))
                continue
            valid_candidate = _valid_measurements(measurements)
            if not valid_candidate:
                records.append(_ineligible(base, "candidate_measurement_not_evaluable"))
                continue
            signatures = {_measurement_signature(item) for item in valid_candidate}
            if len(signatures) != 1:
                records.append(_ineligible(base, "candidate_measurement_context_ambiguous"))
                continue
            signature = next(iter(signatures))
            control_experiments = [
                control_id
                for control_id in controls
                if _control_context_matches(experiment, experiments[control_id])
                and experiments[control_id].quality_status is QualityStatus.VALID
                and _experiment_context_complete(experiments[control_id])
            ]
            control_groups = {
                control_id: [
                    measurement
                    for measurement in _valid_measurements(
                        measurements_by_experiment.get(control_id, [])
                    )
                    if _measurement_signature(measurement) == signature
                ]
                for control_id in control_experiments
            }
            control_groups = {
                control_id: measurements
                for control_id, measurements in control_groups.items()
                if measurements
            }
            if not control_groups:
                records.append(_ineligible(base, "control_match_missing"))
                continue
            candidate_value = fmean(float(item.canonical_value) for item in valid_candidate)
            control_values = [
                fmean(float(item.canonical_value) for item in measurements)
                for measurements in control_groups.values()
            ]
            control_value = fmean(control_values)
            if control_value == 0:
                records.append(_ineligible(base, "control_value_zero"))
                continue
            ratio = candidate_value / control_value
            observed_direction = _observed_direction(ratio, config)
            link = links[0]
            records.append(
                CalibrationRecord(
                    **base,
                    control_experiment_ids=tuple(sorted(control_groups)),
                    control_measurement_ids=tuple(
                        f"{control_id}/{item.measurement_id}"
                        for control_id in sorted(control_groups)
                        for item in control_groups[control_id]
                    ),
                    eligibility_status="eligible",
                    ineligibility_reasons=(),
                    candidate_value=candidate_value,
                    control_value=control_value,
                    observed_ratio=ratio,
                    observed_direction=observed_direction,
                    direction_consistent=(
                        observed_direction == link.predicted_direction
                        if link.predicted_direction
                        else None
                    ),
                    hit=ratio >= config.increase_threshold_ratio,
                )
            )
    return tuple(records)


def _record_base(
    experiment: object,
    intervention: object,
    links: list[PredictionLinkRecord],
    measurements: list[MeasurementRecord],
) -> dict[str, object]:
    link = links[0] if len(links) == 1 else None
    return {
        "experiment_id": experiment.experiment_id,  # type: ignore[attr-defined]
        "intervention_id": intervention.intervention_id,  # type: ignore[attr-defined]
        "target_id": experiment.target_id,  # type: ignore[attr-defined]
        "evidence_id": link.evidence_id if link else "",
        "prediction_run_id": link.prediction_run_id if link else "",
        "prediction_rank": link.prediction_rank if link else None,
        "evidence_tier": link.evidence_tier if link else "",
        "recommendation_tier": link.recommendation_tier if link else "",
        "predicted_direction": link.predicted_direction if link else "",
        "measurement_ids": tuple(item.measurement_id for item in measurements),
        "measurement_statuses": tuple(item.status.value for item in measurements),
    }


def _ineligible(base: dict[str, object], reason: str) -> CalibrationRecord:
    return CalibrationRecord(
        **base,
        control_experiment_ids=(),
        control_measurement_ids=(),
        eligibility_status="ineligible",
        ineligibility_reasons=(reason,),
    )


def _valid_measurements(measurements: list[MeasurementRecord]) -> list[MeasurementRecord]:
    return [
        item
        for item in measurements
        if item.status is MeasurementStatus.VALID
        and item.canonical_value is not None
        and not item.excluded
    ]


def _control_context_matches(candidate: object, control: object) -> bool:
    return (
        candidate.target_id == control.target_id  # type: ignore[attr-defined]
        and candidate.host == control.host  # type: ignore[attr-defined]
        and candidate.batch_id == control.batch_id  # type: ignore[attr-defined]
        and candidate.condition == control.condition  # type: ignore[attr-defined]
    )


def _experiment_context_complete(experiment: object) -> bool:
    host = experiment.host  # type: ignore[attr-defined]
    condition = experiment.condition  # type: ignore[attr-defined]
    text_values = (
        host.species,
        host.strain,
        host.parent_strain,
        condition.medium,
        condition.carbon_source,
        condition.culture_mode,
        condition.oxygen_or_agitation,
    )
    if any(str(value).strip().lower() in {"", "missing", "unknown"} for value in text_values):
        return False
    return all(
        value is not None
        for value in (condition.temperature_c, condition.ph, condition.sampling_time_h)
    )


def _measurement_signature(measurement: MeasurementRecord) -> tuple[str, str, str, str]:
    return (
        measurement.assay_type,
        measurement.assay_method,
        measurement.canonical_unit,
        measurement.compartment,
    )


def _observed_direction(ratio: float, config: CalibrationConfig) -> str:
    if ratio >= config.increase_threshold_ratio:
        return "increase"
    if ratio <= config.decrease_threshold_ratio:
        return "decrease"
    return "neutral"


def _build_target_summaries(
    records: tuple[CalibrationRecord, ...],
    config: CalibrationConfig,
) -> tuple[TargetCalibrationSummary, ...]:
    summaries: list[TargetCalibrationSummary] = []
    for target_id in sorted({record.target_id for record in records}):
        target_records = [record for record in records if record.target_id == target_id]
        eligible = [record for record in target_records if record.eligibility_status == "eligible"]
        direction_values = [
            record.direction_consistent
            for record in eligible
            if record.direction_consistent is not None
        ]
        summaries.append(
            TargetCalibrationSummary(
                target_id=target_id,
                eligible_count=len(eligible),
                ineligible_count=len(target_records) - len(eligible),
                direction_consistency_rate=(
                    sum(value is True for value in direction_values) / len(direction_values)
                    if direction_values
                    else None
                ),
                rank_correlation=_rank_correlation(eligible, config.minimum_rank_pairs),
                top_k_metrics=tuple(_top_k_metrics(eligible, config)),
                evidence_tier_metrics=tuple(_evidence_tier_metrics(eligible)),
            )
        )
    return tuple(summaries)


def _top_k_metrics(
    records: list[CalibrationRecord],
    config: CalibrationConfig,
) -> list[TopKCalibrationMetric]:
    metrics: list[TopKCalibrationMetric] = []
    for k in config.top_k:
        tested = [
            record
            for record in records
            if record.prediction_rank is not None and record.prediction_rank <= k
        ]
        hit_count = sum(record.hit is True for record in tested)
        hit_rate = hit_count / len(tested) if tested else None
        metrics.append(
            TopKCalibrationMetric(
                k=k,
                tested_count=len(tested),
                hit_count=hit_count,
                hit_rate=hit_rate,
                relative_baseline_enrichment=(
                    hit_rate / config.baseline_hit_rate if hit_rate is not None else None
                ),
            )
        )
    return metrics


def _evidence_tier_metrics(
    records: list[CalibrationRecord],
) -> list[EvidenceTierCalibrationMetric]:
    metrics: list[EvidenceTierCalibrationMetric] = []
    for tier in sorted({record.evidence_tier or "unknown" for record in records}):
        tested = [record for record in records if (record.evidence_tier or "unknown") == tier]
        hit_count = sum(record.hit is True for record in tested)
        metrics.append(
            EvidenceTierCalibrationMetric(
                evidence_tier=tier,
                tested_count=len(tested),
                hit_count=hit_count,
                hit_rate=hit_count / len(tested) if tested else None,
            )
        )
    return metrics


def _rank_correlation(
    records: list[CalibrationRecord],
    minimum_pairs: int,
) -> float | None:
    pairs = [
        (float(record.prediction_rank), float(record.observed_ratio))
        for record in records
        if record.prediction_rank is not None and record.observed_ratio is not None
    ]
    if len(pairs) < minimum_pairs:
        return None
    prediction_ranks = _rank_values([pair[0] for pair in pairs], reverse=False)
    observed_ranks = _rank_values([pair[1] for pair in pairs], reverse=True)
    return _pearson(prediction_ranks, observed_ranks)


def _rank_values(values: list[float], *, reverse: bool) -> list[float]:
    indexed = sorted(enumerate(values), key=lambda item: item[1], reverse=reverse)
    ranks = [0.0] * len(values)
    position = 0
    while position < len(indexed):
        end = position + 1
        while end < len(indexed) and indexed[end][1] == indexed[position][1]:
            end += 1
        average_rank = ((position + 1) + end) / 2.0
        for index, _ in indexed[position:end]:
            ranks[index] = average_rank
        position = end
    return ranks


def _pearson(left: list[float], right: list[float]) -> float | None:
    left_mean = fmean(left)
    right_mean = fmean(right)
    numerator = sum((x - left_mean) * (y - right_mean) for x, y in zip(left, right))
    left_scale = sum((x - left_mean) ** 2 for x in left) ** 0.5
    right_scale = sum((y - right_mean) ** 2 for y in right) ** 0.5
    if left_scale == 0 or right_scale == 0:
        return None
    return numerator / (left_scale * right_scale)


def _json_ready(value: object) -> object:
    if isinstance(value, Enum):
        return value.value
    if isinstance(value, tuple):
        return [_json_ready(item) for item in value]
    if isinstance(value, list):
        return [_json_ready(item) for item in value]
    if isinstance(value, Mapping):
        return {str(key): _json_ready(item) for key, item in value.items()}
    return value


__all__ = [
    "CalibrationConfig",
    "CalibrationOutputs",
    "CalibrationRecord",
    "CalibrationSummary",
    "EvidenceTierCalibrationMetric",
    "TargetCalibrationSummary",
    "TopKCalibrationMetric",
    "build_calibration_summary",
    "write_calibration_outputs",
]
