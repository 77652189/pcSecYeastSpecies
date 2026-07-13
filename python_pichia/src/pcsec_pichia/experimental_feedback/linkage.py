from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Iterable, Mapping

from pcsec_pichia.experimental_feedback.schema import (
    ExperimentBundle,
    InterventionType,
    PredictionLinkRecord,
    PredictionLinkStatus,
    SchemaValidationError,
)


@dataclass(frozen=True)
class PredictionRecord:
    prediction_run_id: str
    evidence_id: str
    target_id: str
    gene_id: str
    intervention_type: InterventionType
    context_id: str = ""
    rank: int | None = None
    recommendation_tier: str = ""
    evidence_tier: str = ""
    reaction_id: str = ""
    common_name: str = ""
    predicted_direction: str = ""
    score: float | None = None

    def validate(self) -> None:
        for field_name in (
            "prediction_run_id",
            "evidence_id",
            "target_id",
            "gene_id",
        ):
            if not str(getattr(self, field_name) or "").strip():
                raise SchemaValidationError(f"prediction {field_name} must be non-empty.")
        if not isinstance(self.intervention_type, InterventionType):
            raise SchemaValidationError("prediction intervention_type must be InterventionType.")
        if self.rank is not None and (
            not isinstance(self.rank, int) or isinstance(self.rank, bool) or self.rank < 1
        ):
            raise SchemaValidationError("prediction rank must be a positive integer or None.")
        if self.score is not None and (
            not isinstance(self.score, (int, float))
            or isinstance(self.score, bool)
            or not math.isfinite(float(self.score))
        ):
            raise SchemaValidationError("prediction score must be a finite number or None.")
        if self.predicted_direction not in {"", "increase", "decrease", "neutral"}:
            raise SchemaValidationError("prediction direction must be increase, decrease, neutral, or empty.")


@dataclass(frozen=True)
class PredictionIndex:
    records: tuple[PredictionRecord, ...]

    def validate(self) -> None:
        seen: set[tuple[str, str]] = set()
        for record in self.records:
            record.validate()
            key = (record.prediction_run_id, record.evidence_id)
            if key in seen:
                raise SchemaValidationError(f"duplicate prediction identity: {key}")
            seen.add(key)


@dataclass(frozen=True)
class PredictionLinkageResult:
    links: tuple[PredictionLinkRecord, ...]
    control_count: int = 0

    @property
    def matched_count(self) -> int:
        return self._count(PredictionLinkStatus.MATCHED)

    @property
    def ambiguous_count(self) -> int:
        return self._count(PredictionLinkStatus.AMBIGUOUS)

    @property
    def missing_prediction_count(self) -> int:
        return self._count(PredictionLinkStatus.MISSING_PREDICTION)

    @property
    def context_mismatch_count(self) -> int:
        return self._count(PredictionLinkStatus.CONTEXT_MISMATCH)

    def _count(self, status: PredictionLinkStatus) -> int:
        return sum(1 for link in self.links if link.status is status)


def build_prediction_index(screen_runs: Iterable[Mapping[str, Any]]) -> PredictionIndex:
    records: list[PredictionRecord] = []
    for screen_run in screen_runs:
        if not isinstance(screen_run, Mapping):
            raise SchemaValidationError("screen run must be a mapping.")
        run_id = str(
            screen_run.get("prediction_run_id")
            or screen_run.get("source_run")
            or screen_run.get("run_id")
            or ""
        )
        raw_items = screen_run.get("evidence_items") or screen_run.get("rows") or ()
        if not isinstance(raw_items, Iterable) or isinstance(raw_items, (str, bytes, Mapping)):
            raise SchemaValidationError("screen run evidence_items must be an iterable of mappings.")
        for default_rank, item in enumerate(raw_items, start=1):
            if not isinstance(item, Mapping):
                raise SchemaValidationError("prediction row must be a mapping.")
            record = _prediction_from_mapping(
                item,
                default_run_id=run_id,
                default_rank=default_rank,
            )
            record.validate()
            records.append(record)
    index = PredictionIndex(records=tuple(records))
    index.validate()
    return index


def link_experiments_to_predictions(
    bundle: ExperimentBundle,
    prediction_index: PredictionIndex,
) -> PredictionLinkageResult:
    prediction_index.validate()
    experiments = {}
    for record in bundle.experiments:
        if record.experiment_id in experiments:
            raise SchemaValidationError(f"duplicate experiment_id during linkage: {record.experiment_id}")
        experiments[record.experiment_id] = record
    links: list[PredictionLinkRecord] = []
    control_count = 0
    intervention_keys: set[tuple[str, str]] = set()
    for intervention in bundle.interventions:
        key = (intervention.experiment_id, intervention.intervention_id)
        if key in intervention_keys:
            raise SchemaValidationError(f"duplicate intervention identity during linkage: {key}")
        intervention_keys.add(key)
        if not isinstance(intervention.intervention_type, InterventionType):
            raise SchemaValidationError("intervention_type must be canonical before linkage.")
        if intervention.intervention_type is InterventionType.CONTROL:
            control_count += 1
            continue
        experiment = experiments.get(intervention.experiment_id)
        if experiment is None:
            links.append(_unmatched(intervention, "", PredictionLinkStatus.MISSING_PREDICTION, "missing_experiment"))
            continue
        if not intervention.gene_id:
            reason = "common_name_only" if intervention.common_name else "gene_id_missing"
            links.append(_unmatched(intervention, experiment.target_id, PredictionLinkStatus.AMBIGUOUS, reason))
            continue
        if not intervention.prediction_run_id:
            links.append(
                _unmatched(
                    intervention,
                    experiment.target_id,
                    PredictionLinkStatus.MISSING_PREDICTION,
                    "prediction_run_id_missing",
                )
            )
            continue
        run_gene_type = [
            record
            for record in prediction_index.records
            if record.prediction_run_id == intervention.prediction_run_id
            and record.gene_id == intervention.gene_id
            and record.intervention_type is intervention.intervention_type
        ]
        if not run_gene_type:
            links.append(
                _unmatched(
                    intervention,
                    experiment.target_id,
                    PredictionLinkStatus.MISSING_PREDICTION,
                    "prediction_not_found",
                )
            )
            continue
        target_matches = [record for record in run_gene_type if record.target_id == experiment.target_id]
        if not target_matches:
            links.append(
                _unmatched(
                    intervention,
                    experiment.target_id,
                    PredictionLinkStatus.CONTEXT_MISMATCH,
                    "target_mismatch",
                )
            )
            continue
        candidates = target_matches
        if intervention.evidence_id:
            candidates = [record for record in candidates if record.evidence_id == intervention.evidence_id]
            if not candidates:
                links.append(
                    _unmatched(
                        intervention,
                        experiment.target_id,
                        PredictionLinkStatus.MISSING_PREDICTION,
                        "evidence_id_not_found",
                    )
                )
                continue
        if len(candidates) != 1:
            links.append(
                _unmatched(
                    intervention,
                    experiment.target_id,
                    PredictionLinkStatus.AMBIGUOUS,
                    "multiple_predictions_match",
                )
            )
            continue
        prediction = candidates[0]
        if (experiment.context_id or prediction.context_id) and experiment.context_id != prediction.context_id:
            reason = "context_id_mismatch"
            if not experiment.context_id:
                reason = "experiment_context_missing"
            elif not prediction.context_id:
                reason = "prediction_context_missing"
            links.append(
                _unmatched(
                    intervention,
                    experiment.target_id,
                    PredictionLinkStatus.CONTEXT_MISMATCH,
                    reason,
                    prediction=prediction,
                )
            )
            continue
        links.append(
            PredictionLinkRecord(
                experiment_id=experiment.experiment_id,
                intervention_id=intervention.intervention_id,
                prediction_run_id=prediction.prediction_run_id,
                evidence_id=prediction.evidence_id,
                target_id=prediction.target_id,
                gene_id=prediction.gene_id,
                intervention_type=prediction.intervention_type,
                status=PredictionLinkStatus.MATCHED,
                common_name=intervention.common_name,
                reaction_id=prediction.reaction_id,
            )
        )
    result = PredictionLinkageResult(links=tuple(links), control_count=control_count)
    for link in result.links:
        link.validate()
    return result


def _prediction_from_mapping(
    item: Mapping[str, Any],
    *,
    default_run_id: str,
    default_rank: int,
) -> PredictionRecord:
    intervention_text = str(item.get("intervention_type") or "").strip()
    try:
        intervention_type = InterventionType(intervention_text)
    except ValueError as exc:
        raise SchemaValidationError(f"unsupported prediction intervention_type: {intervention_text}") from exc
    score = item.get("score")
    if score is None:
        score = item.get("secretion_ratio_vs_wildtype")
    return PredictionRecord(
        prediction_run_id=str(item.get("prediction_run_id") or item.get("source_run") or default_run_id),
        evidence_id=str(item.get("evidence_id") or ""),
        target_id=str(item.get("target_id") or item.get("target_key") or ""),
        gene_id=str(item.get("canonical_gene_id") or item.get("gene_id") or ""),
        intervention_type=intervention_type,
        context_id=str(item.get("context_id") or item.get("medium_condition_id") or ""),
        rank=_optional_int(item.get("rank")) or default_rank,
        recommendation_tier=str(item.get("recommendation_tier") or ""),
        evidence_tier=str(item.get("evidence_tier") or item.get("recommendation_tier") or ""),
        reaction_id=str(item.get("reaction_id") or ""),
        common_name=str(item.get("gene_display_name") or item.get("common_name") or ""),
        predicted_direction=_prediction_direction(item),
        score=_optional_float(score),
    )


def _unmatched(
    intervention: Any,
    target_id: str,
    status: PredictionLinkStatus,
    reason: str,
    *,
    prediction: PredictionRecord | None = None,
) -> PredictionLinkRecord:
    return PredictionLinkRecord(
        experiment_id=intervention.experiment_id,
        intervention_id=intervention.intervention_id,
        prediction_run_id=(prediction.prediction_run_id if prediction else intervention.prediction_run_id),
        evidence_id=(prediction.evidence_id if prediction else intervention.evidence_id),
        target_id=target_id,
        gene_id=intervention.gene_id,
        intervention_type=intervention.intervention_type,
        status=status,
        common_name=intervention.common_name,
        reaction_id=prediction.reaction_id if prediction else "",
        reason=reason,
    )


def _optional_int(value: object) -> int | None:
    if value is None or value == "":
        return None
    if isinstance(value, bool):
        raise SchemaValidationError(f"invalid prediction rank: {value}")
    if isinstance(value, float) and not value.is_integer():
        raise SchemaValidationError(f"invalid prediction rank: {value}")
    try:
        return int(value)
    except (TypeError, ValueError) as exc:
        raise SchemaValidationError(f"invalid prediction rank: {value}") from exc


def _optional_float(value: object) -> float | None:
    if value is None or value == "":
        return None
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise SchemaValidationError(f"invalid prediction score: {value}") from exc
    if not math.isfinite(parsed):
        raise SchemaValidationError("prediction score must be finite.")
    return parsed


def _prediction_direction(item: Mapping[str, Any]) -> str:
    explicit = str(item.get("predicted_direction") or item.get("effect_direction") or "").strip()
    if explicit:
        return explicit
    parsed = _optional_float(item.get("secretion_ratio_vs_wildtype"))
    if parsed is None:
        return ""
    if parsed > 1.0:
        return "increase"
    if parsed < 1.0:
        return "decrease"
    return "neutral"


__all__ = [
    "PredictionIndex",
    "PredictionLinkageResult",
    "PredictionRecord",
    "build_prediction_index",
    "link_experiments_to_predictions",
]
