from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
import math


SCHEMA_VERSION = 1


class ExperimentalFeedbackError(ValueError):
    """Base error for invalid experimental-feedback contracts."""


class SchemaValidationError(ExperimentalFeedbackError):
    """Raised when a canonical record violates the frozen schema contract."""


class UnitValidationError(SchemaValidationError):
    """Raised when a measurement uses an unregistered canonical unit."""


CANONICAL_UNITS = {
    "titer": "mg/L",
    "biomass": "gDCW/L",
    "specific_productivity": "mg/gDCW/h",
    "growth_rate": "1/h",
    "time": "h",
    "viability": "%",
    "od600": "OD600",
}

SUPPORTED_COMPARTMENTS = {
    "extracellular",
    "intracellular",
    "whole_culture",
    "not_applicable",
}


class InterventionType(str, Enum):
    CONTROL = "control"
    KO = "KO"
    OE = "OE"


class MeasurementStatus(str, Enum):
    VALID = "valid"
    BELOW_LOD = "below_lod"
    BELOW_LOQ = "below_loq"
    ABOVE_RANGE = "above_range"
    MISSING = "missing"
    ASSAY_FAILED = "assay_failed"
    EXCLUDED = "excluded"


class PredictionLinkStatus(str, Enum):
    MATCHED = "matched"
    AMBIGUOUS = "ambiguous"
    MISSING_PREDICTION = "missing_prediction"
    CONTEXT_MISMATCH = "context_mismatch"


class QualityStatus(str, Enum):
    VALID = "valid"
    WARNING = "warning"
    INVALID = "invalid"
    EXCLUDED = "excluded"


@dataclass(frozen=True)
class HostContext:
    species: str
    strain: str
    parent_strain: str

    def validate(self) -> None:
        _require_text(self.species, "host.species")
        _require_text(self.strain, "host.strain")
        _require_text(self.parent_strain, "host.parent_strain")


@dataclass(frozen=True)
class ConditionContext:
    medium: str
    carbon_source: str
    culture_mode: str
    temperature_c: float | None
    ph: float | None
    oxygen_or_agitation: str
    sampling_time_h: float | None

    def validate(self) -> None:
        for field_name in ("medium", "carbon_source", "culture_mode", "oxygen_or_agitation"):
            _require_text(getattr(self, field_name), f"condition.{field_name}")
        _validate_optional_number(self.temperature_c, "condition.temperature_c")
        _validate_optional_number(self.ph, "condition.ph")
        _validate_optional_number(self.sampling_time_h, "condition.sampling_time_h")
        if self.sampling_time_h is not None and self.sampling_time_h < 0:
            raise SchemaValidationError("condition.sampling_time_h must be non-negative when provided.")


@dataclass(frozen=True)
class ExperimentRecord:
    experiment_id: str
    target_id: str
    host: HostContext
    batch_id: str
    condition: ConditionContext
    target_name: str = ""
    context_id: str = ""
    biological_replicate_id: str = ""
    operator_id: str = ""
    quality_status: QualityStatus = QualityStatus.VALID
    schema_version: int = SCHEMA_VERSION

    def validate(self) -> None:
        _validate_version(self.schema_version)
        for field_name in ("experiment_id", "target_id", "batch_id"):
            _require_text(getattr(self, field_name), field_name)
        if not isinstance(self.quality_status, QualityStatus):
            raise SchemaValidationError("quality_status must be a QualityStatus value.")
        if self.target_id not in {"hLF", "OPN"}:
            _require_text(self.target_name, "target_name")
        self.host.validate()
        self.condition.validate()


@dataclass(frozen=True)
class ExperimentImportManifest:
    source_file: str
    source_sha256: str
    imported_at: str
    record_count: int
    warnings: tuple[str, ...] = ()
    schema_version: int = SCHEMA_VERSION

    def validate(self) -> None:
        _validate_version(self.schema_version)
        _require_text(self.source_file, "source_file")
        _require_text(self.imported_at, "imported_at")
        if len(self.source_sha256) != 64 or any(
            char not in "0123456789abcdef" for char in self.source_sha256.lower()
        ):
            raise SchemaValidationError("source_sha256 must be a 64-character hex digest.")
        if self.record_count < 0:
            raise SchemaValidationError("record_count must be non-negative.")


@dataclass(frozen=True)
class ExperimentImportConflict:
    code: str
    record_type: str
    record_id: str
    first_payload_json: str
    conflicting_payload_json: str

    def validate(self) -> None:
        for field_name in (
            "code",
            "record_type",
            "record_id",
            "first_payload_json",
            "conflicting_payload_json",
        ):
            _require_text(getattr(self, field_name), field_name)


@dataclass(frozen=True)
class InterventionRecord:
    experiment_id: str
    intervention_id: str
    component_index: int
    intervention_type: InterventionType
    gene_id: str = ""
    common_name: str = ""
    construction_method: str = ""
    construct_id: str = ""
    promoter: str = ""
    induction_mode: str = ""
    copy_number: float | None = None
    prediction_run_id: str = ""
    evidence_id: str = ""
    warnings: tuple[str, ...] = ()
    schema_version: int = SCHEMA_VERSION

    def validate(self) -> None:
        _validate_version(self.schema_version)
        _require_text(self.experiment_id, "experiment_id")
        _require_text(self.intervention_id, "intervention_id")
        if not isinstance(self.intervention_type, InterventionType):
            raise SchemaValidationError("intervention_type must be an InterventionType value.")
        if not isinstance(self.component_index, int) or isinstance(self.component_index, bool) or self.component_index < 1:
            raise SchemaValidationError("component_index must be at least 1.")
        if self.intervention_type is InterventionType.KO:
            _require_text(self.gene_id, "gene_id")
            _require_text(self.construction_method, "construction_method")
        elif self.intervention_type is InterventionType.OE:
            _require_text(self.gene_id, "gene_id")
            _require_text(self.construct_id, "construct_id")
            _require_text(self.promoter, "promoter")
            _require_text(self.induction_mode, "induction_mode")
            if self.copy_number is None and "copy_number_unknown" not in self.warnings:
                raise SchemaValidationError(
                    "OE with unknown copy_number requires copy_number_unknown warning."
                )
            if self.copy_number is not None:
                _validate_optional_number(self.copy_number, "copy_number")
                if self.copy_number <= 0:
                    raise SchemaValidationError("copy_number must be positive when provided.")


@dataclass(frozen=True)
class MeasurementRecord:
    experiment_id: str
    measurement_id: str
    assay_type: str
    assay_method: str
    compartment: str
    raw_value: float | None
    raw_unit: str
    canonical_value: float | None
    canonical_unit: str
    status: MeasurementStatus
    technical_replicate_id: str = ""
    status_reason: str = ""
    excluded: bool = False
    exclusion_reason: str = ""
    reviewer_id: str = ""
    schema_version: int = SCHEMA_VERSION

    def validate(self) -> None:
        _validate_version(self.schema_version)
        for field_name in (
            "experiment_id",
            "measurement_id",
            "assay_type",
            "assay_method",
            "compartment",
            "raw_unit",
            "canonical_unit",
        ):
            _require_text(getattr(self, field_name), field_name)
        if not isinstance(self.status, MeasurementStatus):
            raise SchemaValidationError("status must be a MeasurementStatus value.")
        if self.compartment not in SUPPORTED_COMPARTMENTS:
            raise SchemaValidationError(
                f"compartment must be one of {sorted(SUPPORTED_COMPARTMENTS)}."
            )
        _validate_optional_number(self.raw_value, "raw_value")
        _validate_optional_number(self.canonical_value, "canonical_value")
        expected_unit = CANONICAL_UNITS.get(self.assay_type)
        if expected_unit is None:
            raise UnitValidationError(f"unsupported assay_type: {self.assay_type}")
        if self.canonical_unit != expected_unit:
            raise UnitValidationError(
                f"{self.assay_type} canonical_unit must be {expected_unit}; got {self.canonical_unit}."
            )
        if self.raw_unit != expected_unit:
            raise UnitValidationError(
                f"raw_unit {self.raw_unit} requires an explicit conversion registry entry to {expected_unit}."
            )
        if self.status is MeasurementStatus.VALID and (
            self.raw_value is None or self.canonical_value is None
        ):
            raise SchemaValidationError("valid measurements require raw_value and canonical_value.")
        if self.status is not MeasurementStatus.VALID and (
            self.raw_value == 0 or self.canonical_value == 0
        ):
            raise SchemaValidationError("non-valid measurement states must not be encoded as zero.")
        if self.status in {
            MeasurementStatus.BELOW_LOD,
            MeasurementStatus.BELOW_LOQ,
            MeasurementStatus.ABOVE_RANGE,
            MeasurementStatus.MISSING,
            MeasurementStatus.ASSAY_FAILED,
            MeasurementStatus.EXCLUDED,
        }:
            _require_text(self.status_reason, "status_reason")
        if self.excluded and not self.exclusion_reason:
            raise SchemaValidationError("excluded measurements require exclusion_reason.")
        if self.excluded != (self.status is MeasurementStatus.EXCLUDED):
            raise SchemaValidationError("excluded flag and measurement status must agree.")


@dataclass(frozen=True)
class PredictionLinkRecord:
    experiment_id: str
    intervention_id: str
    prediction_run_id: str
    evidence_id: str
    target_id: str
    gene_id: str
    intervention_type: InterventionType
    status: PredictionLinkStatus
    common_name: str = ""
    reaction_id: str = ""
    reason: str = ""
    prediction_rank: int | None = None
    predicted_direction: str = ""
    evidence_tier: str = ""
    recommendation_tier: str = ""
    prediction_score: float | None = None
    prediction_context_id: str = ""
    schema_version: int = SCHEMA_VERSION

    def validate(self) -> None:
        _validate_version(self.schema_version)
        for field_name in ("experiment_id", "intervention_id", "target_id"):
            _require_text(getattr(self, field_name), field_name)
        if not isinstance(self.intervention_type, InterventionType):
            raise SchemaValidationError("intervention_type must be an InterventionType value.")
        if not isinstance(self.status, PredictionLinkStatus):
            raise SchemaValidationError("status must be a PredictionLinkStatus value.")
        if self.status is PredictionLinkStatus.MATCHED:
            for field_name in ("prediction_run_id", "evidence_id", "gene_id"):
                _require_text(getattr(self, field_name), field_name)
        else:
            _require_text(self.reason, "reason")
        if self.prediction_rank is not None and (
            not isinstance(self.prediction_rank, int)
            or isinstance(self.prediction_rank, bool)
            or self.prediction_rank < 1
        ):
            raise SchemaValidationError("prediction_rank must be a positive integer or None.")
        _validate_optional_number(self.prediction_score, "prediction_score")
        if self.predicted_direction not in {"", "increase", "decrease", "neutral"}:
            raise SchemaValidationError("predicted_direction must be increase, decrease, neutral, or empty.")


@dataclass(frozen=True)
class ExperimentBundle:
    experiments: tuple[ExperimentRecord, ...]
    interventions: tuple[InterventionRecord, ...]
    measurements: tuple[MeasurementRecord, ...] = ()
    prediction_links: tuple[PredictionLinkRecord, ...] = ()
    source_file: str = field(default="", compare=False)
    warnings: tuple[str, ...] = field(default=(), compare=False)
    import_manifest: ExperimentImportManifest | None = field(default=None, compare=False)
    import_conflicts: tuple[ExperimentImportConflict, ...] = field(default=(), compare=False)
    schema_version: int = SCHEMA_VERSION

    def validate(self) -> None:
        _validate_version(self.schema_version)
        if not self.experiments:
            raise SchemaValidationError("bundle must contain at least one experiment.")
        experiment_ids: set[str] = set()
        for experiment in self.experiments:
            experiment.validate()
            if experiment.experiment_id in experiment_ids:
                raise SchemaValidationError(f"duplicate experiment_id: {experiment.experiment_id}")
            experiment_ids.add(experiment.experiment_id)
        intervention_keys: set[tuple[str, str]] = set()
        component_keys: set[tuple[str, int]] = set()
        for intervention in self.interventions:
            intervention.validate()
            if intervention.experiment_id not in experiment_ids:
                raise SchemaValidationError(
                    f"intervention references missing experiment_id: {intervention.experiment_id}"
                )
            key = (intervention.experiment_id, intervention.intervention_id)
            if key in intervention_keys:
                raise SchemaValidationError(f"duplicate intervention_id within experiment: {key}")
            intervention_keys.add(key)
            component_key = (intervention.experiment_id, intervention.component_index)
            if component_key in component_keys:
                raise SchemaValidationError(
                    f"duplicate component_index within experiment: {component_key}"
                )
            component_keys.add(component_key)
        measurement_keys: set[tuple[str, str]] = set()
        for measurement in self.measurements:
            measurement.validate()
            if measurement.experiment_id not in experiment_ids:
                raise SchemaValidationError(
                    f"measurement references missing experiment_id: {measurement.experiment_id}"
                )
            key = (measurement.experiment_id, measurement.measurement_id)
            if key in measurement_keys:
                raise SchemaValidationError(f"duplicate measurement_id within experiment: {key}")
            measurement_keys.add(key)
        for link in self.prediction_links:
            link.validate()
            key = (link.experiment_id, link.intervention_id)
            if key not in intervention_keys:
                raise SchemaValidationError(f"prediction link references missing intervention: {key}")


def _validate_version(value: int) -> None:
    if value != SCHEMA_VERSION:
        raise SchemaValidationError(f"schema_version must be {SCHEMA_VERSION}.")


def _require_text(value: object, field_name: str) -> None:
    if not str(value or "").strip():
        raise SchemaValidationError(f"{field_name} must be non-empty.")


def _validate_optional_number(value: object, field_name: str) -> None:
    if value is not None and (not isinstance(value, (int, float)) or isinstance(value, bool)):
        raise SchemaValidationError(f"{field_name} must be numeric or None.")
    if value is not None and not math.isfinite(float(value)):
        raise SchemaValidationError(f"{field_name} must be finite when provided.")


__all__ = [
    "SCHEMA_VERSION",
    "CANONICAL_UNITS",
    "SUPPORTED_COMPARTMENTS",
    "ConditionContext",
    "ExperimentBundle",
    "ExperimentImportManifest",
    "ExperimentImportConflict",
    "ExperimentRecord",
    "ExperimentalFeedbackError",
    "HostContext",
    "InterventionRecord",
    "InterventionType",
    "MeasurementRecord",
    "MeasurementStatus",
    "PredictionLinkRecord",
    "PredictionLinkStatus",
    "QualityStatus",
    "SchemaValidationError",
    "UnitValidationError",
]
