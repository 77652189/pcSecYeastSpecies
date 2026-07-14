from __future__ import annotations

import math
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Mapping

from pcsec_pichia.external_refs.capacity_sources import (
    ExternalCapacitySource,
    ExternalCapacitySourceType,
    RetrievalMode,
)
from pcsec_pichia.external_refs.schema import utc_now_iso
from pcsec_pichia.oe_capacity.schema import OECapacityValidationError




CANDIDATE_SCHEMA_VERSION = 1
CANDIDATE_RECORDS_FILENAME = "external_capacity_candidates.jsonl"
CANDIDATE_MANIFEST_FILENAME = "external_capacity_candidate_manifest.json"
RAW_MEASUREMENTS_FILENAME = "raw_capacity_measurements.jsonl"
PROMOTION_MANIFEST_FILENAME = "capacity_promotion_manifest.json"


class CapacityApplicabilityScope(str, Enum):
    TARGET_SPECIFIC = "target_specific"
    HOST_CONDITION = "host_condition"
    EXTERNAL_MODEL_CALIBRATED = "external_model_calibrated"
    HOMOLOG_TRANSFERRED = "homolog_transferred"


class CapacityParameterKind(str, Enum):
    ABUNDANCE = "abundance"
    KCAT = "kcat"
    BASELINE_CAPACITY = "baseline_capacity"


class CapacityConfidence(str, Enum):
    HIGH = "high"
    MEDIUM = "medium"
    LOW = "low"


class CapacityCandidateStatus(str, Enum):
    RAW = "raw"
    REVIEW_REQUIRED = "review_required"
    REVIEW_READY = "review_ready"
    REJECTED = "rejected"
    PROMOTED = "promoted"


class PromotionDecision(str, Enum):
    PENDING = "pending"
    APPROVED = "approved"
    REJECTED = "rejected"


@dataclass(frozen=True)
class HostCondition:
    species: str
    strain: str
    medium: str
    carbon_source: str
    culture_mode: str
    growth_rate_per_h: float
    temperature_c: float | None = None
    ph: float | None = None
    oxygen_condition: str = ""
    biomass_basis: str = "gDW"

    def validate(self) -> None:
        for name in ("species", "strain", "medium", "carbon_source", "culture_mode", "biomass_basis"):
            _require_text(getattr(self, name), name)
        _require_positive(self.growth_rate_per_h, "growth_rate_per_h")
        if self.temperature_c is not None:
            _require_positive(self.temperature_c, "temperature_c")
        if self.ph is not None:
            _require_positive(self.ph, "ph")


@dataclass(frozen=True)
class RawCapacityMeasurement:
    measurement_id: str
    source_id: str
    parameter_kind: CapacityParameterKind
    nominal_value: float
    lower_bound: float
    upper_bound: float
    unit: str
    condition: HostCondition
    external_gene_id: str = ""
    external_protein_id: str = ""
    external_enzyme_id: str = ""
    biomass_basis: str = ""
    notes: str = ""

    def validate(self) -> None:
        for name in ("measurement_id", "source_id", "unit"):
            _require_text(getattr(self, name), name)
        if not any((self.external_gene_id, self.external_protein_id, self.external_enzyme_id)):
            raise OECapacityValidationError("measurement requires an external gene/protein/enzyme identifier.")
        for name in ("nominal_value", "lower_bound", "upper_bound"):
            _require_positive(getattr(self, name), name)
        if not self.lower_bound <= self.nominal_value <= self.upper_bound:
            raise OECapacityValidationError("measurement interval must contain nominal_value.")
        self.condition.validate()


@dataclass(frozen=True)
class CapacityModelBinding:
    target_id: str
    context_id: str
    mapping_id: str
    model_fingerprint: str
    gene_id: str
    enzyme_id: str
    reaction_id: str
    formation_or_dilution_reaction_id: str
    mapping_evidence: tuple[str, ...]
    external_gene_id: str = ""
    external_protein_id: str = ""
    external_enzyme_id: str = ""

    def validate(self) -> None:
        for name in (
            "target_id",
            "context_id",
            "mapping_id",
            "model_fingerprint",
            "gene_id",
            "enzyme_id",
            "reaction_id",
            "formation_or_dilution_reaction_id",
        ):
            _require_text(getattr(self, name), name)
        if not self.mapping_evidence:
            raise OECapacityValidationError("crosswalk requires mapping_evidence.")


@dataclass(frozen=True)
class CapacityConversionStep:
    step_id: str
    input_value: float
    input_unit: str
    output_value: float
    output_unit: str
    formula: str
    factor: float
    source_ref: str
    missing_metadata: tuple[str, ...] = ()

    def validate(self) -> None:
        for name in ("step_id", "input_unit", "output_unit", "formula", "source_ref"):
            _require_text(getattr(self, name), name)
        for name in ("input_value", "output_value", "factor"):
            _require_positive(getattr(self, name), name)


@dataclass(frozen=True)
class ExternalCapacityCandidate:
    candidate_id: str
    applicability_scope: CapacityApplicabilityScope
    source_ids: tuple[str, ...]
    measurement_ids: tuple[str, ...]
    model_bindings: tuple[CapacityModelBinding, ...]
    condition: HostCondition
    nominal_capacity: float | None
    lower_capacity: float | None
    upper_capacity: float | None
    unit: str
    confidence: CapacityConfidence
    status: CapacityCandidateStatus
    conversion_steps: tuple[CapacityConversionStep, ...]
    target_id: str = ""
    conflicts: tuple[str, ...] = ()
    missing_information: tuple[str, ...] = ()
    rejection_reasons: tuple[str, ...] = ()
    warnings: tuple[str, ...] = ()

    def validate(self) -> None:
        _require_text(self.candidate_id, "candidate_id")
        if not self.source_ids or not self.measurement_ids:
            raise OECapacityValidationError("candidate requires source_ids and measurement_ids.")
        if not self.model_bindings:
            raise OECapacityValidationError("candidate requires at least one current-model binding.")
        for binding in self.model_bindings:
            binding.validate()
        self.condition.validate()
        expected_context_id = _condition_context_id(self.condition)
        if any(
            binding.context_id != expected_context_id
            for binding in self.model_bindings
        ):
            raise OECapacityValidationError(
                "candidate model binding context does not match the structured host condition."
            )
        values = (self.lower_capacity, self.nominal_capacity, self.upper_capacity)
        if any(value is not None for value in values):
            if any(value is None for value in values):
                raise OECapacityValidationError("capacity interval must be fully populated or fully unavailable.")
            for name in ("nominal_capacity", "lower_capacity", "upper_capacity"):
                _require_positive(getattr(self, name), name)
            if not float(self.lower_capacity) <= float(self.nominal_capacity) <= float(self.upper_capacity):
                raise OECapacityValidationError("capacity interval must contain nominal_capacity.")
        if self.unit != "model_flux":
            raise OECapacityValidationError("candidate canonical unit must be model_flux.")
        if self.applicability_scope is CapacityApplicabilityScope.TARGET_SPECIFIC and not self.target_id:
            raise OECapacityValidationError("target_specific candidate requires target_id.")
        if self.applicability_scope is not CapacityApplicabilityScope.TARGET_SPECIFIC and self.target_id:
            raise OECapacityValidationError("only target_specific candidates may set target_id.")
        if self.target_id and any(
            binding.target_id != self.target_id for binding in self.model_bindings
        ):
            raise OECapacityValidationError(
                "target_specific candidate bindings must match target_id."
            )
        if self.applicability_scope is CapacityApplicabilityScope.HOMOLOG_TRANSFERRED:
            if self.confidence is not CapacityConfidence.LOW:
                raise OECapacityValidationError("homolog_transferred candidate confidence must be low.")
            if self.status in {CapacityCandidateStatus.REVIEW_READY, CapacityCandidateStatus.PROMOTED}:
                raise OECapacityValidationError("homolog transfer alone cannot become review-ready or promoted.")
        if self.status in {CapacityCandidateStatus.REVIEW_READY, CapacityCandidateStatus.PROMOTED}:
            if self.nominal_capacity is None:
                raise OECapacityValidationError("review-ready candidate requires a canonical capacity interval.")
            if self.conflicts or self.missing_information or self.rejection_reasons:
                raise OECapacityValidationError("review-ready candidate cannot retain conflicts or missing/rejected evidence.")
            if not self.conversion_steps:
                raise OECapacityValidationError("review-ready candidate requires a traceable conversion chain.")
        for step in self.conversion_steps:
            step.validate()


@dataclass(frozen=True)
class ExternalCapacityCandidateBundle:
    model_fingerprints: tuple[str, ...]
    sources: tuple[ExternalCapacitySource, ...]
    measurements: tuple[RawCapacityMeasurement, ...]
    candidates: tuple[ExternalCapacityCandidate, ...]
    schema_version: int = CANDIDATE_SCHEMA_VERSION
    generated_at: str = field(default_factory=utc_now_iso)

    def validate(self) -> None:
        if self.schema_version != CANDIDATE_SCHEMA_VERSION:
            raise OECapacityValidationError(f"candidate bundle requires schema_version={CANDIDATE_SCHEMA_VERSION}.")
        if not self.model_fingerprints:
            raise OECapacityValidationError("candidate bundle requires model_fingerprints.")
        for fingerprint in self.model_fingerprints:
            _require_text(fingerprint, "model_fingerprint")
        source_ids: set[str] = set()
        for source in self.sources:
            source.validate()
            if source.source_id in source_ids:
                raise OECapacityValidationError("duplicate external capacity source_id.")
            source_ids.add(source.source_id)
        measurement_ids: set[str] = set()
        for measurement in self.measurements:
            measurement.validate()
            if measurement.measurement_id in measurement_ids:
                raise OECapacityValidationError("duplicate capacity measurement_id.")
            if measurement.source_id not in source_ids:
                raise OECapacityValidationError("measurement references an unknown source_id.")
            measurement_ids.add(measurement.measurement_id)
        candidate_ids: set[str] = set()
        for candidate in self.candidates:
            candidate.validate()
            if candidate.candidate_id in candidate_ids:
                raise OECapacityValidationError("duplicate external capacity candidate_id.")
            if not {binding.model_fingerprint for binding in candidate.model_bindings}.issubset(
                set(self.model_fingerprints)
            ):
                raise OECapacityValidationError("candidate binding model_fingerprint mismatch.")
            if not set(candidate.source_ids).issubset(source_ids):
                raise OECapacityValidationError("candidate references an unknown source_id.")
            if not set(candidate.measurement_ids).issubset(measurement_ids):
                raise OECapacityValidationError("candidate references an unknown measurement_id.")
            candidate_ids.add(candidate.candidate_id)


@dataclass(frozen=True)
class CapacityPromotionManifest:
    decision: PromotionDecision
    candidate_ids: tuple[str, ...]
    model_fingerprints: tuple[str, ...]
    candidate_bundle_sha256: str
    asset_path: str
    reviewer: str = ""
    reviewed_at: str = ""
    reason: str = ""
    promoted_asset_sha256: str = ""
    expected_asset_sha256: str = ""
    schema_version: int = 1

    def validate(self) -> None:
        if self.schema_version != 1:
            raise OECapacityValidationError("promotion manifest requires schema_version=1.")
        if not self.candidate_ids:
            raise OECapacityValidationError("promotion manifest requires candidate_ids.")
        if not self.model_fingerprints:
            raise OECapacityValidationError("promotion manifest requires model_fingerprints.")
        _require_sha256(self.candidate_bundle_sha256, "candidate_bundle_sha256")
        _require_text(self.asset_path, "asset_path")
        if self.expected_asset_sha256 != "missing":
            _require_sha256(self.expected_asset_sha256, "expected_asset_sha256")
        if self.decision is PromotionDecision.APPROVED:
            _require_text(self.reviewer, "reviewer")
            _require_text(self.reviewed_at, "reviewed_at")


def validate_external_capacity_candidate_bundle(
    bundle: ExternalCapacityCandidateBundle,
) -> None:
    bundle.validate()


def _condition_context_id(condition: HostCondition) -> str:
    return f"{condition.carbon_source.strip().lower()}_mu_{condition.growth_rate_per_h:g}"


def _as_object(value: object, label: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise OECapacityValidationError(f"{label} must be an object.")
    return value


def _as_object_list(value: object, label: str) -> list[Mapping[str, Any]]:
    if not isinstance(value, list):
        raise OECapacityValidationError(f"{label} must be an array.")
    return [_as_object(item, label) for item in value]


def _require_text(value: object, field_name: str) -> None:
    if not isinstance(value, str) or not value.strip():
        raise OECapacityValidationError(f"{field_name} must be non-empty text.")


def _require_positive(value: object, field_name: str) -> None:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(float(value)) or float(value) <= 0:
        raise OECapacityValidationError(f"{field_name} must be a positive finite number.")


def _require_sha256(value: str, field_name: str) -> None:
    if len(value) != 64 or any(character not in "0123456789abcdef" for character in value.lower()):
        raise OECapacityValidationError(f"{field_name} must be a 64-character hex digest.")


def _float(value: object, field_name: str) -> float:
    if isinstance(value, bool):
        raise OECapacityValidationError(f"{field_name} must be numeric.")
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise OECapacityValidationError(f"{field_name} must be numeric.") from exc
    if not math.isfinite(parsed):
        raise OECapacityValidationError(f"{field_name} must be finite.")
    return parsed


def _optional_float(value: object, field_name: str) -> float | None:
    if value in (None, ""):
        return None
    return _float(value, field_name)
