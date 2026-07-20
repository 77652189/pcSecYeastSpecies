from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum


SCHEMA_VERSION = 1


class SecretoryResourceError(ValueError):
    """Base error for invalid secretory-resource contracts."""


class SecretoryResourceValidationError(SecretoryResourceError):
    """Raised when a catalog entry violates the frozen Round 0 contract."""


class ResourceCategory(str, Enum):
    ER_TRANSLOCATION = "er_translocation"
    FOLDING_CHAPERONE = "folding_chaperone"
    DISULFIDE_BOND_FORMATION = "disulfide_bond_formation"
    GLYCOSYLATION = "glycosylation"
    VESICLE_TRAFFICKING = "vesicle_trafficking"
    ER_QUALITY_CONTROL = "er_quality_control"
    TARGET_SPECIFIC_COST = "target_specific_cost"


IN_SCOPE_CATEGORIES: tuple[ResourceCategory, ...] = (
    ResourceCategory.ER_TRANSLOCATION,
    ResourceCategory.FOLDING_CHAPERONE,
    ResourceCategory.DISULFIDE_BOND_FORMATION,
    ResourceCategory.GLYCOSYLATION,
    ResourceCategory.VESICLE_TRAFFICKING,
    ResourceCategory.ER_QUALITY_CONTROL,
    ResourceCategory.TARGET_SPECIFIC_COST,
)

# All seven pichia_next_plan.md Round 0 categories are now in scope; none
# remain deferred. Kept as an explicit empty tuple (rather than removed)
# because validation.py's out-of-scope gate and external callers key off
# this name existing.
DEFERRED_CATEGORIES: tuple[ResourceCategory, ...] = ()


class ExecutionStatus(str, Enum):
    EXECUTABLE = "executable"
    EVIDENCE_ONLY = "evidence_only"
    UNAVAILABLE = "unavailable"
    NOT_APPLICABLE = "not_applicable"
    CONFLICT = "conflict"


class CalibrationMode(str, Enum):
    """Mirrors the direction-2 relative/absolute split (ADR-002).

    Round 0 never produces ABSOLUTE_CALIBRATED; the member exists so a later
    round can add an absolute layer without redefining this enum.
    """

    RELATIVE_UNCALIBRATED = "relative_uncalibrated"
    ABSOLUTE_CALIBRATED = "absolute_calibrated"


class EvidenceClass(str, Enum):
    CURRENT_MODEL_HANDLE = "current_model_handle"
    TARGET_STRUCTURAL_PROFILE = "target_structural_profile"
    CLASSIFIER_INFERRED = "classifier_inferred"


@dataclass(frozen=True)
class ResourceSource:
    source_ref: str
    version: str
    evidence_class: EvidenceClass
    license: str = "internal_repo_asset"
    reviewed: bool = False

    def validate(self) -> None:
        _require_text(self.source_ref, "source.source_ref")
        _require_text(self.version, "source.version")
        if not isinstance(self.evidence_class, EvidenceClass):
            raise SecretoryResourceValidationError("source.evidence_class must be an EvidenceClass value.")
        _require_text(self.license, "source.license")


@dataclass(frozen=True)
class ResourceApplicability:
    host: str
    target_id: str
    condition: str
    model_fingerprint: str | None = None

    def validate(self) -> None:
        _require_text(self.host, "applicability.host")
        _require_text(self.target_id, "applicability.target_id")
        _require_text(self.condition, "applicability.condition")


@dataclass(frozen=True)
class SecretoryResource:
    resource_id: str
    category: ResourceCategory
    canonical_unit: str
    model_handles: tuple[str, ...]
    source: ResourceSource
    applicability: ResourceApplicability
    status: ExecutionStatus
    calibration_mode: CalibrationMode | None = None
    nominal: float | None = None
    lower: float | None = None
    upper: float | None = None
    uncertainty_note: str = ""
    warnings: tuple[str, ...] = ()
    limitations: tuple[str, ...] = ()
    conflicts: tuple[str, ...] = ()
    feature_flag: str = "secretory_resources_round0_enabled"
    baseline_behavior: str = (
        "feature flag off: catalog builder returns no resources for this "
        "target; no existing pipeline references this layer, so behavior "
        "is identical to before Round 0 existed."
    )
    reviewed_by: str = ""
    schema_version: int = SCHEMA_VERSION

    def validate(self) -> None:
        _require_text(self.resource_id, "resource_id")
        if not isinstance(self.category, ResourceCategory):
            raise SecretoryResourceValidationError("category must be a ResourceCategory value.")
        _require_text(self.canonical_unit, "canonical_unit")
        if not isinstance(self.status, ExecutionStatus):
            raise SecretoryResourceValidationError("status must be an ExecutionStatus value.")
        self.source.validate()
        self.applicability.validate()

        numeric_value_set = any(value is not None for value in (self.nominal, self.lower, self.upper))
        if numeric_value_set:
            # Round 0 freezes architecture and status only; no category has a
            # computed value yet. A later round must relax this, not Round 0.
            raise SecretoryResourceValidationError(
                f"{self.resource_id}: Round 0 does not compute nominal/lower/upper values; "
                "leave them None regardless of status."
            )

        if self.status is ExecutionStatus.EXECUTABLE:
            if not self.model_handles:
                raise SecretoryResourceValidationError(
                    f"{self.resource_id}: executable status requires at least one model_handle."
                )
            if self.calibration_mode is not CalibrationMode.RELATIVE_UNCALIBRATED:
                raise SecretoryResourceValidationError(
                    f"{self.resource_id}: Round 0 executable status must use "
                    "calibration_mode=relative_uncalibrated; absolute calibration is not "
                    "authorized this round."
                )
        elif self.status is ExecutionStatus.CONFLICT:
            if not self.conflicts:
                raise SecretoryResourceValidationError(
                    f"{self.resource_id}: conflict status requires at least one conflict entry."
                )
            if self.calibration_mode is not None:
                raise SecretoryResourceValidationError(
                    f"{self.resource_id}: conflict status must leave calibration_mode unset."
                )
        else:
            if self.calibration_mode is not None:
                raise SecretoryResourceValidationError(
                    f"{self.resource_id}: {self.status.value} status must leave calibration_mode unset."
                )
            if self.status in (ExecutionStatus.UNAVAILABLE, ExecutionStatus.NOT_APPLICABLE) and self.model_handles:
                raise SecretoryResourceValidationError(
                    f"{self.resource_id}: {self.status.value} status must not carry model_handles "
                    "(that would silently look executable)."
                )


@dataclass(frozen=True)
class SecretoryResourceCatalog:
    target_id: str
    host: str
    feature_enabled: bool
    resources: tuple[SecretoryResource, ...] = ()
    schema_version: int = SCHEMA_VERSION

    def validate(self) -> None:
        _require_text(self.target_id, "target_id")
        _require_text(self.host, "host")
        if not self.feature_enabled and self.resources:
            raise SecretoryResourceValidationError(
                "feature_enabled=False must produce an empty catalog (baseline/feature-off contract)."
            )
        seen_ids: set[str] = set()
        for resource in self.resources:
            resource.validate()
            if resource.applicability.target_id != self.target_id:
                raise SecretoryResourceValidationError(
                    f"{resource.resource_id}: applicability.target_id must match catalog.target_id "
                    f"({resource.applicability.target_id!r} != {self.target_id!r})."
                )
            if resource.resource_id in seen_ids:
                raise SecretoryResourceValidationError(f"duplicate resource_id: {resource.resource_id}")
            seen_ids.add(resource.resource_id)

    def by_category(self, category: ResourceCategory) -> tuple[SecretoryResource, ...]:
        return tuple(resource for resource in self.resources if resource.category is category)


def _require_text(value: str, field_name: str) -> None:
    if not isinstance(value, str) or not value.strip():
        raise SecretoryResourceValidationError(f"{field_name} must be a non-empty string.")


__all__ = [
    "SCHEMA_VERSION",
    "CalibrationMode",
    "DEFERRED_CATEGORIES",
    "EvidenceClass",
    "ExecutionStatus",
    "IN_SCOPE_CATEGORIES",
    "ResourceApplicability",
    "ResourceCategory",
    "ResourceSource",
    "SecretoryResource",
    "SecretoryResourceCatalog",
    "SecretoryResourceError",
    "SecretoryResourceValidationError",
]
