from __future__ import annotations

import math
from dataclasses import dataclass
from enum import Enum

from pcsec_pichia.errors import OECapacityError, OECapacityValidationError


class OECapacityParameterConflictError(OECapacityValidationError):
    """Raised when equally preferred parameter evidence conflicts."""


class EvidenceSourceType(str, Enum):
    CURRENT_MODEL = "current_model"
    LOCAL_ENZYME_DATA = "local_enzyme_data"
    REVIEWED_PICHIA_MAPPING = "reviewed_pichia_mapping"
    EXTERNAL_PICHIA_MODEL = "external_pichia_model"
    PICHIA_LITERATURE = "pichia_literature"
    HOMOLOGY_TRANSFER = "homology_transfer"
    SMOKE_FIXTURE = "smoke_fixture"


class ConfidenceLevel(str, Enum):
    HIGH = "high"
    MEDIUM = "medium"
    LOW = "low"


class OEDoseMode(str, Enum):
    EXPLICIT_MULTIPLIER = "explicit_multiplier"
    PROMOTER_COPY_MAPPING = "promoter_copy_mapping"
    CATEGORICAL_ONLY = "categorical_only"


class GPRRole(str, Enum):
    SINGLE_GENE = "single_gene"
    ISOENZYME = "isoenzyme"
    COMPLEX_SUBUNIT = "complex_subunit"
    MIXED = "mixed"
    UNRESOLVED = "unresolved"


class OEExecutionStatus(str, Enum):
    GENE_LEVEL_EXECUTABLE = "gene_level_executable"
    PARTIAL_MAPPING = "partial_mapping"
    ISOENZYME_AMBIGUOUS = "isoenzyme_ambiguous"
    COMPLEX_LIMITED = "complex_limited"
    EXTERNAL_EVIDENCE_ONLY = "external_evidence_only"
    CATEGORICAL_DOSE_ONLY = "categorical_dose_only"
    PROXY_ONLY = "proxy_only"
    UNRESOLVED = "unresolved"


class OEExecutionMode(str, Enum):
    GENE_CAPACITY = "gene_capacity"
    RELATIVE_GENE_CAPACITY = "relative_gene_capacity"
    REACTION_PROXY = "reaction_proxy"
    COMPARISON = "comparison"
    NOT_EXECUTABLE = "not_executable"


class OEProductMode(str, Enum):
    REACTION_PROXY = "reaction_proxy"
    RELATIVE_UNCALIBRATED = "relative_uncalibrated"
    ABSOLUTE_CAPACITY = "absolute_capacity"
    NOT_EXECUTABLE = "not_executable"


class OEProductState(str, Enum):
    REACTION_PROXY = "reaction_proxy"
    RELATIVE_UNCALIBRATED = "relative_uncalibrated"
    ABSOLUTE_AVAILABLE = "absolute_available"
    ABSOLUTE_UNAVAILABLE = "absolute_unavailable"
    NOT_EXECUTABLE = "not_executable"


class AbsoluteCapacityAvailability(str, Enum):
    AVAILABLE_REVIEWED = "available_reviewed"
    UNAVAILABLE_MISSING_REVIEWED_ANCHOR = "unavailable_missing_reviewed_anchor"
    UNAVAILABLE_INCOMPATIBLE_ANCHOR = "unavailable_incompatible_anchor"
    NOT_APPLICABLE = "not_applicable"


class OECalibrationStatus(str, Enum):
    PROXY_ONLY = "proxy_only"
    RELATIVE_UNCALIBRATED = "relative_uncalibrated"
    REVIEWED_ABSOLUTE = "reviewed_absolute"
    UNAVAILABLE = "unavailable"
    NOT_APPLICABLE = "not_applicable"


class ParameterScenario(str, Enum):
    LOW = "low"
    NOMINAL = "nominal"
    HIGH = "high"


class ResourceCostMode(str, Enum):
    CURRENT_PROTEIN_POOL = "current_protein_pool"
    FORMATION_DILUTION = "formation_dilution"
    NOT_AVAILABLE = "not_available"


class ConstraintChangeKind(str, Enum):
    ENZYME_CAPACITY_COEFFICIENT = "enzyme_capacity_coefficient"
    FORMATION_DILUTION_BOUND = "formation_dilution_bound"
    PROTEIN_RESOURCE_COEFFICIENT = "protein_resource_coefficient"
    REACTION_BOUND_PROXY = "reaction_bound_proxy"


_STRICTLY_POSITIVE_PARAMETERS = {
    "kcat",
    "molecular_weight",
    "baseline_enzyme_amount",
    "baseline_capacity",
    "complex_stoichiometry",
    "expression_multiplier",
}


@dataclass(frozen=True)
class ParameterEstimate:
    parameter_name: str
    nominal_value: float
    lower_bound: float
    upper_bound: float
    unit: str
    source_type: EvidenceSourceType
    source_ref: str
    source_version: str
    confidence: ConfidenceLevel
    is_transferred: bool = False
    warnings: tuple[str, ...] = ()

    def validate(self) -> None:
        _require_text(self.parameter_name, "parameter_name")
        _require_text(self.unit, "unit")
        _require_text(self.source_ref, "source_ref")
        _require_text(self.source_version, "source_version")
        for name, value in (
            ("nominal_value", self.nominal_value),
            ("lower_bound", self.lower_bound),
            ("upper_bound", self.upper_bound),
        ):
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
                raise OECapacityValidationError(f"{name} must be a finite number.")
        if not self.lower_bound <= self.nominal_value <= self.upper_bound:
            raise OECapacityValidationError(
                "parameter interval must satisfy lower_bound <= nominal_value <= upper_bound."
            )
        if self.parameter_name in _STRICTLY_POSITIVE_PARAMETERS:
            for name, value in (
                ("nominal_value", self.nominal_value),
                ("lower_bound", self.lower_bound),
                ("upper_bound", self.upper_bound),
            ):
                _require_positive_number(value, name)
        if self.is_transferred and self.source_type is not EvidenceSourceType.HOMOLOGY_TRANSFER:
            raise OECapacityValidationError(
                "is_transferred parameters must use source_type=homology_transfer."
            )

    def value_for_scenario(self, scenario: ParameterScenario) -> float:
        if scenario is ParameterScenario.LOW:
            return float(self.lower_bound)
        if scenario is ParameterScenario.HIGH:
            return float(self.upper_bound)
        return float(self.nominal_value)


@dataclass(frozen=True)
class CapacityAnchor:
    anchor_id: str
    target_id: str
    context_id: str
    gene_id: str
    enzyme_id: str
    formation_or_dilution_reaction_id: str
    model_fingerprint: str
    baseline_capacity: float
    unit: str
    source_ref: str
    source_version: str
    reviewed_by: str
    reviewed_at: str

    def validate(self) -> None:
        for field_name, value in (
            ("anchor_id", self.anchor_id),
            ("target_id", self.target_id),
            ("context_id", self.context_id),
            ("gene_id", self.gene_id),
            ("enzyme_id", self.enzyme_id),
            (
                "formation_or_dilution_reaction_id",
                self.formation_or_dilution_reaction_id,
            ),
            ("model_fingerprint", self.model_fingerprint),
            ("source_ref", self.source_ref),
            ("source_version", self.source_version),
            ("reviewed_by", self.reviewed_by),
            ("reviewed_at", self.reviewed_at),
        ):
            _require_text(value, field_name)
        _require_positive_number(self.baseline_capacity, "baseline_capacity")
        if self.unit != "model_flux":
            raise OECapacityValidationError(
                "CapacityAnchor requires canonical unit=model_flux."
            )

    def as_parameter_estimate(self) -> ParameterEstimate:
        self.validate()
        return ParameterEstimate(
            parameter_name="baseline_enzyme_amount",
            nominal_value=float(self.baseline_capacity),
            lower_bound=float(self.baseline_capacity),
            upper_bound=float(self.baseline_capacity),
            unit=self.unit,
            source_type=EvidenceSourceType.REVIEWED_PICHIA_MAPPING,
            source_ref=self.source_ref,
            source_version=self.source_version,
            confidence=ConfidenceLevel.HIGH,
        )

    def binding(
        self,
        *,
        asset_version: str,
        asset_sha256: str,
    ) -> "CapacityAnchorBinding":
        self.validate()
        binding = CapacityAnchorBinding(
            anchor_id=self.anchor_id,
            target_id=self.target_id,
            context_id=self.context_id,
            gene_id=self.gene_id,
            enzyme_id=self.enzyme_id,
            formation_or_dilution_reaction_id=self.formation_or_dilution_reaction_id,
            model_fingerprint=self.model_fingerprint,
            asset_version=asset_version,
            asset_sha256=asset_sha256,
            source_ref=self.source_ref,
            reviewed_by=self.reviewed_by,
            reviewed_at=self.reviewed_at,
        )
        binding.validate()
        return binding


@dataclass(frozen=True)
class CapacityAnchorBinding:
    anchor_id: str
    target_id: str
    context_id: str
    gene_id: str
    enzyme_id: str
    formation_or_dilution_reaction_id: str
    model_fingerprint: str
    asset_version: str
    asset_sha256: str
    source_ref: str
    reviewed_by: str
    reviewed_at: str

    def validate(self) -> None:
        for field_name, value in (
            ("anchor_id", self.anchor_id),
            ("target_id", self.target_id),
            ("context_id", self.context_id),
            ("gene_id", self.gene_id),
            ("enzyme_id", self.enzyme_id),
            (
                "formation_or_dilution_reaction_id",
                self.formation_or_dilution_reaction_id,
            ),
            ("model_fingerprint", self.model_fingerprint),
            ("asset_version", self.asset_version),
            ("source_ref", self.source_ref),
            ("reviewed_by", self.reviewed_by),
            ("reviewed_at", self.reviewed_at),
        ):
            _require_text(value, field_name)
        if len(self.asset_sha256) != 64 or any(
            character not in "0123456789abcdef"
            for character in self.asset_sha256.lower()
        ):
            raise OECapacityValidationError(
                "capacity anchor binding asset_sha256 must be a 64-character hex digest."
            )


@dataclass(frozen=True)
class CapacityAnchorCatalog:
    model_fingerprint: str
    anchors: tuple[CapacityAnchor, ...]
    schema_version: int = 1
    source_ref: str = ""
    asset_version: str = ""
    source_sha256: str = ""

    def validate(self) -> None:
        if self.schema_version != 1:
            raise OECapacityValidationError(
                "CapacityAnchorCatalog requires schema_version=1."
            )
        anchor_ids: set[str] = set()
        identities: set[tuple[str, str, str, str, str, str]] = set()
        for anchor in self.anchors:
            anchor.validate()
            identity = (
                anchor.model_fingerprint,
                anchor.target_id,
                anchor.context_id,
                anchor.gene_id,
                anchor.enzyme_id,
                anchor.formation_or_dilution_reaction_id,
            )
            if anchor.anchor_id in anchor_ids or identity in identities:
                raise OECapacityValidationError("duplicate capacity anchor identity.")
            anchor_ids.add(anchor.anchor_id)
            identities.add(identity)
        if self.anchors:
            _require_text(self.source_ref, "capacity anchor catalog source_ref")
            _require_text(self.asset_version, "capacity anchor catalog asset_version")
            if len(self.source_sha256) != 64 or any(
                character not in "0123456789abcdef"
                for character in self.source_sha256.lower()
            ):
                raise OECapacityValidationError(
                    "capacity anchor catalog source_sha256 must be a 64-character hex digest."
                )


@dataclass(frozen=True)
class OEDoseSpec:
    dose_id: str
    dose_mode: OEDoseMode
    expression_multiplier: float | None = None
    promoter: str = ""
    copy_number: float | None = None
    induction_mode: str = ""
    mapping_source: str = ""
    uncertainty: ParameterEstimate | None = None
    warnings: tuple[str, ...] = ()

    def validate(self) -> None:
        _require_text(self.dose_id, "dose_id")
        if self.expression_multiplier is not None:
            _require_positive_number(self.expression_multiplier, "expression_multiplier")
        if self.copy_number is not None:
            _require_positive_number(self.copy_number, "copy_number")
        if self.uncertainty is not None:
            self.uncertainty.validate()
            if self.uncertainty.parameter_name != "expression_multiplier":
                raise OECapacityValidationError(
                    "dose uncertainty must describe expression_multiplier."
                )
        if self.dose_mode is OEDoseMode.CATEGORICAL_ONLY:
            if self.expression_multiplier is not None:
                raise OECapacityValidationError(
                    "categorical_only dose must not define expression_multiplier."
                )
            if not self.promoter and not self.induction_mode:
                raise OECapacityValidationError(
                    "categorical_only dose requires promoter or induction_mode."
                )
        elif self.dose_mode is OEDoseMode.EXPLICIT_MULTIPLIER:
            if self.expression_multiplier is None:
                raise OECapacityValidationError(
                    "explicit_multiplier dose requires expression_multiplier."
                )
        elif self.dose_mode is OEDoseMode.PROMOTER_COPY_MAPPING:
            if self.expression_multiplier is None or not self.mapping_source:
                raise OECapacityValidationError(
                    "promoter_copy_mapping requires expression_multiplier and mapping_source."
                )


@dataclass(frozen=True)
class GeneEnzymeReactionMapping:
    mapping_id: str
    model_fingerprint: str
    gene_id: str
    enzyme_id: str
    reaction_id: str
    gpr_rule: str
    gpr_role: GPRRole
    mapping_source: EvidenceSourceType
    mapping_confidence: ConfidenceLevel
    execution_status: OEExecutionStatus
    complex_id: str = ""
    subunit_ids: tuple[str, ...] = ()
    subunit_stoichiometry: tuple[tuple[str, float], ...] = ()
    enzyme_variable_id: str = ""
    formation_or_dilution_reaction_id: str = ""
    source_ref: str = ""
    source_version: str = ""
    assembly_evidence_ref: str = ""
    warnings: tuple[str, ...] = ()
    missing_information: tuple[str, ...] = ()

    def validate(self) -> None:
        _require_text(self.mapping_id, "mapping_id")
        _require_text(self.gene_id, "gene_id")
        expected_status = derive_mapping_execution_status(
            role=self.gpr_role,
            mapping_source=self.mapping_source,
            model_mapping_complete=all(
                (
                    self.model_fingerprint,
                    self.enzyme_id,
                    self.reaction_id,
                    self.enzyme_variable_id,
                    self.formation_or_dilution_reaction_id,
                )
            ),
        )
        if self.execution_status is not expected_status:
            raise OECapacityValidationError(
                f"{self.gpr_role.value} mapping must use "
                f"execution_status={expected_status.value}."
            )
        if self.gpr_role is GPRRole.COMPLEX_SUBUNIT:
            if not self.complex_id or not self.subunit_stoichiometry:
                if self.execution_status is OEExecutionStatus.GENE_LEVEL_EXECUTABLE:
                    raise OECapacityValidationError(
                        "complex_subunit mapping requires complex_id and subunit_stoichiometry "
                        "before gene_level_executable."
                    )
            if (
                self.execution_status is OEExecutionStatus.GENE_LEVEL_EXECUTABLE
                and not self.assembly_evidence_ref
            ):
                raise OECapacityValidationError(
                    "complex_subunit gene_level_executable mapping requires "
                    "assembly_evidence_ref."
                )
        for subunit_id, coefficient in self.subunit_stoichiometry:
            _require_text(subunit_id, "subunit_id")
            _require_positive_number(coefficient, "subunit_stoichiometry")
        if self.execution_status is OEExecutionStatus.GENE_LEVEL_EXECUTABLE:
            for field_name, value in (
                ("model_fingerprint", self.model_fingerprint),
                ("enzyme_id", self.enzyme_id),
                ("reaction_id", self.reaction_id),
                ("enzyme_variable_id", self.enzyme_variable_id),
                (
                    "formation_or_dilution_reaction_id",
                    self.formation_or_dilution_reaction_id,
                ),
            ):
                _require_text(value, field_name)
            if self.mapping_source not in {
                EvidenceSourceType.CURRENT_MODEL,
                EvidenceSourceType.LOCAL_ENZYME_DATA,
                EvidenceSourceType.REVIEWED_PICHIA_MAPPING,
            }:
                raise OECapacityValidationError(
                    "gene_level_executable mapping requires a reviewed current-model source."
                )


@dataclass(frozen=True)
class GeneCapacitySpec:
    mapping: GeneEnzymeReactionMapping
    kcat: ParameterEstimate | None
    molecular_weight: ParameterEstimate | None
    baseline_enzyme_amount: ParameterEstimate | None
    complex_stoichiometry: ParameterEstimate | None
    dose: OEDoseSpec
    parameter_scenario: ParameterScenario
    resource_cost_mode: ResourceCostMode
    capacity_anchor_binding: CapacityAnchorBinding | None = None
    warnings: tuple[str, ...] = ()
    missing_information: tuple[str, ...] = ()

    def validate(self) -> None:
        self.mapping.validate()
        self.dose.validate()
        if self.mapping.execution_status is not OEExecutionStatus.GENE_LEVEL_EXECUTABLE:
            raise OECapacityValidationError(
                "GeneCapacitySpec requires a gene_level_executable mapping."
            )
        for estimate in (
            self.kcat,
            self.molecular_weight,
            self.baseline_enzyme_amount,
            self.complex_stoichiometry,
        ):
            if estimate is not None:
                estimate.validate()
        required = {
            "kcat": self.kcat,
            "molecular_weight": self.molecular_weight,
            "baseline_enzyme_amount": self.baseline_enzyme_amount,
        }
        absent = tuple(name for name, value in required.items() if value is None)
        if absent:
            raise OECapacityValidationError(
                "GeneCapacitySpec missing required parameters: " + ", ".join(absent)
            )
        if (
            self.baseline_enzyme_amount is not None
            and self.baseline_enzyme_amount.unit != "model_flux"
        ):
            raise OECapacityValidationError(
                "baseline_enzyme_amount requires canonical unit=model_flux."
            )
        if self.capacity_anchor_binding is None:
            raise OECapacityValidationError(
                "absolute GeneCapacitySpec requires reviewed capacity_anchor_binding."
            )
        self.capacity_anchor_binding.validate()
        binding = self.capacity_anchor_binding
        if (
            binding.gene_id != self.mapping.gene_id
            or binding.enzyme_id != self.mapping.enzyme_id
            or binding.formation_or_dilution_reaction_id
            != self.mapping.formation_or_dilution_reaction_id
            or binding.model_fingerprint != self.mapping.model_fingerprint
        ):
            raise OECapacityValidationError(
                "capacity anchor binding does not match the GeneCapacitySpec mapping."
            )
        if self.mapping.gpr_role is GPRRole.COMPLEX_SUBUNIT and self.complex_stoichiometry is None:
            raise OECapacityValidationError(
                "complex_subunit GeneCapacitySpec requires complex_stoichiometry."
            )
        if self.resource_cost_mode is ResourceCostMode.NOT_AVAILABLE:
            raise OECapacityValidationError(
                "executable GeneCapacitySpec requires an auditable resource cost mode."
            )


@dataclass(frozen=True)
class GeneCapacityParameterSet:
    parameter_set_id: str
    mapping_id: str
    gene_id: str
    enzyme_id: str
    kcat: ParameterEstimate | None
    molecular_weight: ParameterEstimate | None
    baseline_enzyme_amount: ParameterEstimate | None
    complex_stoichiometry: ParameterEstimate | None = None
    capacity_anchor_binding: CapacityAnchorBinding | None = None
    warnings: tuple[str, ...] = ()

    def validate(self) -> None:
        _require_text(self.parameter_set_id, "parameter_set_id")
        _require_text(self.gene_id, "gene_id")
        _require_text(self.enzyme_id, "enzyme_id")
        if not self.mapping_id:
            raise OECapacityValidationError("mapping_id must be non-empty.")
        expected_names = (
            ("kcat", self.kcat),
            ("molecular_weight", self.molecular_weight),
            ("baseline_enzyme_amount", self.baseline_enzyme_amount),
            ("complex_stoichiometry", self.complex_stoichiometry),
        )
        for expected_name, estimate in expected_names:
            if estimate is None:
                continue
            estimate.validate()
            if estimate.parameter_name != expected_name:
                raise OECapacityValidationError(
                    f"{expected_name} estimate must use parameter_name={expected_name}."
                )
        if self.capacity_anchor_binding is not None:
            self.capacity_anchor_binding.validate()
        if self.baseline_enzyme_amount is None and self.capacity_anchor_binding is not None:
            raise OECapacityValidationError(
                "capacity_anchor_binding requires baseline_enzyme_amount."
            )

    @property
    def missing_information(self) -> tuple[str, ...]:
        values = (
            ("kcat", self.kcat),
            ("molecular_weight", self.molecular_weight),
            ("baseline_enzyme_amount", self.baseline_enzyme_amount),
        )
        missing = [name for name, value in values if value is None]
        if self.baseline_enzyme_amount is not None and self.capacity_anchor_binding is None:
            missing.append("reviewed_baseline_capacity")
        return tuple(missing)

    @property
    def uses_smoke_fixture(self) -> bool:
        return any(
            estimate is not None
            and estimate.source_type is EvidenceSourceType.SMOKE_FIXTURE
            for estimate in (
                self.kcat,
                self.molecular_weight,
                self.baseline_enzyme_amount,
                self.complex_stoichiometry,
            )
        )


@dataclass(frozen=True)
class ParameterPolicy:
    parameter_sets: tuple[GeneCapacityParameterSet, ...]
    scenarios: tuple[ParameterScenario, ...] = (
        ParameterScenario.LOW,
        ParameterScenario.NOMINAL,
        ParameterScenario.HIGH,
    )
    strict_conflicts: bool = True
    test_only_allow_smoke_fixture: bool = False

    def validate(self) -> None:
        if not self.scenarios or len(set(self.scenarios)) != len(self.scenarios):
            raise OECapacityValidationError(
                "ParameterPolicy scenarios must be non-empty and unique."
            )
        ids: set[str] = set()
        for parameter_set in self.parameter_sets:
            parameter_set.validate()
            if (
                not self.test_only_allow_smoke_fixture
                and parameter_set.uses_smoke_fixture
            ):
                raise OECapacityValidationError(
                    "smoke_fixture parameters are test-only and rejected by production policy."
                )
            if parameter_set.parameter_set_id in ids:
                raise OECapacityValidationError(
                    f"duplicate parameter_set_id: {parameter_set.parameter_set_id}"
                )
            ids.add(parameter_set.parameter_set_id)


@dataclass(frozen=True)
class RelativeOEScenarioSpec:
    mapping: GeneEnzymeReactionMapping
    dose: OEDoseSpec
    parameter_scenario: ParameterScenario
    relative_capacity_factor: float
    kcat: ParameterEstimate | None = None
    molecular_weight: ParameterEstimate | None = None
    parameter_sources: tuple[str, ...] = ()
    warnings: tuple[str, ...] = ()
    limitations: tuple[str, ...] = ()

    def validate(self) -> None:
        self.mapping.validate()
        self.dose.validate()
        if self.mapping.execution_status is not OEExecutionStatus.GENE_LEVEL_EXECUTABLE:
            raise OECapacityValidationError(
                "RelativeOEScenarioSpec requires a complete current-model mapping."
            )
        if self.dose.expression_multiplier is None:
            raise OECapacityValidationError(
                "RelativeOEScenarioSpec requires a numeric dimensionless dose."
            )
        _require_positive_number(
            self.relative_capacity_factor,
            "relative_capacity_factor",
        )
        for estimate in (self.kcat, self.molecular_weight):
            if estimate is not None:
                estimate.validate()
        if not self.parameter_sources:
            raise OECapacityValidationError(
                "RelativeOEScenarioSpec requires auditable parameter_sources."
            )
        if not self.limitations:
            raise OECapacityValidationError(
                "RelativeOEScenarioSpec requires explicit limitations."
            )


@dataclass(frozen=True)
class OECapacityPlan:
    gene_id: str
    target_id: str
    context_id: str
    requested_dose: OEDoseSpec
    execution_mode: OEExecutionMode
    execution_status: OEExecutionStatus
    executable_capacity_specs: tuple[GeneCapacitySpec, ...] = ()
    explain_only_mappings: tuple[GeneEnzymeReactionMapping, ...] = ()
    proxy_reaction_ids: tuple[str, ...] = ()
    constraint_changes: tuple[CapacityConstraintChange, ...] = ()
    uncertainty_scenarios: tuple[ParameterScenario, ...] = ()
    missing_information: tuple[str, ...] = ()
    warnings: tuple[str, ...] = ()
    structural_mappings: tuple[GeneEnzymeReactionMapping, ...] = ()
    available_relative_scenario_specs: tuple[RelativeOEScenarioSpec, ...] = ()
    relative_scenario_specs: tuple[RelativeOEScenarioSpec, ...] = ()
    product_mode: OEProductMode = OEProductMode.RELATIVE_UNCALIBRATED
    product_state: OEProductState = OEProductState.NOT_EXECUTABLE
    absolute_capacity_availability: AbsoluteCapacityAvailability = (
        AbsoluteCapacityAvailability.UNAVAILABLE_MISSING_REVIEWED_ANCHOR
    )
    calibration_status: OECalibrationStatus = OECalibrationStatus.UNAVAILABLE
    absolute_solver_allowed: bool = False
    model_fingerprint: str = ""
    limitations: tuple[str, ...] = ()

    def validate(self) -> None:
        _require_text(self.gene_id, "gene_id")
        _require_text(self.target_id, "target_id")
        _require_text(self.context_id, "context_id")
        self.requested_dose.validate()
        for spec in self.executable_capacity_specs:
            spec.validate()
            if spec.mapping.gene_id != self.gene_id:
                raise OECapacityValidationError(
                    "all executable capacity specs must match plan gene_id."
                )
        for mapping in self.explain_only_mappings:
            mapping.validate()
        for mapping in self.structural_mappings:
            mapping.validate()
        for spec in self.available_relative_scenario_specs:
            spec.validate()
            if spec.mapping.gene_id != self.gene_id:
                raise OECapacityValidationError(
                    "available relative scenario specs must match plan gene_id."
                )
        for spec in self.relative_scenario_specs:
            spec.validate()
            if spec.mapping.gene_id != self.gene_id:
                raise OECapacityValidationError(
                    "relative scenario specs must match plan gene_id."
                )
        _require_text(self.model_fingerprint, "model_fingerprint")
        if not self.limitations:
            raise OECapacityValidationError("OE product plan requires explicit limitations.")
        if self.execution_mode is OEExecutionMode.GENE_CAPACITY:
            if self.execution_status is OEExecutionStatus.PROXY_ONLY:
                raise OECapacityValidationError(
                    "proxy_only status must use execution_mode=reaction_proxy, not gene_capacity."
                )
            if not self.executable_capacity_specs:
                raise OECapacityValidationError(
                    "gene_capacity mode requires executable_capacity_specs."
                )
            if self.requested_dose.dose_mode is OEDoseMode.CATEGORICAL_ONLY:
                raise OECapacityValidationError(
                    "categorical_only dose cannot execute gene_capacity mode."
                )
        elif self.execution_mode is OEExecutionMode.REACTION_PROXY:
            if not self.proxy_reaction_ids:
                raise OECapacityValidationError(
                    "reaction_proxy mode requires proxy_reaction_ids."
                )
            if (
                self.product_state is OEProductState.REACTION_PROXY
                and self.execution_status is not OEExecutionStatus.PROXY_ONLY
            ):
                raise OECapacityValidationError(
                    "reaction_proxy mode must use execution_status=proxy_only."
                )
        elif self.execution_mode is OEExecutionMode.RELATIVE_GENE_CAPACITY:
            if not self.relative_scenario_specs:
                raise OECapacityValidationError(
                    "relative_gene_capacity requires relative scenario specs."
                )
            if self.executable_capacity_specs or self.absolute_solver_allowed:
                raise OECapacityValidationError(
                    "relative_gene_capacity cannot contain absolute capacity specs."
                )
        elif self.execution_mode is OEExecutionMode.COMPARISON:
            if not self.executable_capacity_specs or not self.proxy_reaction_ids:
                raise OECapacityValidationError(
                    "comparison mode requires both executable_capacity_specs and "
                    "proxy_reaction_ids."
                )
            if self.requested_dose.dose_mode is OEDoseMode.CATEGORICAL_ONLY:
                raise OECapacityValidationError(
                    "categorical_only dose cannot execute comparison mode."
                )
        elif self.execution_mode is OEExecutionMode.NOT_EXECUTABLE:
            if self.executable_capacity_specs or self.constraint_changes:
                raise OECapacityValidationError(
                    "not_executable plan cannot contain executable specs or constraint changes."
                )
        if self.absolute_solver_allowed:
            if self.product_state is not OEProductState.ABSOLUTE_AVAILABLE:
                raise OECapacityValidationError(
                    "absolute_solver_allowed requires product_state=absolute_available."
                )
            if not self.executable_capacity_specs:
                raise OECapacityValidationError(
                    "absolute_solver_allowed requires executable capacity specs."
                )
            for spec in self.executable_capacity_specs:
                binding = spec.capacity_anchor_binding
                if binding is None:
                    raise OECapacityValidationError(
                        "absolute solver requires reviewed capacity anchor bindings."
                    )
                if (
                    binding.target_id != self.target_id
                    or binding.context_id != self.context_id
                    or binding.model_fingerprint != self.model_fingerprint
                ):
                    raise OECapacityValidationError(
                        "capacity anchor binding does not match plan target/context/model."
                    )
        elif self.executable_capacity_specs:
            raise OECapacityValidationError(
                "executable capacity specs require absolute_solver_allowed=True."
            )
        if self.product_state is OEProductState.RELATIVE_UNCALIBRATED:
            if not self.relative_scenario_specs:
                raise OECapacityValidationError(
                    "relative_uncalibrated requires relative scenario specs."
                )
            if self.calibration_status is not OECalibrationStatus.RELATIVE_UNCALIBRATED:
                raise OECapacityValidationError(
                    "relative_uncalibrated requires matching calibration_status."
                )
            if self.executable_capacity_specs or self.absolute_solver_allowed:
                raise OECapacityValidationError(
                    "relative_uncalibrated cannot contain absolute executable specs."
                )
            if self.execution_mode is not OEExecutionMode.RELATIVE_GENE_CAPACITY:
                raise OECapacityValidationError(
                    "relative_uncalibrated requires relative_gene_capacity execution."
                )
            if self.product_mode is not OEProductMode.RELATIVE_UNCALIBRATED:
                raise OECapacityValidationError(
                    "relative_uncalibrated requires matching product_mode."
                )
        if self.product_state is OEProductState.REACTION_PROXY:
            if self.relative_scenario_specs:
                raise OECapacityValidationError(
                    "reaction_proxy cannot retain relative scenario specs."
                )
            if self.execution_mode is not OEExecutionMode.REACTION_PROXY:
                raise OECapacityValidationError(
                    "reaction_proxy product state requires reaction_proxy execution."
                )
            if self.calibration_status is not OECalibrationStatus.PROXY_ONLY:
                raise OECapacityValidationError(
                    "reaction_proxy requires proxy_only calibration status."
                )
        if self.product_state is OEProductState.ABSOLUTE_AVAILABLE:
            if self.relative_scenario_specs:
                raise OECapacityValidationError(
                    "absolute_available cannot retain relative scenario specs."
                )
            if not self.absolute_solver_allowed:
                raise OECapacityValidationError(
                    "absolute_available requires the reviewed absolute solver gate."
                )
            if self.calibration_status is not OECalibrationStatus.REVIEWED_ABSOLUTE:
                raise OECapacityValidationError(
                    "absolute_available requires reviewed_absolute calibration status."
                )
            if self.absolute_capacity_availability is not AbsoluteCapacityAvailability.AVAILABLE_REVIEWED:
                raise OECapacityValidationError(
                    "absolute_available requires available_reviewed availability."
                )
            if self.product_mode is not OEProductMode.ABSOLUTE_CAPACITY:
                raise OECapacityValidationError(
                    "absolute_available requires absolute_capacity product_mode."
                )
        if self.product_state is OEProductState.ABSOLUTE_UNAVAILABLE:
            if self.relative_scenario_specs:
                raise OECapacityValidationError(
                    "absolute_unavailable cannot retain relative scenario specs."
                )
            if self.execution_mode is not OEExecutionMode.NOT_EXECUTABLE:
                raise OECapacityValidationError(
                    "absolute_unavailable must not execute a solver mode."
                )
            if self.absolute_solver_allowed:
                raise OECapacityValidationError(
                    "absolute_unavailable cannot allow the absolute solver."
                )
            if self.absolute_capacity_availability is AbsoluteCapacityAvailability.AVAILABLE_REVIEWED:
                raise OECapacityValidationError(
                    "absolute_unavailable requires an unavailable capacity state."
                )
            if self.calibration_status is not OECalibrationStatus.UNAVAILABLE:
                raise OECapacityValidationError(
                    "absolute_unavailable requires unavailable calibration status."
                )
            if self.product_mode is not OEProductMode.ABSOLUTE_CAPACITY:
                raise OECapacityValidationError(
                    "absolute_unavailable requires absolute_capacity product_mode."
                )
        if self.product_state is OEProductState.NOT_EXECUTABLE:
            if self.relative_scenario_specs:
                raise OECapacityValidationError(
                    "not_executable cannot retain relative scenario specs."
                )
            if self.execution_mode is not OEExecutionMode.NOT_EXECUTABLE:
                raise OECapacityValidationError(
                    "not_executable product state requires not_executable execution."
                )
            if self.product_mode is not OEProductMode.NOT_EXECUTABLE:
                raise OECapacityValidationError(
                    "not_executable product state requires matching product_mode."
                )


@dataclass(frozen=True)
class GeneCapacityCatalog:
    model_fingerprint: str
    mappings: tuple[GeneEnzymeReactionMapping, ...]
    source_refs: tuple[str, ...] = ()
    warnings: tuple[str, ...] = ()

    def validate(self) -> None:
        _require_text(self.model_fingerprint, "model_fingerprint")
        mapping_ids: set[str] = set()
        for mapping in self.mappings:
            mapping.validate()
            if mapping.mapping_id in mapping_ids:
                raise OECapacityValidationError(
                    f"duplicate mapping_id: {mapping.mapping_id}"
                )
            mapping_ids.add(mapping.mapping_id)
            if mapping.model_fingerprint and mapping.model_fingerprint != self.model_fingerprint:
                raise OECapacityValidationError(
                    "all current-model mappings must match catalog model_fingerprint."
                )


@dataclass(frozen=True)
class GeneCapacityCoverage:
    total_mappings: int
    gene_count: int
    reaction_count: int
    enzyme_count: int
    by_role: tuple[tuple[str, int], ...]
    by_status: tuple[tuple[str, int], ...]

    def validate(self) -> None:
        for field_name, value in (
            ("total_mappings", self.total_mappings),
            ("gene_count", self.gene_count),
            ("reaction_count", self.reaction_count),
            ("enzyme_count", self.enzyme_count),
        ):
            if not isinstance(value, int) or isinstance(value, bool) or value < 0:
                raise OECapacityValidationError(
                    f"{field_name} must be a non-negative integer."
                )


@dataclass(frozen=True)
class GeneCapacityValidationIssue:
    code: str
    message: str
    severity: str
    mapping_id: str = ""
    field_name: str = ""

    def validate(self) -> None:
        _require_text(self.code, "code")
        _require_text(self.message, "message")
        if self.severity not in {"error", "warning"}:
            raise OECapacityValidationError("severity must be error or warning.")


@dataclass(frozen=True)
class GeneCapacityValidationResult:
    issues: tuple[GeneCapacityValidationIssue, ...] = ()

    @property
    def errors(self) -> tuple[GeneCapacityValidationIssue, ...]:
        return tuple(issue for issue in self.issues if issue.severity == "error")

    @property
    def warnings(self) -> tuple[GeneCapacityValidationIssue, ...]:
        return tuple(issue for issue in self.issues if issue.severity == "warning")

    @property
    def is_valid(self) -> bool:
        return not self.errors

    def validate(self) -> None:
        for issue in self.issues:
            issue.validate()


@dataclass(frozen=True)
class CapacityConstraintChange:
    change_id: str
    scenario: ParameterScenario
    change_kind: ConstraintChangeKind
    constraint_block: str
    variable_id: str
    reaction_id: str
    old_value: float
    new_value: float
    unit: str
    source_ref: str
    resource_cost_mode: ResourceCostMode
    metadata: tuple[tuple[str, str], ...] = ()
    warnings: tuple[str, ...] = ()

    def validate(self) -> None:
        for field_name, value in (
            ("change_id", self.change_id),
            ("constraint_block", self.constraint_block),
            ("variable_id", self.variable_id),
            ("unit", self.unit),
            ("source_ref", self.source_ref),
        ):
            _require_text(value, field_name)
        for field_name, value in (("old_value", self.old_value), ("new_value", self.new_value)):
            if (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not math.isfinite(value)
            ):
                raise OECapacityValidationError(f"{field_name} must be a finite number.")


@dataclass(frozen=True)
class OECapacityConstraintBundle:
    model_fingerprint: str
    plan: OECapacityPlan
    changes: tuple[CapacityConstraintChange, ...]
    warnings: tuple[str, ...] = ()

    def validate(self) -> None:
        _require_text(self.model_fingerprint, "model_fingerprint")
        self.plan.validate()
        if self.changes and self.plan.execution_mode is OEExecutionMode.NOT_EXECUTABLE:
            raise OECapacityValidationError(
                "not_executable plan cannot produce constraint changes."
            )
        change_ids: set[str] = set()
        for change in self.changes:
            change.validate()
            if change.change_id in change_ids:
                raise OECapacityValidationError(
                    f"duplicate constraint change_id: {change.change_id}"
                )
            change_ids.add(change.change_id)
            if (
                change.change_kind is ConstraintChangeKind.REACTION_BOUND_PROXY
                and self.plan.execution_mode is not OEExecutionMode.REACTION_PROXY
            ):
                raise OECapacityValidationError(
                    "reaction_bound_proxy changes are only valid for reaction_proxy plans."
                )


@dataclass(frozen=True)
class SolverSnapshot:
    execution_mode: OEExecutionMode
    backend: str
    solver_status: str
    success: bool
    secretion_objective: float | None
    growth_retention: float | None
    max_feasible_growth_rate: float | None
    protein_resource_cost: float | None
    parameter_scenario: ParameterScenario | None = None
    constraint_counts: tuple[tuple[str, int], ...] = ()
    key_fluxes: tuple[tuple[str, float], ...] = ()
    message: str = ""
    warnings: tuple[str, ...] = ()
    attempt_id: str = ""

    def validate(self) -> None:
        _require_text(self.backend, "backend")
        _require_text(self.solver_status, "solver_status")
        for field_name, value in (
            ("secretion_objective", self.secretion_objective),
            ("growth_retention", self.growth_retention),
            ("max_feasible_growth_rate", self.max_feasible_growth_rate),
            ("protein_resource_cost", self.protein_resource_cost),
        ):
            if value is not None and (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not math.isfinite(value)
            ):
                raise OECapacityValidationError(f"{field_name} must be finite when present.")
        for name, count in self.constraint_counts:
            _require_text(name, "constraint_count name")
            if not isinstance(count, int) or isinstance(count, bool) or count < 0:
                raise OECapacityValidationError(
                    "constraint counts must be non-negative integers."
                )


@dataclass(frozen=True)
class OECapacityScenarioResult:
    parameter_scenario: ParameterScenario
    baseline: SolverSnapshot
    perturbed: SolverSnapshot
    objective_delta: float | None = None
    protein_resource_cost_delta: float | None = None
    failure_reason: str = ""

    def validate(self) -> None:
        self.baseline.validate()
        self.perturbed.validate()
        if self.perturbed.execution_mode not in {
            OEExecutionMode.GENE_CAPACITY,
            OEExecutionMode.RELATIVE_GENE_CAPACITY,
        }:
            raise OECapacityValidationError(
                "scenario perturbed snapshot must use a gene-capacity execution mode."
            )
        if self.baseline.parameter_scenario is not self.parameter_scenario:
            raise OECapacityValidationError(
                "scenario result must match baseline parameter_scenario."
            )
        if self.perturbed.parameter_scenario is not self.parameter_scenario:
            raise OECapacityValidationError(
                "scenario result must match perturbed parameter_scenario."
            )
        for field_name, value in (
            ("objective_delta", self.objective_delta),
            ("protein_resource_cost_delta", self.protein_resource_cost_delta),
        ):
            if value is not None and (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not math.isfinite(value)
            ):
                raise OECapacityValidationError(
                    f"{field_name} must be finite when present."
                )
        if (not self.baseline.success or not self.perturbed.success) and not self.failure_reason:
            raise OECapacityValidationError(
                "failed scenario result requires failure_reason."
            )


@dataclass(frozen=True)
class OECapacityComparisonResult:
    gene_id: str
    target_id: str
    context_id: str
    execution_status: OEExecutionStatus
    baseline: SolverSnapshot
    proxy: SolverSnapshot | None
    gene_capacity_scenarios: tuple[SolverSnapshot, ...]
    gene_capacity_vs_baseline_delta: float | None
    gene_capacity_vs_proxy_delta: float | None
    protein_resource_cost_delta: float | None
    skipped_reason: str = ""
    missing_information: tuple[str, ...] = ()
    traceability: tuple[tuple[str, str], ...] = ()
    warnings: tuple[str, ...] = ()
    scenario_results: tuple[OECapacityScenarioResult, ...] = ()
    proxy_attempts: tuple[SolverSnapshot, ...] = ()
    relative_scenarios: tuple[SolverSnapshot, ...] = ()
    relative_scenario_results: tuple[OECapacityScenarioResult, ...] = ()
    relative_vs_baseline_delta: float | None = None
    relative_vs_proxy_delta: float | None = None

    def validate(self) -> None:
        _require_text(self.gene_id, "gene_id")
        _require_text(self.target_id, "target_id")
        _require_text(self.context_id, "context_id")
        self.baseline.validate()
        if self.proxy is not None:
            self.proxy.validate()
            if self.proxy.execution_mode is not OEExecutionMode.REACTION_PROXY:
                raise OECapacityValidationError(
                    "proxy snapshot must use execution_mode=reaction_proxy."
                )
        seen_scenarios: set[ParameterScenario] = set()
        for snapshot in self.gene_capacity_scenarios:
            snapshot.validate()
            if snapshot.execution_mode is not OEExecutionMode.GENE_CAPACITY:
                raise OECapacityValidationError(
                    "gene capacity snapshots must use execution_mode=gene_capacity."
                )
            if snapshot.parameter_scenario is None:
                raise OECapacityValidationError(
                    "gene capacity snapshots require parameter_scenario."
                )
            if snapshot.parameter_scenario in seen_scenarios:
                raise OECapacityValidationError(
                    f"duplicate parameter scenario: {snapshot.parameter_scenario.value}"
                )
            seen_scenarios.add(snapshot.parameter_scenario)
        for field_name, value in (
            ("gene_capacity_vs_baseline_delta", self.gene_capacity_vs_baseline_delta),
            ("gene_capacity_vs_proxy_delta", self.gene_capacity_vs_proxy_delta),
            ("protein_resource_cost_delta", self.protein_resource_cost_delta),
        ):
            if value is not None and (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not math.isfinite(value)
            ):
                raise OECapacityValidationError(f"{field_name} must be finite when present.")
        scenario_ids: set[ParameterScenario] = set()
        for result in self.scenario_results:
            result.validate()
            if result.parameter_scenario in scenario_ids:
                raise OECapacityValidationError(
                    f"duplicate scenario result: {result.parameter_scenario.value}"
                )
            scenario_ids.add(result.parameter_scenario)
        attempt_ids: set[str] = set()
        for snapshot in self.proxy_attempts:
            snapshot.validate()
            if snapshot.execution_mode is not OEExecutionMode.REACTION_PROXY:
                raise OECapacityValidationError(
                    "proxy attempts must use execution_mode=reaction_proxy."
                )
            _require_text(snapshot.attempt_id, "proxy attempt_id")
            if snapshot.attempt_id in attempt_ids:
                raise OECapacityValidationError(
                    f"duplicate proxy attempt_id: {snapshot.attempt_id}"
                )
            attempt_ids.add(snapshot.attempt_id)
        relative_ids: set[ParameterScenario] = set()
        for snapshot in self.relative_scenarios:
            snapshot.validate()
            if snapshot.execution_mode is not OEExecutionMode.RELATIVE_GENE_CAPACITY:
                raise OECapacityValidationError(
                    "relative snapshots must use relative_gene_capacity execution."
                )
        for result in self.relative_scenario_results:
            result.validate()
            if result.perturbed.execution_mode is not OEExecutionMode.RELATIVE_GENE_CAPACITY:
                raise OECapacityValidationError(
                    "relative scenario results require relative_gene_capacity snapshots."
                )
            if result.parameter_scenario in relative_ids:
                raise OECapacityValidationError(
                    f"duplicate relative scenario: {result.parameter_scenario.value}"
                )
            relative_ids.add(result.parameter_scenario)
        for field_name, value in (
            ("relative_vs_baseline_delta", self.relative_vs_baseline_delta),
            ("relative_vs_proxy_delta", self.relative_vs_proxy_delta),
        ):
            if value is not None and (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not math.isfinite(value)
            ):
                raise OECapacityValidationError(f"{field_name} must be finite when present.")


@dataclass(frozen=True)
class OECapacityScreenRequest:
    gene_id: str
    target_id: str
    context_id: str
    dose: OEDoseSpec
    execution_mode: OEExecutionMode
    product_mode: OEProductMode | None = None

    def validate(self) -> None:
        _require_text(self.gene_id, "gene_id")
        _require_text(self.target_id, "target_id")
        _require_text(self.context_id, "context_id")
        self.dose.validate()
        if self.product_mode is not None and not isinstance(
            self.product_mode, OEProductMode
        ):
            raise OECapacityValidationError("product_mode must be an OEProductMode value.")
        # Categorical dose is retained as an auditable not-executable product
        # result instead of being rejected at the facade boundary.


@dataclass(frozen=True)
class OECapacityScreenConfig:
    feature_enabled: bool
    compare_proxy: bool
    parameter_scenarios: tuple[ParameterScenario, ...]
    growth_rate: float
    solver_options: tuple[tuple[str, str], ...] = ()

    def validate(self) -> None:
        if not isinstance(self.feature_enabled, bool):
            raise OECapacityValidationError("feature_enabled must be boolean.")
        if not isinstance(self.compare_proxy, bool):
            raise OECapacityValidationError("compare_proxy must be boolean.")
        _require_positive_number(self.growth_rate, "growth_rate")
        if not self.parameter_scenarios:
            raise OECapacityValidationError(
                "parameter_scenarios must contain at least one scenario."
            )
        if len(set(self.parameter_scenarios)) != len(self.parameter_scenarios):
            raise OECapacityValidationError("parameter_scenarios must be unique.")
        option_names: set[str] = set()
        for name, value in self.solver_options:
            _require_text(name, "solver option name")
            _require_text(value, "solver option value")
            if name in option_names:
                raise OECapacityValidationError(f"duplicate solver option: {name}")
            option_names.add(name)


@dataclass(frozen=True)
class OECapacityScreenRow:
    gene_id: str
    target_id: str
    context_id: str
    execution_mode: OEExecutionMode
    execution_status: OEExecutionStatus
    product_mode: OEProductMode
    product_state: OEProductState
    absolute_capacity_availability: AbsoluteCapacityAvailability
    calibration_status: OECalibrationStatus
    absolute_solver_allowed: bool
    model_fingerprint: str
    dose_id: str
    dose_mode: OEDoseMode
    expression_multiplier: float | None
    mapping_ids: tuple[str, ...]
    parameter_sources: tuple[str, ...]
    parameter_confidence: ConfidenceLevel | None
    uncertainty_scenarios: tuple[ParameterScenario, ...]
    baseline_objective: float | None
    proxy_objective: float | None
    gene_capacity_objective: float | None
    gene_capacity_vs_baseline_delta: float | None
    gene_capacity_vs_proxy_delta: float | None
    protein_resource_cost_delta: float | None
    missing_information: tuple[str, ...] = ()
    warnings: tuple[str, ...] = ()
    screen_status: str = "completed"
    scenario_results: tuple[OECapacityScenarioResult, ...] = ()
    proxy_attempts: tuple[SolverSnapshot, ...] = ()
    summary_source: str = ""
    mapping_sources: tuple[str, ...] = ()
    dose_source: str = ""
    relative_capacity_factors: tuple[tuple[ParameterScenario, float], ...] = ()
    nominal_capacity: float | None = None
    nominal_capacities: tuple[tuple[str, float], ...] = ()
    limitations: tuple[str, ...] = ()
    relative_objective: float | None = None
    relative_vs_baseline_delta: float | None = None
    relative_vs_proxy_delta: float | None = None
    relative_scenario_results: tuple[OECapacityScenarioResult, ...] = ()

    def validate(self) -> None:
        for field_name, value in (
            ("gene_id", self.gene_id),
            ("target_id", self.target_id),
            ("context_id", self.context_id),
            ("dose_id", self.dose_id),
            ("model_fingerprint", self.model_fingerprint),
        ):
            _require_text(value, field_name)
        if self.expression_multiplier is not None:
            _require_positive_number(self.expression_multiplier, "expression_multiplier")
        if len(set(self.mapping_ids)) != len(self.mapping_ids):
            raise OECapacityValidationError("mapping_ids must be unique.")
        if len(set(self.parameter_sources)) != len(self.parameter_sources):
            raise OECapacityValidationError("parameter_sources must be unique.")
        if len(set(self.uncertainty_scenarios)) != len(self.uncertainty_scenarios):
            raise OECapacityValidationError("uncertainty_scenarios must be unique.")
        for field_name, value in (
            ("baseline_objective", self.baseline_objective),
            ("proxy_objective", self.proxy_objective),
            ("gene_capacity_objective", self.gene_capacity_objective),
            (
                "gene_capacity_vs_baseline_delta",
                self.gene_capacity_vs_baseline_delta,
            ),
            ("gene_capacity_vs_proxy_delta", self.gene_capacity_vs_proxy_delta),
            ("protein_resource_cost_delta", self.protein_resource_cost_delta),
            ("relative_objective", self.relative_objective),
            ("relative_vs_baseline_delta", self.relative_vs_baseline_delta),
            ("relative_vs_proxy_delta", self.relative_vs_proxy_delta),
        ):
            if value is not None and (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not math.isfinite(value)
            ):
                raise OECapacityValidationError(
                    f"{field_name} must be finite when present."
                )
        if self.screen_status not in {
            "completed",
            "partial_failure",
            "failed",
            "skipped",
            "not_executable",
        }:
            raise OECapacityValidationError(
                "screen_status must be completed, partial_failure, failed, skipped, "
                "or not_executable."
            )
        for result in self.scenario_results:
            result.validate()
        for snapshot in self.proxy_attempts:
            snapshot.validate()
        for result in self.relative_scenario_results:
            result.validate()
        if not self.limitations:
            raise OECapacityValidationError("screen row requires explicit limitations.")
        for scenario, factor in self.relative_capacity_factors:
            if not isinstance(scenario, ParameterScenario):
                raise OECapacityValidationError(
                    "relative capacity factor scenario must be ParameterScenario."
                )
            _require_positive_number(factor, "relative_capacity_factor")
        for handle, value in self.nominal_capacities:
            _require_text(handle, "nominal capacity handle")
            _require_positive_number(value, "nominal capacity value")
        if self.product_state is OEProductState.RELATIVE_UNCALIBRATED:
            if self.absolute_solver_allowed or self.nominal_capacity is not None:
                raise OECapacityValidationError(
                    "relative_uncalibrated cannot expose absolute capacity."
                )
            if self.gene_capacity_objective is not None or self.scenario_results:
                raise OECapacityValidationError(
                    "relative_uncalibrated cannot contain absolute scenario solver output."
                )
            if not self.relative_capacity_factors:
                raise OECapacityValidationError(
                    "relative_uncalibrated requires relative capacity factors."
                )
            if not self.relative_scenario_results:
                raise OECapacityValidationError(
                    "relative_uncalibrated requires relative scenario solver evidence."
                )
        if self.product_state is OEProductState.ABSOLUTE_UNAVAILABLE:
            if self.execution_mode is not OEExecutionMode.NOT_EXECUTABLE:
                raise OECapacityValidationError(
                    "absolute_unavailable row must be not_executable."
                )
            if any(
                value is not None
                for value in (
                    self.baseline_objective,
                    self.proxy_objective,
                    self.gene_capacity_objective,
                    self.nominal_capacity,
                )
            ):
                raise OECapacityValidationError(
                    "absolute_unavailable row must not contain solver objectives or nominal capacity."
                )
        if self.product_state is not OEProductState.RELATIVE_UNCALIBRATED:
            if (
                self.relative_capacity_factors
                or self.relative_scenario_results
                or self.relative_objective is not None
                or self.relative_vs_baseline_delta is not None
                or self.relative_vs_proxy_delta is not None
            ):
                raise OECapacityValidationError(
                    "non-relative product rows cannot contain relative solver evidence."
                )
        if self.product_state is OEProductState.NOT_EXECUTABLE:
            if self.execution_mode is not OEExecutionMode.NOT_EXECUTABLE:
                raise OECapacityValidationError(
                    "not_executable product state requires execution_mode=not_executable."
                )


def derive_mapping_execution_status(
    *,
    role: GPRRole,
    mapping_source: EvidenceSourceType,
    model_mapping_complete: bool,
) -> OEExecutionStatus:
    if mapping_source in {
        EvidenceSourceType.EXTERNAL_PICHIA_MODEL,
        EvidenceSourceType.PICHIA_LITERATURE,
        EvidenceSourceType.HOMOLOGY_TRANSFER,
        EvidenceSourceType.SMOKE_FIXTURE,
    }:
        return OEExecutionStatus.EXTERNAL_EVIDENCE_ONLY
    if role is GPRRole.ISOENZYME:
        return OEExecutionStatus.ISOENZYME_AMBIGUOUS
    if role is GPRRole.COMPLEX_SUBUNIT:
        return OEExecutionStatus.COMPLEX_LIMITED
    if role is GPRRole.MIXED:
        return OEExecutionStatus.PARTIAL_MAPPING
    if role is GPRRole.UNRESOLVED:
        return OEExecutionStatus.UNRESOLVED
    if not model_mapping_complete:
        return OEExecutionStatus.PARTIAL_MAPPING
    return OEExecutionStatus.GENE_LEVEL_EXECUTABLE


@dataclass(frozen=True)
class OECapacityScreenResult:
    model_fingerprint: str
    config: OECapacityScreenConfig
    rows: tuple[OECapacityScreenRow, ...]
    failures: tuple[OECapacityScreenRow, ...] = ()
    warnings: tuple[str, ...] = ()

    def validate(self) -> None:
        _require_text(self.model_fingerprint, "model_fingerprint")
        self.config.validate()
        for row in (*self.rows, *self.failures):
            row.validate()
        identities = [
            (row.gene_id, row.target_id, row.context_id, row.dose_id)
            for row in (*self.rows, *self.failures)
        ]
        if len(set(identities)) != len(identities):
            raise OECapacityValidationError(
                "screen rows and failures must have unique request identities."
            )


@dataclass(frozen=True)
class OECapacityOutputs:
    output_dir: str
    rows_path: str
    manifest_path: str
    report_path: str

    def validate(self) -> None:
        for field_name, value in (
            ("output_dir", self.output_dir),
            ("rows_path", self.rows_path),
            ("manifest_path", self.manifest_path),
            ("report_path", self.report_path),
        ):
            _require_text(value, field_name)


def _require_text(value: str, field_name: str) -> None:
    if not isinstance(value, str) or not value.strip():
        raise OECapacityValidationError(f"{field_name} must be non-empty.")


def _require_positive_number(value: float, field_name: str) -> None:
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(value)
        or value <= 0
    ):
        raise OECapacityValidationError(f"{field_name} must be a positive finite number.")


__all__ = [
    "CapacityAnchor",
    "CapacityAnchorCatalog",
    "ConfidenceLevel",
    "ConstraintChangeKind",
    "CapacityConstraintChange",
    "EvidenceSourceType",
    "GeneCapacitySpec",
    "GeneCapacityCatalog",
    "GeneCapacityCoverage",
    "GeneCapacityParameterSet",
    "GeneCapacityValidationIssue",
    "GeneCapacityValidationResult",
    "GeneEnzymeReactionMapping",
    "GPRRole",
    "OEDoseMode",
    "OEDoseSpec",
    "OECapacityError",
    "OECapacityParameterConflictError",
    "OECapacityConstraintBundle",
    "OECapacityComparisonResult",
    "OECapacityScenarioResult",
    "OECapacityPlan",
    "OECapacityOutputs",
    "OECapacityScreenConfig",
    "OECapacityScreenRequest",
    "OECapacityScreenResult",
    "OECapacityScreenRow",
    "OECapacityValidationError",
    "OEExecutionMode",
    "OEExecutionStatus",
    "ParameterScenario",
    "ParameterEstimate",
    "ParameterPolicy",
    "ResourceCostMode",
    "SolverSnapshot",
    "derive_mapping_execution_status",
]
