from __future__ import annotations

import math
from dataclasses import dataclass
from enum import Enum


class OECapacityError(RuntimeError):
    """Base error for gene-level OE capacity workflows."""


class OECapacityValidationError(OECapacityError, ValueError):
    """Raised when an OE capacity contract violates a frozen invariant."""


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
    REACTION_PROXY = "reaction_proxy"
    COMPARISON = "comparison"
    NOT_EXECUTABLE = "not_executable"


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
        if self.is_transferred and self.source_type is not EvidenceSourceType.HOMOLOGY_TRANSFER:
            raise OECapacityValidationError(
                "is_transferred parameters must use source_type=homology_transfer."
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
        external_only_sources = {
            EvidenceSourceType.EXTERNAL_PICHIA_MODEL,
            EvidenceSourceType.PICHIA_LITERATURE,
            EvidenceSourceType.HOMOLOGY_TRANSFER,
            EvidenceSourceType.SMOKE_FIXTURE,
        }
        if (
            self.mapping_source in external_only_sources
            and self.execution_status is OEExecutionStatus.GENE_LEVEL_EXECUTABLE
        ):
            raise OECapacityValidationError(
                "external evidence must remain external_evidence_only until current-model "
                "gene, enzyme, and reaction mappings are confirmed."
            )
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
        if self.mapping.gpr_role is GPRRole.COMPLEX_SUBUNIT and self.complex_stoichiometry is None:
            raise OECapacityValidationError(
                "complex_subunit GeneCapacitySpec requires complex_stoichiometry."
            )
        if self.resource_cost_mode is ResourceCostMode.NOT_AVAILABLE:
            raise OECapacityValidationError(
                "executable GeneCapacitySpec requires an auditable resource cost mode."
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
            if self.execution_status is not OEExecutionStatus.PROXY_ONLY:
                raise OECapacityValidationError(
                    "reaction_proxy mode must use execution_status=proxy_only."
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


@dataclass(frozen=True)
class OECapacityScreenRequest:
    gene_id: str
    target_id: str
    context_id: str
    dose: OEDoseSpec
    execution_mode: OEExecutionMode

    def validate(self) -> None:
        _require_text(self.gene_id, "gene_id")
        _require_text(self.target_id, "target_id")
        _require_text(self.context_id, "context_id")
        self.dose.validate()
        if (
            self.execution_mode in {OEExecutionMode.GENE_CAPACITY, OEExecutionMode.COMPARISON}
            and self.dose.dose_mode is OEDoseMode.CATEGORICAL_ONLY
        ):
            raise OECapacityValidationError(
                "categorical_only dose cannot run gene_capacity or comparison execution."
            )


@dataclass(frozen=True)
class OECapacityScreenConfig:
    feature_enabled: bool
    compare_proxy: bool
    parameter_scenarios: tuple[ParameterScenario, ...]
    growth_rate: float
    solver_options: tuple[tuple[str, str], ...] = ()

    def validate(self) -> None:
        _require_positive_number(self.growth_rate, "growth_rate")
        if not self.parameter_scenarios:
            raise OECapacityValidationError(
                "parameter_scenarios must contain at least one scenario."
            )
        if len(set(self.parameter_scenarios)) != len(self.parameter_scenarios):
            raise OECapacityValidationError("parameter_scenarios must be unique.")


@dataclass(frozen=True)
class OECapacityScreenRow:
    gene_id: str
    target_id: str
    context_id: str
    execution_mode: OEExecutionMode
    execution_status: OEExecutionStatus
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
    gene_capacity_vs_proxy_delta: float | None
    protein_resource_cost_delta: float | None
    missing_information: tuple[str, ...] = ()
    warnings: tuple[str, ...] = ()

    def validate(self) -> None:
        for field_name, value in (
            ("gene_id", self.gene_id),
            ("target_id", self.target_id),
            ("context_id", self.context_id),
            ("dose_id", self.dose_id),
        ):
            _require_text(value, field_name)


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
    "ConfidenceLevel",
    "ConstraintChangeKind",
    "CapacityConstraintChange",
    "EvidenceSourceType",
    "GeneCapacitySpec",
    "GeneCapacityCatalog",
    "GeneCapacityValidationIssue",
    "GeneCapacityValidationResult",
    "GeneEnzymeReactionMapping",
    "GPRRole",
    "OEDoseMode",
    "OEDoseSpec",
    "OECapacityError",
    "OECapacityConstraintBundle",
    "OECapacityComparisonResult",
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
    "ResourceCostMode",
    "SolverSnapshot",
]
