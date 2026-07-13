from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from scipy import sparse

from pcsec_pichia.core.pichia_enzymes import (
    CombinedEnzymeData,
    MetabolicEnzymeData,
    SecretoryEnzymeData,
)
from pcsec_pichia.oe_capacity import (
    ConfidenceLevel,
    ConstraintChangeKind,
    EvidenceSourceType,
    GeneCapacitySpec,
    GeneEnzymeReactionMapping,
    GPRRole,
    OEDoseMode,
    OEDoseSpec,
    OECapacityPlan,
    OECapacityValidationError,
    OEExecutionMode,
    OEExecutionStatus,
    ParameterEstimate,
    ParameterScenario,
    ResourceCostMode,
    build_oe_capacity_constraints,
    run_gene_level_oe_comparison,
)
from pcsec_pichia.probe import CobraModel


def test_constraint_builder_targets_formation_kcat_and_resource_terms_not_reaction_bound() -> None:
    plan = _plan(expression_multiplier=2.0)

    bundle = build_oe_capacity_constraints(_prepared_tiny_model(), plan)

    bundle.validate()
    kinds = {change.change_kind for change in bundle.changes}
    assert ConstraintChangeKind.FORMATION_DILUTION_BOUND in kinds
    assert ConstraintChangeKind.ENZYME_CAPACITY_COEFFICIENT in kinds
    assert ConstraintChangeKind.PROTEIN_RESOURCE_COEFFICIENT in kinds
    assert ConstraintChangeKind.REACTION_BOUND_PROXY not in kinds
    formation_change = next(
        change
        for change in bundle.changes
        if change.change_kind is ConstraintChangeKind.FORMATION_DILUTION_BOUND
    )
    assert formation_change.variable_id == "R1_complex_formation"
    assert formation_change.new_value == 2.0


def test_constraint_builder_rejects_mapping_not_present_in_prepared_model() -> None:
    prepared = _prepared_tiny_model()
    plan = _plan(expression_multiplier=2.0)
    spec = plan.executable_capacity_specs[0]
    missing_mapping = replace(
        spec.mapping,
        formation_or_dilution_reaction_id="missing_formation",
    )
    plan = replace(plan, executable_capacity_specs=(replace(spec, mapping=missing_mapping),))

    with pytest.raises(
        OECapacityValidationError,
        match="formation/dilution reaction is missing",
    ):
        build_oe_capacity_constraints(prepared, plan)


def test_constraint_builder_rejects_noncanonical_parameter_units() -> None:
    plan = _plan(expression_multiplier=2.0)
    spec = plan.executable_capacity_specs[0]
    plan = replace(
        plan,
        executable_capacity_specs=(
            replace(spec, kcat=replace(spec.kcat, unit="1/s")),
        ),
    )

    with pytest.raises(OECapacityValidationError, match="canonical unit 1/h"):
        build_oe_capacity_constraints(_prepared_tiny_model(), plan)


def test_single_candidate_solver_preserves_1x_and_distinguishes_gene_capacity_from_proxy() -> None:
    prepared = _prepared_tiny_model()

    one_x = run_gene_level_oe_comparison(
        prepared,
        _plan(expression_multiplier=1.0),
    )
    two_x = run_gene_level_oe_comparison(
        prepared,
        _plan(expression_multiplier=2.0),
    )

    assert one_x.baseline.success is True
    assert one_x.gene_capacity_scenarios[0].success is True
    assert one_x.gene_capacity_scenarios[0].secretion_objective == pytest.approx(
        one_x.baseline.secretion_objective,
        rel=1e-9,
        abs=1e-12,
    )
    assert one_x.gene_capacity_vs_baseline_delta == pytest.approx(0.0, abs=1e-12)

    assert two_x.proxy is not None and two_x.proxy.success is True
    assert two_x.proxy.secretion_objective == pytest.approx(
        two_x.baseline.secretion_objective,
        rel=1e-9,
    )
    assert two_x.gene_capacity_scenarios[0].secretion_objective is not None
    assert two_x.baseline.secretion_objective is not None
    assert (
        two_x.gene_capacity_scenarios[0].secretion_objective
        > two_x.baseline.secretion_objective
    )
    assert two_x.protein_resource_cost_delta is not None
    assert two_x.protein_resource_cost_delta > 0
    assert two_x.proxy.key_fluxes == ()
    assert "Selected legacy proxy reaction: R1." in two_x.proxy.warnings


def test_single_candidate_solver_applies_kcat_to_active_metabolic_coupling() -> None:
    prepared = _prepared_tiny_model()

    reference = run_gene_level_oe_comparison(
        prepared,
        _plan(expression_multiplier=2.0, kcat_nominal=100.0),
    )
    lower_kcat = run_gene_level_oe_comparison(
        prepared,
        _plan(expression_multiplier=2.0, kcat_nominal=50.0),
    )

    assert reference.gene_capacity_scenarios[0].secretion_objective is not None
    assert lower_kcat.gene_capacity_scenarios[0].secretion_objective is not None
    assert (
        lower_kcat.gene_capacity_scenarios[0].secretion_objective
        < reference.gene_capacity_scenarios[0].secretion_objective
    )


def test_zero_baseline_formation_flux_is_explicitly_reported() -> None:
    result = run_gene_level_oe_comparison(
        _prepared_tiny_model(use_direct_target_path=True),
        _plan(expression_multiplier=2.0),
    )

    assert result.baseline.success is True
    assert "nonzero_baseline_formation_flux" in result.missing_information
    assert any(
        "zero baseline flux" in warning
        for warning in result.warnings
    )
    assert dict(result.traceability)[
        "baseline_formation_flux:R1_complex_formation"
    ] == "0"


def test_proxy_only_plan_preserves_legacy_comparison_without_capacity_bundle() -> None:
    plan = _plan(expression_multiplier=2.0)
    plan = replace(
        plan,
        execution_mode=OEExecutionMode.REACTION_PROXY,
        execution_status=OEExecutionStatus.PROXY_ONLY,
        executable_capacity_specs=(),
        uncertainty_scenarios=(),
        missing_information=("capacity_parameters",),
    )

    result = run_gene_level_oe_comparison(_prepared_tiny_model(), plan)

    assert result.baseline.success is True
    assert result.proxy is not None and result.proxy.success is True
    assert result.gene_capacity_scenarios == ()
    assert result.skipped_reason == "gene_capacity_not_executable"
    assert dict(result.traceability)["constraint_change_count"] == "0"


def _plan(
    *,
    expression_multiplier: float,
    kcat_nominal: float = 100.0,
) -> OECapacityPlan:
    mapping = GeneEnzymeReactionMapping(
        mapping_id="map-G1-R1",
        model_fingerprint="model-v1",
        gene_id="G1",
        enzyme_id="R1_complex",
        reaction_id="R1",
        gpr_rule="G1",
        gpr_role=GPRRole.SINGLE_GENE,
        enzyme_variable_id="R1_complex_formation",
        formation_or_dilution_reaction_id="R1_complex_formation",
        mapping_source=EvidenceSourceType.CURRENT_MODEL,
        mapping_confidence=ConfidenceLevel.HIGH,
        execution_status=OEExecutionStatus.GENE_LEVEL_EXECUTABLE,
    )
    dose = OEDoseSpec(
        dose_id=f"{expression_multiplier:g}x",
        dose_mode=OEDoseMode.EXPLICIT_MULTIPLIER,
        expression_multiplier=expression_multiplier,
    )
    spec = GeneCapacitySpec(
        mapping=mapping,
        kcat=ParameterEstimate(
            "kcat",
            kcat_nominal,
            kcat_nominal,
            kcat_nominal,
            "1/h",
            EvidenceSourceType.LOCAL_ENZYME_DATA,
            "tiny",
            "v1",
            ConfidenceLevel.HIGH,
        ),
        molecular_weight=ParameterEstimate(
            "molecular_weight",
            60000.0,
            60000.0,
            60000.0,
            "g/mol",
            EvidenceSourceType.LOCAL_ENZYME_DATA,
            "tiny",
            "v1",
            ConfidenceLevel.HIGH,
        ),
        baseline_enzyme_amount=ParameterEstimate(
            "baseline_enzyme_amount",
            1.0,
            1.0,
            1.0,
            "relative_capacity",
            EvidenceSourceType.CURRENT_MODEL,
            "R1_complex_formation",
            "model-v1",
            ConfidenceLevel.HIGH,
        ),
        complex_stoichiometry=None,
        dose=dose,
        parameter_scenario=ParameterScenario.NOMINAL,
        resource_cost_mode=ResourceCostMode.CURRENT_PROTEIN_POOL,
    )
    return OECapacityPlan(
        gene_id="G1",
        target_id="hLF",
        context_id="tiny",
        requested_dose=dose,
        execution_mode=OEExecutionMode.COMPARISON,
        execution_status=OEExecutionStatus.GENE_LEVEL_EXECUTABLE,
        executable_capacity_specs=(spec,),
        proxy_reaction_ids=("R1",),
        uncertainty_scenarios=(ParameterScenario.NOMINAL,),
    )


def _prepared_tiny_model(*, use_direct_target_path: bool = False) -> SimpleNamespace:
    reactions = (
        "SUB_IN",
        "R1",
        "DIRECT_TARGET",
        "EX_TARGET",
        "R1_complex_formation",
        "R1_complex_dilution",
        "dilute_dummy",
        "BIOMASS",
    )
    matrix = np.zeros((3, len(reactions)))
    matrix[0, reactions.index("SUB_IN")] = 1.0
    matrix[0, reactions.index("R1")] = -1.0
    matrix[1, reactions.index("R1")] = 1.0
    matrix[0, reactions.index("DIRECT_TARGET")] = -0.5
    matrix[1, reactions.index("DIRECT_TARGET")] = 1.0
    matrix[1, reactions.index("EX_TARGET")] = -1.0
    matrix[2, reactions.index("R1_complex_formation")] = 1.0
    matrix[2, reactions.index("R1_complex_dilution")] = -1.0
    lower = np.zeros(len(reactions))
    upper = np.full(len(reactions), 1000.0)
    upper[reactions.index("SUB_IN")] = 10.0
    upper[reactions.index("DIRECT_TARGET")] = (
        1000.0 if use_direct_target_path else 0.0
    )
    upper[reactions.index("R1_complex_formation")] = 0.0001
    lower[reactions.index("BIOMASS")] = 0.1
    upper[reactions.index("BIOMASS")] = 0.1
    model = CobraModel(
        source_file="tiny-oe-capacity",
        rxns=list(reactions),
        mets=["A", "TARGET", "R1_complex"],
        genes=["G1"],
        lb=lower,
        ub=upper,
        b=np.zeros(3),
        s_matrix=sparse.csc_matrix(matrix),
        rules=["", "x(1)", "", "", "", "", "", ""],
        gr_rules=["", "G1", "", "", "", "", "", ""],
    )
    metabolic = MetabolicEnzymeData(
        source_file=Path("tiny-metabolic"),
        enzymes=["R1_complex"],
        kcat=np.array([100.0]),
    )
    secretory = SecretoryEnzymeData(
        source_file=Path("tiny-secretory"),
        reaction_coefficient_sources=(),
        complexes=[],
        compartments=[],
        kcat=np.array([]),
        coefficient_refs=[],
        reaction_coefficients={},
    )
    combined = CombinedEnzymeData(
        source_files=(),
        enzymes=["R1_complex"],
        kcat=np.array([100.0]),
        enzyme_mw=np.array([60000.0]),
        proteins=[],
        protein_length=np.array([]),
        protein_mw=np.array([]),
    )
    return SimpleNamespace(
        target_id="hLF",
        fixed_model=model,
        exchange_reaction_id="EX_TARGET",
        metabolic=metabolic,
        secretory=secretory,
        combined=combined,
    )
