from __future__ import annotations

import pytest
from dataclasses import FrozenInstanceError

from pcsec_pichia.oe_capacity import (
    ConfidenceLevel,
    ConstraintChangeKind,
    CapacityConstraintChange,
    EvidenceSourceType,
    GeneCapacityCatalog,
    GeneCapacitySpec,
    GeneEnzymeReactionMapping,
    GPRRole,
    OEDoseMode,
    OEDoseSpec,
    OECapacityPlan,
    OECapacityConstraintBundle,
    OECapacityComparisonResult,
    OECapacityPhaseError,
    OECapacityScreenConfig,
    OECapacityScreenRequest,
    OECapacityValidationError,
    OEExecutionMode,
    OEExecutionStatus,
    SolverSnapshot,
    ParameterScenario,
    ParameterEstimate,
    ResourceCostMode,
    build_gene_enzyme_reaction_catalog,
    build_oe_dose_spec,
    build_oe_capacity_constraints,
)


def test_parameter_estimate_preserves_interval_unit_and_provenance() -> None:
    estimate = ParameterEstimate(
        parameter_name="kcat",
        nominal_value=120.0,
        lower_bound=80.0,
        upper_bound=180.0,
        unit="1/s",
        source_type=EvidenceSourceType.LOCAL_ENZYME_DATA,
        source_ref="Enzymedata/current",
        source_version="model-v1",
        confidence=ConfidenceLevel.HIGH,
    )

    estimate.validate()
    assert estimate.nominal_value == 120.0
    assert estimate.unit == "1/s"
    assert estimate.source_type is EvidenceSourceType.LOCAL_ENZYME_DATA

    with pytest.raises(OECapacityValidationError, match="lower_bound <= nominal_value"):
        ParameterEstimate(
            parameter_name="kcat",
            nominal_value=50.0,
            lower_bound=80.0,
            upper_bound=180.0,
            unit="1/s",
            source_type=EvidenceSourceType.LOCAL_ENZYME_DATA,
            source_ref="Enzymedata/current",
            source_version="model-v1",
            confidence=ConfidenceLevel.HIGH,
        ).validate()


def test_capacity_spec_and_plan_are_frozen_and_reject_proxy_disguised_as_gene_capacity() -> None:
    mapping = GeneEnzymeReactionMapping(
        mapping_id="map-G1-R1",
        model_fingerprint="model-v1",
        gene_id="G1",
        enzyme_id="R1_complex",
        reaction_id="R1",
        gpr_rule="x(1)",
        gpr_role=GPRRole.SINGLE_GENE,
        enzyme_variable_id="R1_complex",
        formation_or_dilution_reaction_id="R1_complex_formation",
        mapping_source=EvidenceSourceType.CURRENT_MODEL,
        mapping_confidence=ConfidenceLevel.HIGH,
        execution_status=OEExecutionStatus.GENE_LEVEL_EXECUTABLE,
    )
    dose = OEDoseSpec(
        dose_id="2x",
        dose_mode=OEDoseMode.EXPLICIT_MULTIPLIER,
        expression_multiplier=2.0,
    )
    kcat = ParameterEstimate(
        "kcat",
        120.0,
        80.0,
        180.0,
        "1/s",
        EvidenceSourceType.LOCAL_ENZYME_DATA,
        "Enzymedata/current",
        "model-v1",
        ConfidenceLevel.HIGH,
    )
    molecular_weight = ParameterEstimate(
        "molecular_weight",
        60.0,
        60.0,
        60.0,
        "kDa",
        EvidenceSourceType.LOCAL_ENZYME_DATA,
        "Enzymedata/current",
        "model-v1",
        ConfidenceLevel.HIGH,
    )
    baseline = ParameterEstimate(
        "baseline_enzyme_amount",
        1.0,
        0.5,
        1.5,
        "relative_amount",
        EvidenceSourceType.CURRENT_MODEL,
        "model formation variable",
        "model-v1",
        ConfidenceLevel.MEDIUM,
    )
    spec = GeneCapacitySpec(
        mapping=mapping,
        kcat=kcat,
        molecular_weight=molecular_weight,
        baseline_enzyme_amount=baseline,
        complex_stoichiometry=None,
        dose=dose,
        parameter_scenario=ParameterScenario.NOMINAL,
        resource_cost_mode=ResourceCostMode.CURRENT_PROTEIN_POOL,
    )
    spec.validate()

    plan = OECapacityPlan(
        gene_id="G1",
        target_id="hLF",
        context_id="ctx-hlf",
        requested_dose=dose,
        execution_mode=OEExecutionMode.GENE_CAPACITY,
        execution_status=OEExecutionStatus.GENE_LEVEL_EXECUTABLE,
        executable_capacity_specs=(spec,),
    )
    plan.validate()
    with pytest.raises(FrozenInstanceError):
        plan.gene_id = "G2"  # type: ignore[misc]

    with pytest.raises(OECapacityValidationError, match="reaction_proxy"):
        OECapacityPlan(
            gene_id="G1",
            target_id="hLF",
            context_id="ctx-hlf",
            requested_dose=dose,
            execution_mode=OEExecutionMode.GENE_CAPACITY,
            execution_status=OEExecutionStatus.PROXY_ONLY,
            proxy_reaction_ids=("R1",),
        ).validate()

    with pytest.raises(OECapacityValidationError, match="both executable_capacity_specs"):
        OECapacityPlan(
            gene_id="G1",
            target_id="hLF",
            context_id="ctx-hlf",
            requested_dose=dose,
            execution_mode=OEExecutionMode.COMPARISON,
            execution_status=OEExecutionStatus.GENE_LEVEL_EXECUTABLE,
            executable_capacity_specs=(spec,),
        ).validate()


def test_catalog_requires_one_model_fingerprint_and_unique_mapping_ids() -> None:
    mapping = GeneEnzymeReactionMapping(
        mapping_id="map-G1-R1",
        model_fingerprint="model-v1",
        gene_id="G1",
        enzyme_id="R1_complex",
        reaction_id="R1",
        gpr_rule="x(1)",
        gpr_role=GPRRole.SINGLE_GENE,
        enzyme_variable_id="R1_complex",
        formation_or_dilution_reaction_id="R1_complex_formation",
        mapping_source=EvidenceSourceType.CURRENT_MODEL,
        mapping_confidence=ConfidenceLevel.HIGH,
        execution_status=OEExecutionStatus.GENE_LEVEL_EXECUTABLE,
    )
    catalog = GeneCapacityCatalog(
        model_fingerprint="model-v1",
        mappings=(mapping,),
        source_refs=("Model/current", "Enzymedata/current"),
    )
    catalog.validate()

    with pytest.raises(OECapacityValidationError, match="duplicate mapping_id"):
        GeneCapacityCatalog(
            model_fingerprint="model-v1",
            mappings=(mapping, mapping),
        ).validate()


def test_constraint_bundle_keeps_gene_capacity_changes_separate_from_proxy_bounds() -> None:
    plan = _gene_capacity_plan()
    change = CapacityConstraintChange(
        change_id="formation-R1-2x",
        scenario=ParameterScenario.NOMINAL,
        change_kind=ConstraintChangeKind.FORMATION_DILUTION_BOUND,
        constraint_block="metabolic_coupling",
        variable_id="R1_complex_formation",
        reaction_id="R1",
        old_value=1.0,
        new_value=2.0,
        unit="relative_capacity",
        source_ref="current model formation variable",
        resource_cost_mode=ResourceCostMode.CURRENT_PROTEIN_POOL,
    )
    change.validate()
    bundle = OECapacityConstraintBundle(
        model_fingerprint="model-v1",
        plan=plan,
        changes=(change,),
    )
    bundle.validate()

    with pytest.raises(OECapacityValidationError, match="reaction_bound_proxy"):
        OECapacityConstraintBundle(
            model_fingerprint="model-v1",
            plan=plan,
            changes=(
                CapacityConstraintChange(
                    change_id="proxy-R1",
                    scenario=ParameterScenario.NOMINAL,
                    change_kind=ConstraintChangeKind.REACTION_BOUND_PROXY,
                    constraint_block="reaction_bounds",
                    variable_id="R1",
                    reaction_id="R1",
                    old_value=1000.0,
                    new_value=2000.0,
                    unit="mmol/gDCW/h",
                    source_ref="legacy proxy",
                    resource_cost_mode=ResourceCostMode.NOT_AVAILABLE,
                ),
            ),
        ).validate()


def test_phase_error_remains_explicit_for_future_contracts() -> None:
    error = OECapacityPhaseError("future_api", 5)

    assert error.api_name == "future_api"
    assert error.required_round == 5
    assert "Round 5" in str(error)


def test_execution_statuses_and_comparison_snapshots_are_explicit() -> None:
    assert {status.value for status in OEExecutionStatus} == {
        "gene_level_executable",
        "partial_mapping",
        "isoenzyme_ambiguous",
        "complex_limited",
        "external_evidence_only",
        "categorical_dose_only",
        "proxy_only",
        "unresolved",
    }
    baseline = SolverSnapshot(
        execution_mode=OEExecutionMode.NOT_EXECUTABLE,
        backend="scipy_highs_reference",
        solver_status="optimal",
        success=True,
        secretion_objective=1.0,
        growth_retention=1.0,
        max_feasible_growth_rate=0.2,
        protein_resource_cost=0.37,
        constraint_counts=(("equalities", 10), ("inequalities", 4)),
    )
    proxy = SolverSnapshot(
        execution_mode=OEExecutionMode.REACTION_PROXY,
        backend="scipy_highs_reference",
        solver_status="optimal",
        success=True,
        secretion_objective=1.2,
        growth_retention=0.95,
        max_feasible_growth_rate=0.19,
        protein_resource_cost=0.37,
    )
    gene_capacity = SolverSnapshot(
        execution_mode=OEExecutionMode.GENE_CAPACITY,
        backend="scipy_highs_reference",
        solver_status="optimal",
        success=True,
        secretion_objective=1.1,
        growth_retention=0.93,
        max_feasible_growth_rate=0.18,
        protein_resource_cost=0.39,
        parameter_scenario=ParameterScenario.NOMINAL,
    )
    result = OECapacityComparisonResult(
        gene_id="G1",
        target_id="hLF",
        context_id="ctx-hlf",
        execution_status=OEExecutionStatus.GENE_LEVEL_EXECUTABLE,
        baseline=baseline,
        proxy=proxy,
        gene_capacity_scenarios=(gene_capacity,),
        gene_capacity_vs_baseline_delta=0.1,
        gene_capacity_vs_proxy_delta=-0.1,
        protein_resource_cost_delta=0.02,
    )
    result.validate()
    assert result.proxy is not result.gene_capacity_scenarios[0]


def test_screen_contract_rejects_categorical_dose_for_gene_capacity_execution() -> None:
    categorical = OEDoseSpec(
        dose_id="strong-promoter",
        dose_mode=OEDoseMode.CATEGORICAL_ONLY,
        promoter="pGAP",
    )
    with pytest.raises(OECapacityValidationError, match="categorical_only"):
        OECapacityScreenRequest(
            gene_id="G1",
            target_id="hLF",
            context_id="ctx-hlf",
            dose=categorical,
            execution_mode=OEExecutionMode.GENE_CAPACITY,
        ).validate()

    config = OECapacityScreenConfig(
        feature_enabled=True,
        compare_proxy=True,
        parameter_scenarios=(
            ParameterScenario.LOW,
            ParameterScenario.NOMINAL,
            ParameterScenario.HIGH,
        ),
        growth_rate=0.1,
    )
    config.validate()


def test_dose_and_mapping_contracts_keep_proxy_and_gene_capacity_distinct() -> None:
    categorical = OEDoseSpec(
        dose_id="strong-promoter",
        dose_mode=OEDoseMode.CATEGORICAL_ONLY,
        promoter="pGAP",
    )
    categorical.validate()
    assert categorical.expression_multiplier is None

    with pytest.raises(OECapacityValidationError, match="categorical_only"):
        OEDoseSpec(
            dose_id="invalid-category",
            dose_mode=OEDoseMode.CATEGORICAL_ONLY,
            promoter="pGAP",
            expression_multiplier=2.0,
        ).validate()


def _gene_capacity_plan() -> OECapacityPlan:
    mapping = GeneEnzymeReactionMapping(
        mapping_id="map-G1-R1",
        model_fingerprint="model-v1",
        gene_id="G1",
        enzyme_id="R1_complex",
        reaction_id="R1",
        gpr_rule="x(1)",
        gpr_role=GPRRole.SINGLE_GENE,
        enzyme_variable_id="R1_complex",
        formation_or_dilution_reaction_id="R1_complex_formation",
        mapping_source=EvidenceSourceType.CURRENT_MODEL,
        mapping_confidence=ConfidenceLevel.HIGH,
        execution_status=OEExecutionStatus.GENE_LEVEL_EXECUTABLE,
    )
    dose = OEDoseSpec(
        dose_id="2x",
        dose_mode=OEDoseMode.EXPLICIT_MULTIPLIER,
        expression_multiplier=2.0,
    )
    parameters = {
        "kcat": ParameterEstimate(
            "kcat",
            120.0,
            80.0,
            180.0,
            "1/s",
            EvidenceSourceType.LOCAL_ENZYME_DATA,
            "Enzymedata/current",
            "model-v1",
            ConfidenceLevel.HIGH,
        ),
        "molecular_weight": ParameterEstimate(
            "molecular_weight",
            60.0,
            60.0,
            60.0,
            "kDa",
            EvidenceSourceType.LOCAL_ENZYME_DATA,
            "Enzymedata/current",
            "model-v1",
            ConfidenceLevel.HIGH,
        ),
        "baseline": ParameterEstimate(
            "baseline_enzyme_amount",
            1.0,
            0.5,
            1.5,
            "relative_amount",
            EvidenceSourceType.CURRENT_MODEL,
            "model formation variable",
            "model-v1",
            ConfidenceLevel.MEDIUM,
        ),
    }
    spec = GeneCapacitySpec(
        mapping=mapping,
        kcat=parameters["kcat"],
        molecular_weight=parameters["molecular_weight"],
        baseline_enzyme_amount=parameters["baseline"],
        complex_stoichiometry=None,
        dose=dose,
        parameter_scenario=ParameterScenario.NOMINAL,
        resource_cost_mode=ResourceCostMode.CURRENT_PROTEIN_POOL,
    )
    return OECapacityPlan(
        gene_id="G1",
        target_id="hLF",
        context_id="ctx-hlf",
        requested_dose=dose,
        execution_mode=OEExecutionMode.GENE_CAPACITY,
        execution_status=OEExecutionStatus.GENE_LEVEL_EXECUTABLE,
        executable_capacity_specs=(spec,),
    )

    executable = GeneEnzymeReactionMapping(
        mapping_id="map-G1-R1",
        model_fingerprint="model-v1",
        gene_id="G1",
        enzyme_id="R1_complex",
        reaction_id="R1",
        gpr_rule="x(1)",
        gpr_role=GPRRole.SINGLE_GENE,
        enzyme_variable_id="R1_complex",
        formation_or_dilution_reaction_id="R1_complex_formation",
        mapping_source=EvidenceSourceType.CURRENT_MODEL,
        mapping_confidence=ConfidenceLevel.HIGH,
        execution_status=OEExecutionStatus.GENE_LEVEL_EXECUTABLE,
    )
    executable.validate()

    with pytest.raises(OECapacityValidationError, match="external_evidence_only"):
        GeneEnzymeReactionMapping(
            mapping_id="external-G2",
            model_fingerprint="model-v1",
            gene_id="G2",
            enzyme_id="external-enzyme",
            reaction_id="external-reaction",
            gpr_rule="G2",
            gpr_role=GPRRole.SINGLE_GENE,
            mapping_source=EvidenceSourceType.EXTERNAL_PICHIA_MODEL,
            mapping_confidence=ConfidenceLevel.LOW,
            execution_status=OEExecutionStatus.GENE_LEVEL_EXECUTABLE,
        ).validate()

    with pytest.raises(OECapacityValidationError, match="assembly_evidence_ref"):
        GeneEnzymeReactionMapping(
            mapping_id="complex-G3",
            model_fingerprint="model-v1",
            gene_id="G3",
            enzyme_id="C1",
            reaction_id="R3",
            gpr_rule="x(3) & x(4)",
            gpr_role=GPRRole.COMPLEX_SUBUNIT,
            complex_id="C1",
            subunit_ids=("G3", "G4"),
            subunit_stoichiometry=(("G3", 1.0), ("G4", 1.0)),
            enzyme_variable_id="C1",
            formation_or_dilution_reaction_id="C1_formation",
            mapping_source=EvidenceSourceType.CURRENT_MODEL,
            mapping_confidence=ConfidenceLevel.HIGH,
            execution_status=OEExecutionStatus.GENE_LEVEL_EXECUTABLE,
        ).validate()
