from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from pcsec_pichia.oe_capacity import (
    AbsoluteCapacityAvailability,
    ConfidenceLevel,
    EvidenceSourceType,
    GeneCapacityCatalog,
    GeneCapacityParameterSet,
    GeneEnzymeReactionMapping,
    GPRRole,
    OECapacityScreenConfig,
    OECapacityScreenRequest,
    OECapacityValidationError,
    OECalibrationStatus,
    OEExecutionMode,
    OEExecutionStatus,
    OEProductMode,
    OEProductState,
    ParameterEstimate,
    ParameterPolicy,
    ParameterScenario,
    build_oe_dose_spec,
    fingerprint_oe_capacity_model,
    plan_gene_level_overexpression,
    resolve_oe_product_plan,
    run_gene_level_oe_screen,
    summarize_oe_product_candidate,
    write_oe_capacity_outputs,
)


@pytest.mark.parametrize(
    ("requested", "feature_enabled", "expected_state", "expected_execution"),
    (
        (
            OEProductMode.REACTION_PROXY,
            True,
            OEProductState.REACTION_PROXY,
            OEExecutionMode.REACTION_PROXY,
        ),
        (
            OEProductMode.RELATIVE_UNCALIBRATED,
            True,
            OEProductState.RELATIVE_UNCALIBRATED,
            OEExecutionMode.RELATIVE_GENE_CAPACITY,
        ),
        (
            OEProductMode.ABSOLUTE_CAPACITY,
            True,
            OEProductState.ABSOLUTE_UNAVAILABLE,
            OEExecutionMode.NOT_EXECUTABLE,
        ),
        (
            OEProductMode.RELATIVE_UNCALIBRATED,
            False,
            OEProductState.REACTION_PROXY,
            OEExecutionMode.REACTION_PROXY,
        ),
        (
            OEProductMode.NOT_EXECUTABLE,
            True,
            OEProductState.NOT_EXECUTABLE,
            OEExecutionMode.NOT_EXECUTABLE,
        ),
    ),
)
def test_core_resolver_derives_mutually_distinct_product_states(
    requested: OEProductMode,
    feature_enabled: bool,
    expected_state: OEProductState,
    expected_execution: OEExecutionMode,
) -> None:
    plan, _prepared = _unanchored_plan()

    resolved = resolve_oe_product_plan(
        plan,
        requested_mode=requested,
        feature_enabled=feature_enabled,
        compare_proxy=True,
    )

    assert resolved.product_state is expected_state
    assert resolved.execution_mode is expected_execution
    assert not resolved.absolute_solver_allowed
    assert not resolved.executable_capacity_specs
    if expected_state is not OEProductState.RELATIVE_UNCALIBRATED:
        assert not resolved.relative_scenario_specs


def test_relative_product_keeps_mapping_dose_uncertainty_sources_and_limits() -> None:
    plan, _prepared = _unanchored_plan()
    resolved = resolve_oe_product_plan(
        plan,
        requested_mode=OEProductMode.RELATIVE_UNCALIBRATED,
        feature_enabled=True,
        compare_proxy=True,
    )

    summary = summarize_oe_product_candidate(resolved)
    assert summary["product_state"] == "relative_uncalibrated"
    assert resolved.execution_status is OEExecutionStatus.GENE_LEVEL_EXECUTABLE
    assert summary["absolute_capacity_availability"] == "unavailable_missing_reviewed_anchor"
    assert summary["absolute_solver_allowed"] is False
    assert summary["mapping_ids"] == ["map-G1-R1"]
    scenarios = summary["relative_scenarios"]
    assert isinstance(scenarios, list) and len(scenarios) == 3
    assert {item["scenario"] for item in scenarios} == {"low", "nominal", "high"}
    assert all(item["relative_capacity_factor"] == 2.0 for item in scenarios)
    assert all(item["parameter_sources"] for item in scenarios)
    assert all(item["limitations"] for item in scenarios)


def test_explicit_absolute_without_anchor_never_calls_comparison_solver(
    monkeypatch,
    tmp_path: Path,
) -> None:
    _plan, prepared = _unanchored_plan()
    calls = {"comparison": 0}

    def forbidden(*args, **kwargs):
        calls["comparison"] += 1
        raise AssertionError("absolute-unavailable must stop before solver dispatch")

    monkeypatch.setattr(
        "pcsec_pichia.oe_capacity.simulation.run_gene_level_oe_comparison",
        forbidden,
    )
    request = OECapacityScreenRequest(
        gene_id="G1",
        target_id="hLF",
        context_id="ctx-hlf",
        dose=build_oe_dose_spec(
            {
                "dose_id": "2x",
                "dose_mode": "explicit_multiplier",
                "expression_multiplier": 2.0,
            }
        ),
        execution_mode=OEExecutionMode.GENE_CAPACITY,
        product_mode=OEProductMode.ABSOLUTE_CAPACITY,
    )
    result = run_gene_level_oe_screen(
        prepared,
        (request,),
        OECapacityScreenConfig(
            feature_enabled=True,
            compare_proxy=True,
            parameter_scenarios=(
                ParameterScenario.LOW,
                ParameterScenario.NOMINAL,
                ParameterScenario.HIGH,
            ),
            growth_rate=0.1,
        ),
    )

    assert calls["comparison"] == 0
    assert not result.rows
    row = result.failures[0]
    assert row.product_state is OEProductState.ABSOLUTE_UNAVAILABLE
    assert row.absolute_capacity_availability is (
        AbsoluteCapacityAvailability.UNAVAILABLE_MISSING_REVIEWED_ANCHOR
    )
    assert row.calibration_status is OECalibrationStatus.UNAVAILABLE
    assert row.baseline_objective is None
    assert row.proxy_objective is None
    assert row.gene_capacity_objective is None
    assert row.nominal_capacity is None
    outputs = write_oe_capacity_outputs(result, tmp_path / "absolute-unavailable")
    payload = json.loads(Path(outputs.rows_path).read_text(encoding="utf-8"))
    report = Path(outputs.report_path).read_text(encoding="utf-8")
    manifest = json.loads(Path(outputs.manifest_path).read_text(encoding="utf-8"))
    assert payload["product_state"] == "absolute_unavailable"
    assert payload["absolute_solver_allowed"] is False
    assert manifest["coverage"]["by_product_state"] == {"absolute_unavailable": 1}
    assert "absolute_unavailable" in report
    assert "no solver was called" in report


def test_manual_unreviewed_baseline_cannot_restore_absolute_execution() -> None:
    plan, _prepared = _unanchored_plan()
    parameter_set = plan.relative_scenario_specs[0]
    fake_baseline = ParameterEstimate(
        parameter_name="baseline_enzyme_amount",
        nominal_value=1.0,
        lower_bound=1.0,
        upper_bound=1.0,
        unit="model_flux",
        source_type=EvidenceSourceType.LOCAL_ENZYME_DATA,
        source_ref="caller-supplied-1.0",
        source_version="unreviewed",
        confidence=ConfidenceLevel.HIGH,
    )
    with pytest.raises(Exception, match="reviewed capacity_anchor_binding"):
        from pcsec_pichia.oe_capacity import GeneCapacitySpec, ResourceCostMode

        GeneCapacitySpec(
            mapping=parameter_set.mapping,
            kcat=parameter_set.kcat,
            molecular_weight=parameter_set.molecular_weight,
            baseline_enzyme_amount=fake_baseline,
            complex_stoichiometry=None,
            dose=parameter_set.dose,
            parameter_scenario=ParameterScenario.NOMINAL,
            resource_cost_mode=ResourceCostMode.CURRENT_PROTEIN_POOL,
        ).validate()


@pytest.mark.parametrize(
    "mutation",
    (
        {"execution_mode": OEExecutionMode.REACTION_PROXY},
        {"calibration_status": OECalibrationStatus.PROXY_ONLY},
        {"product_mode": OEProductMode.REACTION_PROXY},
        {"absolute_solver_allowed": True},
    ),
)
def test_relative_state_matrix_rejects_contradictory_manual_plans(mutation) -> None:
    plan, _prepared = _unanchored_plan()
    relative = resolve_oe_product_plan(
        plan,
        requested_mode=OEProductMode.RELATIVE_UNCALIBRATED,
        feature_enabled=True,
        compare_proxy=True,
    )
    with pytest.raises(OECapacityValidationError):
        replace(relative, **mutation).validate()


def _unanchored_plan():
    model = SimpleNamespace(
        rxns=["BIOMASS", "R1", "R1_complex_formation"],
        rules=["", "x(1)", ""],
        gr_rules=["", "G1", ""],
        genes=["G1"],
        gene_index={"G1": 0},
        reaction_index={"BIOMASS": 0, "R1": 1, "R1_complex_formation": 2},
        lb=np.array([0.1, 0.0, 0.0]),
        ub=np.array([0.1, 1000.0, 1000.0]),
    )
    fingerprint = fingerprint_oe_capacity_model(model)
    mapping = GeneEnzymeReactionMapping(
        mapping_id="map-G1-R1",
        model_fingerprint=fingerprint,
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
    parameter_set = GeneCapacityParameterSet(
        parameter_set_id="current-map-G1-R1",
        mapping_id=mapping.mapping_id,
        gene_id=mapping.gene_id,
        enzyme_id=mapping.enzyme_id,
        kcat=_estimate("kcat", 100.0, "1/h"),
        molecular_weight=_estimate("molecular_weight", 60000.0, "g/mol"),
        baseline_enzyme_amount=None,
    )
    policy = ParameterPolicy(parameter_sets=(parameter_set,))
    catalog = GeneCapacityCatalog(fingerprint, (mapping,))
    dose = build_oe_dose_spec(
        {
            "dose_id": "2x",
            "dose_mode": "explicit_multiplier",
            "expression_multiplier": 2.0,
        }
    )
    plan = plan_gene_level_overexpression(
        model,
        "G1",
        "hLF",
        "ctx-hlf",
        dose,
        catalog,
        policy,
    )
    prepared = SimpleNamespace(
        fixed_model=model,
        target_id="hLF",
        exchange_reaction_id="EX_TARGET",
        metabolic=SimpleNamespace(),
        secretory=SimpleNamespace(),
        combined=SimpleNamespace(enzymes=("R1_complex",)),
        gene_capacity_catalog=catalog,
        parameter_policy=policy,
    )
    return plan, prepared


def _estimate(name: str, value: float, unit: str) -> ParameterEstimate:
    return ParameterEstimate(
        parameter_name=name,
        nominal_value=value,
        lower_bound=value * 0.8,
        upper_bound=value * 1.2,
        unit=unit,
        source_type=EvidenceSourceType.LOCAL_ENZYME_DATA,
        source_ref="current-model",
        source_version="v1",
        confidence=ConfidenceLevel.HIGH,
    )
