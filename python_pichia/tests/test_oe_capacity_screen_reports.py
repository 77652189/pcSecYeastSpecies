from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from pcsec_pichia.core.pichia_enzymes import (
    CombinedEnzymeData,
    MetabolicEnzymeData,
    SecretoryEnzymeData,
)
from pcsec_pichia.oe_capacity import (
    AbsoluteCapacityAvailability,
    OECapacityComparisonResult,
    OECapacityScenarioResult,
    OECapacityScreenConfig,
    OECapacityScreenRequest,
    OEDoseMode,
    OEDoseSpec,
    OEExecutionMode,
    OEExecutionStatus,
    OECalibrationStatus,
    OEProductMode,
    OEProductState,
    ParameterScenario,
    SolverSnapshot,
    run_gene_level_oe_screen,
    verify_oe_capacity_run,
    write_oe_capacity_outputs,
)


def test_real_boundary_manifest_verifies_with_incomplete_scenarios(
    tmp_path: Path,
) -> None:
    result = run_gene_level_oe_screen(
        _prepared_model(),
        (_request("missing_gene", OEExecutionMode.COMPARISON),),
        _config(),
    )
    run_dir = tmp_path / "boundary-run"
    write_oe_capacity_outputs(
        result,
        run_dir,
        run_identity={
            "run_id": "boundary-run",
            "target_ids": ["hLF"],
            "context_ids": ["ctx-missing_gene"],
            "case_kind": "boundary",
        },
        capacity_asset={
            "path": "Enzymedata/oe_capacity_baseline_capacity.json",
            "version": "reviewed-v1",
            "sha256": "b" * 64,
            "reviewed": True,
        },
    )

    verification = verify_oe_capacity_run(
        run_dir,
        {"target_id": "hLF", "case_kind": "boundary"},
    )

    assert verification["passed"] is True


def test_small_batch_screen_preserves_completed_and_non_executable_rows(
    monkeypatch,
) -> None:
    prepared = _prepared_model()
    captured_modes: list[OEExecutionMode] = []

    def fake_comparison(prepared_model, plan, solver_options):
        assert prepared_model is prepared
        assert solver_options == {"time_limit_seconds": "10"}
        captured_modes.append(plan.execution_mode)
        return _comparison(plan)

    monkeypatch.setattr(
        "pcsec_pichia.oe_capacity.simulation.run_gene_level_oe_comparison",
        fake_comparison,
    )
    result = run_gene_level_oe_screen(
        prepared,
        (
            _request("G1", OEExecutionMode.COMPARISON),
            _request("missing_gene", OEExecutionMode.COMPARISON),
        ),
        _config(),
    )
    result.validate()
    assert captured_modes == [OEExecutionMode.RELATIVE_GENE_CAPACITY]
    assert len(result.rows) == 1
    assert len(result.failures) == 1
    row = result.rows[0]
    assert row.gene_id == "G1"
    assert row.baseline_objective == 1.0
    assert row.proxy_objective == 1.1
    assert row.gene_capacity_objective is None
    assert row.gene_capacity_vs_baseline_delta is None
    assert row.gene_capacity_vs_proxy_delta is None
    assert "reviewed_baseline_capacity" in row.missing_information
    assert result.failures[0].execution_mode is OEExecutionMode.NOT_EXECUTABLE


def test_feature_off_screen_executes_only_legacy_proxy(monkeypatch) -> None:
    prepared = _prepared_model()
    captured_modes: list[OEExecutionMode] = []

    def fake_comparison(_prepared_model, plan, _solver_options):
        captured_modes.append(plan.execution_mode)
        return _comparison(plan)

    monkeypatch.setattr(
        "pcsec_pichia.oe_capacity.simulation.run_gene_level_oe_comparison",
        fake_comparison,
    )
    result = run_gene_level_oe_screen(
        prepared,
        (_request("G1", OEExecutionMode.COMPARISON),),
        OECapacityScreenConfig(
            feature_enabled=False,
            compare_proxy=True,
            parameter_scenarios=(ParameterScenario.NOMINAL,),
            growth_rate=0.1,
        ),
    )

    assert captured_modes == [OEExecutionMode.REACTION_PROXY]
    assert result.rows[0].execution_mode is OEExecutionMode.REACTION_PROXY
    assert result.rows[0].proxy_objective == 1.1
    assert result.rows[0].gene_capacity_objective is None
    assert "only legacy reaction proxy" in result.warnings[0]


def test_screen_outputs_write_jsonl_manifest_and_boundary_report(
    tmp_path: Path,
    monkeypatch,
) -> None:
    monkeypatch.setattr(
        "pcsec_pichia.oe_capacity.simulation.run_gene_level_oe_comparison",
        lambda _prepared, plan, _options: _comparison(
            plan, include_failed_proxy=True
        ),
    )
    result = run_gene_level_oe_screen(
        _prepared_model(),
        (_request("G1", OEExecutionMode.COMPARISON),),
        _config(),
    )
    failed_row = result.failures[0]
    outputs = write_oe_capacity_outputs(result, tmp_path / "oe_capacity")

    outputs.validate()
    rows = [
        json.loads(line)
        for line in Path(outputs.rows_path).read_text(encoding="utf-8").splitlines()
    ]
    manifest = json.loads(Path(outputs.manifest_path).read_text(encoding="utf-8"))
    report = Path(outputs.report_path).read_text(encoding="utf-8")
    assert rows[0]["screen_status"] == "partial_failure"
    assert rows[0]["execution_mode"] == "relative_gene_capacity"
    assert rows[0]["proxy_objective"] == 1.1
    assert rows[0]["gene_capacity_objective"] is None
    assert rows[0]["relative_objective"] == 1.15
    assert rows[0]["gene_capacity_vs_baseline_delta"] is None
    assert rows[0]["proxy_attempts"][1]["success"] is False
    assert rows[0]["proxy_attempts"][1]["message"] == "infeasible proxy"
    assert manifest["schema_version"] == 2
    assert manifest["completed_count"] == 0
    assert manifest["failure_count"] == 1
    assert manifest["predicts_absolute_yield"] is False
    assert len(manifest["files"]["rows"]["sha256"]) == 64
    assert len(manifest["files"]["report"]["sha256"]) == 64
    assert manifest["run_identity"]["target_ids"] == ["hLF"]
    assert manifest["run_identity"]["context_ids"] == ["ctx-G1"]
    assert manifest["run_identity"]["gene_ids"] == ["G1"]
    assert manifest["status"]["state"] == "failed"
    assert manifest["scenario_completeness"]["required"] == ["nominal"]
    assert manifest["capacity_asset"]["reviewed"] is False
    assert "does not predict mg/L" in report
    assert "Legacy proxy comparison enabled: True" in report
    assert "Mapping and parameter traceability" in report
    assert "reviewed_baseline_capacity" in report
    assert "Reaction proxy attempts" in report
    assert "Relative scenario solver evidence" in report
    assert "infeasible proxy" in report
    assert "Scenario solver evidence" not in report


def test_absolute_partial_failure_preserves_scenario_message_in_all_outputs(
    tmp_path: Path,
    monkeypatch,
) -> None:
    monkeypatch.setattr(
        "pcsec_pichia.oe_capacity.simulation.run_gene_level_oe_comparison",
        lambda _prepared, plan, _options: _comparison(
            plan, include_failed_proxy=True
        ),
    )
    relative_result = run_gene_level_oe_screen(
        _prepared_model(),
        (_request("G1", OEExecutionMode.COMPARISON),),
        _config(),
    )
    source_row = relative_result.failures[0]
    baseline = replace(
        source_row.proxy_attempts[0],
        execution_mode=OEExecutionMode.NOT_EXECUTABLE,
        parameter_scenario=ParameterScenario.NOMINAL,
        attempt_id="nominal-absolute-baseline",
    )
    perturbed = replace(
        source_row.proxy_attempts[1],
        execution_mode=OEExecutionMode.GENE_CAPACITY,
        parameter_scenario=ParameterScenario.NOMINAL,
        attempt_id="nominal-absolute-perturbed",
        message="high scenario infeasible",
    )
    absolute_row = replace(
        source_row,
        execution_mode=OEExecutionMode.COMPARISON,
        execution_status=OEExecutionStatus.GENE_LEVEL_EXECUTABLE,
        product_mode=OEProductMode.ABSOLUTE_CAPACITY,
        product_state=OEProductState.ABSOLUTE_AVAILABLE,
        absolute_capacity_availability=AbsoluteCapacityAvailability.AVAILABLE_REVIEWED,
        calibration_status=OECalibrationStatus.REVIEWED_ABSOLUTE,
        absolute_solver_allowed=True,
        scenario_results=(
            OECapacityScenarioResult(
                parameter_scenario=ParameterScenario.NOMINAL,
                baseline=baseline,
                perturbed=perturbed,
                failure_reason="scenario_perturbation_failed",
            ),
        ),
        relative_scenario_results=(),
        relative_capacity_factors=(),
        relative_objective=None,
        relative_vs_baseline_delta=None,
        relative_vs_proxy_delta=None,
        nominal_capacity=1.0,
        nominal_capacities=(("R1_complex_formation", 1.0),),
        screen_status="partial_failure",
    )
    result = replace(relative_result, rows=(), failures=(absolute_row,))
    outputs = write_oe_capacity_outputs(result, tmp_path / "absolute-failure")
    row = json.loads(Path(outputs.rows_path).read_text(encoding="utf-8"))
    manifest = json.loads(Path(outputs.manifest_path).read_text(encoding="utf-8"))
    report = Path(outputs.report_path).read_text(encoding="utf-8")

    assert row["scenario_results"][0]["failure_reason"] == "scenario_perturbation_failed"
    assert row["scenario_results"][0]["perturbed"]["message"] == "high scenario infeasible"
    assert manifest["status"]["state"] == "failed"
    assert "Scenario solver evidence" in report
    assert "scenario_perturbation_failed" in report
    assert "high scenario infeasible" in report


def _request(gene_id: str, mode: OEExecutionMode) -> OECapacityScreenRequest:
    return OECapacityScreenRequest(
        gene_id=gene_id,
        target_id="hLF",
        context_id=f"ctx-{gene_id}",
        dose=OEDoseSpec(
            dose_id="2x",
            dose_mode=OEDoseMode.EXPLICIT_MULTIPLIER,
            expression_multiplier=2.0,
        ),
        execution_mode=mode,
    )


def _config() -> OECapacityScreenConfig:
    return OECapacityScreenConfig(
        feature_enabled=True,
        compare_proxy=True,
        parameter_scenarios=(ParameterScenario.NOMINAL,),
        growth_rate=0.1,
        solver_options=(("time_limit_seconds", "10"),),
    )


def _prepared_model() -> SimpleNamespace:
    model = SimpleNamespace(
        source_file="screen-fixture",
        rxns=["R1", "R1_complex_formation", "BIOMASS"],
        rules=["x(1)", "", ""],
        gr_rules=["G1", "", ""],
        genes=["G1"],
        gene_index={"G1": 0},
        reaction_index={"R1": 0, "R1_complex_formation": 1, "BIOMASS": 2},
        lb=np.array([0.0, 0.0, 0.1]),
        ub=np.array([1000.0, 1.0, 0.1]),
    )
    metabolic = MetabolicEnzymeData(
        source_file=Path("metabolic.mat"),
        enzymes=["R1_complex"],
        kcat=np.array([100.0]),
    )
    combined = CombinedEnzymeData(
        source_files=(Path("combined.mat"),),
        enzymes=["R1_complex"],
        kcat=np.array([100.0]),
        enzyme_mw=np.array([60000.0]),
        proteins=[],
        protein_length=np.array([]),
        protein_mw=np.array([]),
    )
    secretory = SecretoryEnzymeData(
        source_file=Path("secretory.mat"),
        reaction_coefficient_sources=(),
        complexes=[],
        compartments=[],
        kcat=np.array([]),
        coefficient_refs=[],
        reaction_coefficients={},
    )
    return SimpleNamespace(
        target_id="hLF",
        fixed_model=model,
        exchange_reaction_id="EX_TARGET",
        metabolic=metabolic,
        secretory=secretory,
        combined=combined,
    )


def _comparison(
    plan, *, include_failed_proxy: bool = False
) -> OECapacityComparisonResult:
    baseline = SolverSnapshot(
        execution_mode=OEExecutionMode.NOT_EXECUTABLE,
        backend="scipy_highs_reference",
        solver_status="optimal",
        success=True,
        secretion_objective=1.0,
        growth_retention=1.0,
        max_feasible_growth_rate=None,
        protein_resource_cost=0.1,
    )
    proxy = None
    if plan.proxy_reaction_ids or plan.execution_mode in {
        OEExecutionMode.REACTION_PROXY,
        OEExecutionMode.COMPARISON,
    }:
        proxy = SolverSnapshot(
            execution_mode=OEExecutionMode.REACTION_PROXY,
            backend="scipy_highs_reference",
            solver_status="optimal",
            success=True,
            secretion_objective=1.1,
            growth_retention=None,
            max_feasible_growth_rate=None,
            protein_resource_cost=None,
            attempt_id="R1",
        )
    proxy_attempts = ()
    if proxy is not None:
        proxy_attempts = (proxy,)
        if include_failed_proxy:
            proxy_attempts = proxy_attempts + (
                SolverSnapshot(
                execution_mode=OEExecutionMode.REACTION_PROXY,
                backend="scipy_highs_reference",
                solver_status="infeasible",
                success=False,
                secretion_objective=None,
                growth_retention=None,
                max_feasible_growth_rate=None,
                protein_resource_cost=None,
                message="infeasible proxy",
                attempt_id="R2",
                ),
            )
    scenarios = ()
    if plan.execution_mode in {
        OEExecutionMode.GENE_CAPACITY,
        OEExecutionMode.COMPARISON,
    }:
        scenarios = (
            SolverSnapshot(
                execution_mode=OEExecutionMode.GENE_CAPACITY,
                backend="scipy_highs_reference",
                solver_status="optimal",
                success=True,
                secretion_objective=1.2,
                growth_retention=1.0,
                max_feasible_growth_rate=None,
                protein_resource_cost=0.2,
                parameter_scenario=ParameterScenario.NOMINAL,
            ),
        )
    relative_scenarios = ()
    relative_results = ()
    if plan.execution_mode is OEExecutionMode.RELATIVE_GENE_CAPACITY:
        relative_baseline = replace(
            baseline,
            parameter_scenario=ParameterScenario.NOMINAL,
            attempt_id="nominal-relative-baseline",
        )
        relative_perturbed = SolverSnapshot(
            execution_mode=OEExecutionMode.RELATIVE_GENE_CAPACITY,
            backend="scipy_highs_reference",
            solver_status="optimal",
            success=True,
            secretion_objective=1.15,
            growth_retention=1.0,
            max_feasible_growth_rate=None,
            protein_resource_cost=0.15,
            parameter_scenario=ParameterScenario.NOMINAL,
            attempt_id="nominal-relative-2x",
        )
        relative_scenarios = (relative_perturbed,)
        relative_results = (
            OECapacityScenarioResult(
                parameter_scenario=ParameterScenario.NOMINAL,
                baseline=relative_baseline,
                perturbed=relative_perturbed,
                objective_delta=0.15,
            ),
        )
    return OECapacityComparisonResult(
        gene_id=plan.gene_id,
        target_id=plan.target_id,
        context_id=plan.context_id,
        execution_status=plan.execution_status,
        baseline=baseline,
        proxy=proxy,
        gene_capacity_scenarios=scenarios,
        gene_capacity_vs_baseline_delta=0.2 if scenarios else None,
        gene_capacity_vs_proxy_delta=0.1 if scenarios and proxy else None,
        protein_resource_cost_delta=0.1 if scenarios else None,
        missing_information=plan.missing_information,
        warnings=plan.warnings,
        proxy_attempts=proxy_attempts,
        relative_scenarios=relative_scenarios,
        relative_scenario_results=relative_results,
        relative_vs_baseline_delta=0.15 if relative_results else None,
        relative_vs_proxy_delta=0.05 if relative_results and proxy else None,
    )
