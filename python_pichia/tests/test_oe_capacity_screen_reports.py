from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from pcsec_pichia.core.pichia_enzymes import (
    CombinedEnzymeData,
    MetabolicEnzymeData,
    SecretoryEnzymeData,
)
from pcsec_pichia.oe_capacity import (
    OECapacityComparisonResult,
    OECapacityScreenConfig,
    OECapacityScreenRequest,
    OEDoseMode,
    OEDoseSpec,
    OEExecutionMode,
    OEExecutionStatus,
    ParameterScenario,
    SolverSnapshot,
    run_gene_level_oe_screen,
    write_oe_capacity_outputs,
)


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
    assert captured_modes == [OEExecutionMode.COMPARISON]
    assert len(result.rows) == 1
    assert len(result.failures) == 1
    row = result.rows[0]
    assert row.gene_id == "G1"
    assert row.baseline_objective == 1.0
    assert row.proxy_objective == 1.1
    assert row.gene_capacity_objective == 1.2
    assert row.gene_capacity_vs_baseline_delta == 0.2
    assert row.gene_capacity_vs_proxy_delta == 0.1
    assert row.mapping_ids
    assert row.parameter_sources
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
        lambda _prepared, plan, _options: _comparison(plan),
    )
    result = run_gene_level_oe_screen(
        _prepared_model(),
        (_request("G1", OEExecutionMode.COMPARISON),),
        _config(),
    )

    outputs = write_oe_capacity_outputs(result, tmp_path / "oe_capacity")

    outputs.validate()
    rows = [
        json.loads(line)
        for line in Path(outputs.rows_path).read_text(encoding="utf-8").splitlines()
    ]
    manifest = json.loads(Path(outputs.manifest_path).read_text(encoding="utf-8"))
    report = Path(outputs.report_path).read_text(encoding="utf-8")
    assert rows[0]["screen_status"] == "completed"
    assert rows[0]["execution_mode"] == "comparison"
    assert rows[0]["proxy_objective"] == 1.1
    assert rows[0]["gene_capacity_objective"] == 1.2
    assert rows[0]["gene_capacity_vs_baseline_delta"] == 0.2
    assert manifest["completed_count"] == 1
    assert manifest["failure_count"] == 0
    assert manifest["predicts_absolute_yield"] is False
    assert len(manifest["rows_sha256"]) == 64
    assert "does not predict mg/L" in report
    assert "Legacy proxy comparison enabled: True" in report
    assert "Mapping and parameter traceability" in report
    assert "local_enzyme_data" in report


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


def _comparison(plan) -> OECapacityComparisonResult:
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
    if plan.execution_mode in {
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
    )
