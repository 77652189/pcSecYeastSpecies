from __future__ import annotations

import importlib.util
import json
from pathlib import Path

from pcsec_pichia.screens.genome_wide_tradeoff import (
    _classify_solve_outcome,
    _representative_failure_row,
    _summarize_catalog_row,
    _tradeoff_point,
)

SOLVER_OUTCOME_CSV_FIELDS = (
    "solve_outcome_counts",
    "has_timeout",
    "timeout_mu_points",
    "proven_infeasible_mu_points",
    "other_solver_failure_mu_points",
    "feasibility_interpretation",
)


def _load_tool_module(module_name: str):
    tool_path = Path(__file__).resolve().parents[1] / "tools" / f"{module_name}.py"
    spec = importlib.util.spec_from_file_location(module_name, tool_path)
    assert spec is not None
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_timeout_above_success_keeps_max_mu_but_marks_inconclusive() -> None:
    points = [
        _tradeoff_point(0.05, True, "0", 1.25, "Optimization terminated successfully."),
        _tradeoff_point(0.10, False, "1", None, "Time limit reached"),
    ]

    row = _summarize_catalog_row("R_TIMEOUT", "Timeout candidate", "test", "KO", points, None)

    assert row["max_feasible_mu"] == 0.05
    assert row["secretion_at_max_feasible_mu"] == 1.25
    assert row["tradeoff_points"] == tuple(points)
    assert row["solve_outcome_counts"]["success"] == 1
    assert row["solve_outcome_counts"]["time_limit_reached"] == 1
    assert row["has_timeout"] is True
    assert row["timeout_mu_points"] == (0.10,)
    assert row["proven_infeasible_mu_points"] == ()
    assert row["feasibility_interpretation"] == "inconclusive_due_to_timeout"


def test_status_two_marks_proven_infeasible_not_timeout() -> None:
    points = [
        _tradeoff_point(0.05, True, "0", 1.25),
        _tradeoff_point(0.10, False, "2", None, "The problem is infeasible."),
    ]

    row = _summarize_catalog_row("R_INFEASIBLE", "Infeasible candidate", "test", "KO", points, None)

    assert row["max_feasible_mu"] == 0.05
    assert row["has_timeout"] is False
    assert row["timeout_mu_points"] == ()
    assert row["proven_infeasible_mu_points"] == (0.10,)
    assert row["solve_outcome_counts"]["proven_infeasible"] == 1
    assert row["feasibility_interpretation"] == "definitive"


def test_timeout_classification_accepts_status_one_and_highs_status_message() -> None:
    assert _classify_solve_outcome(False, 1, "") == "time_limit_reached"
    assert _classify_solve_outcome(False, 4, "HiGHS Status 13: Time limit reached") == "time_limit_reached"
    assert _classify_solve_outcome(False, 2, "infeasible") == "proven_infeasible"
    assert _classify_solve_outcome(False, 4, "numerical failure") == "other_solver_failure"


def test_oe_failure_selection_preserves_timeout_over_infeasible() -> None:
    failure_row = _representative_failure_row(
        [
            {"success": False, "status": 2, "objective_value": None},
            {"success": False, "status": 1, "objective_value": None},
        ]
    )

    assert failure_row is not None
    assert failure_row["status"] == 1


def test_summary_preserves_existing_row_fields() -> None:
    points = [_tradeoff_point(0.05, True, "0", 1.25)]

    row = _summarize_catalog_row("R_OK", "Compatible candidate", "test", "OE", points, None)

    for key in (
        "gene_id",
        "intervention_type",
        "affected_reactions",
        "support_status",
        "max_feasible_mu",
        "secretion_at_max_feasible_mu",
        "tradeoff_points",
        "skipped_reason",
    ):
        assert key in row
    assert row["tradeoff_points"][0]["solve_outcome"] == "success"
    assert row["feasibility_interpretation"] == "definitive"


def test_cli_csv_records_include_solver_outcome_summary_fields() -> None:
    row = {
        "target_id": "hLF",
        "gene_id": "PAS_chr1-1_0001",
        "common_name": "candidate",
        "candidate_kind": "gene",
        "intervention_type": "KO",
        "support_status": "gpr_supported",
        "secretory_process": "folding",
        "gpr_role": "single_gene",
        "mapping_confidence": "model_gpr",
        "max_feasible_mu": 0.05,
        "secretion_at_max_feasible_mu": 1.25,
        "wildtype_max_feasible_mu": 0.10,
        "wildtype_secretion_at_max_feasible_mu": 2.0,
        "growth_retention_ratio": 0.5,
        "secretion_ratio_vs_wildtype": 0.625,
        "solve_outcome_counts": {"success": 1, "time_limit_reached": 1},
        "has_timeout": True,
        "timeout_mu_points": (0.10,),
        "proven_infeasible_mu_points": (),
        "other_solver_failure_mu_points": (),
        "feasibility_interpretation": "inconclusive_due_to_timeout",
        "affected_reactions": ("R_TIMEOUT",),
        "skipped_reason": None,
        "hypothesis_note": "",
    }

    for module_name in ("run_genome_wide_ko_oe_screen", "run_genome_wide_ko_oe_screen_parallel"):
        module = _load_tool_module(module_name)
        for field in SOLVER_OUTCOME_CSV_FIELDS:
            assert field in module.CSV_FIELDS

        record = module._row_to_csv_record(row)

        assert record["has_timeout"] is True
        assert record["feasibility_interpretation"] == "inconclusive_due_to_timeout"
        assert json.loads(record["solve_outcome_counts"]) == {"success": 1, "time_limit_reached": 1}
        assert json.loads(record["timeout_mu_points"]) == [0.10]
        assert json.loads(record["proven_infeasible_mu_points"]) == []
        assert json.loads(record["other_solver_failure_mu_points"]) == []
