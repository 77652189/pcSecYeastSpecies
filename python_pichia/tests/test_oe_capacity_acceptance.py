from __future__ import annotations

import json
from dataclasses import replace

from pcsec_pichia.oe_capacity import (
    OECapacityAcceptanceObservation,
    OECapacityRegressionCheck,
    build_oe_capacity_acceptance_summary,
    write_oe_capacity_acceptance_outputs,
)


def _observation(target_id: str, case_kind: str) -> OECapacityAcceptanceObservation:
    missing = () if case_kind == "executable" else ("nonzero_baseline_formation_flux",)
    return OECapacityAcceptanceObservation(
        target_id=target_id,
        gene_id="G1" if case_kind == "executable" else "G_BOUNDARY",
        case_kind=case_kind,
        elapsed_seconds=12.5,
        screen_status="completed",
        execution_status="gene_level_executable",
        baseline_objective=1.0,
        proxy_objective=1.0,
        gene_capacity_objective=1.0,
        protein_resource_cost_delta=0.0,
        missing_information=missing,
        output_dir=f"local_runs/oe_capacity/round6/{target_id}-{case_kind}",
    )


def test_acceptance_requires_executable_and_boundary_smoke_for_each_target() -> None:
    observations = tuple(
        _observation(target_id, case_kind)
        for target_id in ("hLF", "OPN_ALPHA_FULL_PROJECT")
        for case_kind in ("executable", "boundary")
    )
    summary = build_oe_capacity_acceptance_summary(
        observations,
        coverage_by_target={
            "hLF": {"gene_level_executable": 10},
            "OPN_ALPHA_FULL_PROJECT": {"gene_level_executable": 10},
        },
        regression_checks=(
            OECapacityRegressionCheck("feature_off", True, "focused test passed"),
            OECapacityRegressionCheck("baseline_1x", True, "focused test passed"),
            OECapacityRegressionCheck("legacy_proxy", True, "focused test passed"),
        ),
    )

    assert summary["passed"] is True
    assert all(row["executable_passed"] for row in summary["targets"])
    assert all(row["boundary_passed"] for row in summary["targets"])
    assert summary["model_relative_only"] is True
    assert summary["predicts_absolute_yield"] is False


def test_acceptance_fails_when_one_target_has_no_clean_executable_case() -> None:
    observations = (
        _observation("hLF", "executable"),
        _observation("hLF", "boundary"),
        _observation("OPN_ALPHA_FULL_PROJECT", "boundary"),
    )
    summary = build_oe_capacity_acceptance_summary(
        observations,
        coverage_by_target={},
        regression_checks=(OECapacityRegressionCheck("feature_off", True, "passed"),),
    )

    assert summary["passed"] is False
    opn = next(row for row in summary["targets"] if row["target_id"] == "OPN_ALPHA_FULL_PROJECT")
    assert opn["executable_passed"] is False


def test_acceptance_fails_when_target_coverage_is_missing() -> None:
    observations = tuple(
        _observation(target_id, case_kind)
        for target_id in ("hLF", "OPN_ALPHA_FULL_PROJECT")
        for case_kind in ("executable", "boundary")
    )
    summary = build_oe_capacity_acceptance_summary(
        observations,
        coverage_by_target={"hLF": {"gene_level_executable": 10}},
        regression_checks=(OECapacityRegressionCheck("feature_off", True, "passed"),),
    )

    assert summary["passed"] is False
    opn = next(row for row in summary["targets"] if row["target_id"] == "OPN_ALPHA_FULL_PROJECT")
    assert opn["coverage_present"] is False


def test_acceptance_requires_complete_baseline_proxy_capacity_comparison() -> None:
    incomplete = replace(
        _observation("hLF", "executable"),
        proxy_objective=None,
    )
    observations = (
        incomplete,
        _observation("hLF", "boundary"),
        _observation("OPN_ALPHA_FULL_PROJECT", "executable"),
        _observation("OPN_ALPHA_FULL_PROJECT", "boundary"),
    )
    summary = build_oe_capacity_acceptance_summary(
        observations,
        coverage_by_target={
            "hLF": {"gene_level_executable": 10},
            "OPN_ALPHA_FULL_PROJECT": {"gene_level_executable": 10},
        },
        regression_checks=(OECapacityRegressionCheck("legacy_proxy", True, "passed"),),
    )

    assert summary["passed"] is False
    hlf = next(row for row in summary["targets"] if row["target_id"] == "hLF")
    assert hlf["executable_passed"] is False


def test_acceptance_outputs_are_local_auditable_json_and_markdown(tmp_path) -> None:
    summary = build_oe_capacity_acceptance_summary(
        tuple(
            _observation(target_id, case_kind)
            for target_id in ("hLF", "OPN_ALPHA_FULL_PROJECT")
            for case_kind in ("executable", "boundary")
        ),
        coverage_by_target={
            "hLF": {"gene_level_executable": 10},
            "OPN_ALPHA_FULL_PROJECT": {"gene_level_executable": 10},
        },
        regression_checks=(OECapacityRegressionCheck("protected_paths", True, "empty diff"),),
    )

    json_path, markdown_path = write_oe_capacity_acceptance_outputs(summary, tmp_path)

    payload = json.loads(json_path.read_text(encoding="utf-8"))
    assert payload["passed"] is True
    assert markdown_path.read_text(encoding="utf-8").startswith(
        "# Phase 2 gene-level OE capacity acceptance"
    )
