from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import pytest

from pcsec_pichia.loading import PcSecPichiaInputs, load_pcsec_pichia_inputs
from pcsec_pichia.probe import TargetSpec
from pcsec_pichia.simulation import (
    GrowthTradeoffResult,
    SecretionSimulationResult,
    run_growth_tradeoff,
    solve_secretion_capacity,
    summarize_simulation_result,
)
from pcsec_pichia.targets import load_builtin_targets


REPO_ROOT = Path(__file__).resolve().parents[2]


@lru_cache(maxsize=1)
def _inputs() -> PcSecPichiaInputs:
    return load_pcsec_pichia_inputs(REPO_ROOT)


@lru_cache(maxsize=1)
def _builtin_targets() -> dict[str, TargetSpec]:
    return {target.target_id: target for target in load_builtin_targets(REPO_ROOT)}


@pytest.mark.parametrize(
    ("target_id", "expected_objective", "expected_status"),
    (
        ("OPN_ALPHA_FULL_PROJECT", 0.006572021526431409, "draft"),
        ("hLF", 0.0032850100270232106, "draft_matlab_alignment_pending"),
    ),
)
def test_builtin_targets_solve_default_secretion_capacity(
    target_id: str,
    expected_objective: float,
    expected_status: str,
) -> None:
    inputs = _inputs()
    target = _builtin_targets()[target_id]

    result = solve_secretion_capacity(
        inputs.prepared_model,
        target,
        inputs.amino_acids,
        inputs.metabolic,
        inputs.secretory,
        inputs.combined,
    )
    summary = summarize_simulation_result(result)

    assert isinstance(result, SecretionSimulationResult)
    assert result.success is True
    assert result.status == "0"
    assert result.objective_value == pytest.approx(expected_objective)
    assert result.secretion_flux == pytest.approx(expected_objective)
    assert result.growth_rate == pytest.approx(0.10)
    assert result.result_status == "draft"
    assert result.target_parameter_status == expected_status
    assert result.matlab_alignment_status == "pending"
    assert result.constraint_counts["eq_total"] > 0
    assert result.constraint_counts["ub_total"] == 1
    assert summary["target_id"] == target_id
    assert summary["objective_value"] == pytest.approx(expected_objective)


@pytest.mark.parametrize(
    ("target_id", "expected_objective", "expected_eq_total"),
    (
        ("OPN_ALPHA_FULL_PROJECT", 0.0021305196599992996, 24434),
        ("hLF", 0.001112112385054876, 24443),
    ),
)
def test_builtin_targets_solve_optional_constraint_secretion_capacity(
    target_id: str,
    expected_objective: float,
    expected_eq_total: int,
) -> None:
    inputs = _inputs()
    target = _builtin_targets()[target_id]

    result = solve_secretion_capacity(
        inputs.prepared_model,
        target,
        inputs.amino_acids,
        inputs.metabolic,
        inputs.secretory,
        inputs.combined,
        write_ribosome_translation_constraint=True,
        write_misfolding_constraints=True,
    )

    assert result.success is True
    assert result.objective_value == pytest.approx(expected_objective)
    assert result.constraint_counts["ribosome_translation"] == 1
    assert result.constraint_counts["misfolding"] == 1418
    assert result.constraint_counts["eq_total"] == expected_eq_total
    assert result.result_status == "draft"
    assert result.matlab_alignment_status == "pending"


def test_opn_growth_tradeoff_runs_small_grid_without_screens() -> None:
    inputs = _inputs()
    opn = _builtin_targets()["OPN_ALPHA_FULL_PROJECT"]

    result = run_growth_tradeoff(
        inputs.prepared_model,
        opn,
        inputs.amino_acids,
        inputs.metabolic,
        inputs.secretory,
        inputs.combined,
        growth_points=(0.05, 0.10),
    )
    summary = summarize_simulation_result(result)

    assert isinstance(result, GrowthTradeoffResult)
    assert result.success is True
    assert result.growth_points == (0.05, 0.10)
    assert len(result.tradeoff_rows) == 2
    assert result.result_status == "draft"
    assert result.matlab_alignment_status == "pending"
    for row in result.tradeoff_rows:
        assert row["success"] is True
        assert row["secretion_flux"] is not None
        assert row["constraint_counts"]["eq_total"] > 0
    assert summary["tradeoff_rows"] == result.tradeoff_rows
