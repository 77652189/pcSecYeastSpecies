from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import pytest

from pcsec_pichia.loading import PcSecPichiaInputs, load_pcsec_pichia_inputs
from pcsec_pichia.probe import TargetSpec
from pcsec_pichia.simulation import (
    GrowthTradeoffResult,
    MATLAB_LEGACY_COST_MEDIUM_EXCHANGES,
    ProteinCostSlopeCompatibilityResult,
    SecretionSimulationResult,
    _fixed_growth_bounds,
    _growth_context_warnings,
    _growth_reaction_context,
    _uptake_cost,
    _uptake_cost_status,
    run_growth_tradeoff,
    run_mixed_carbon_objective_probe,
    run_protein_cost_slope_compatibility,
    solve_secretion_capacity,
    summarize_mixed_carbon_objective_result,
    summarize_simulation_result,
)
from pcsec_pichia.targets import load_builtin_targets


REPO_ROOT = Path(__file__).resolve().parents[2]


@lru_cache(maxsize=1)
def _inputs() -> PcSecPichiaInputs:
    return load_pcsec_pichia_inputs(REPO_ROOT)


@lru_cache(maxsize=None)
def _inputs_for_carbon_source(carbon_source_id: str) -> PcSecPichiaInputs:
    return load_pcsec_pichia_inputs(REPO_ROOT, carbon_source_id=carbon_source_id)


@lru_cache(maxsize=1)
def _builtin_targets() -> dict[str, TargetSpec]:
    return {target.target_id: target for target in load_builtin_targets(REPO_ROOT)}


def test_mixed_carbon_objective_probe_requires_explicit_positive_weights() -> None:
    inputs = _inputs()
    opn = _builtin_targets()["OPN_ALPHA_FULL_PROJECT"]

    result = run_mixed_carbon_objective_probe(
        inputs.prepared_model,
        opn,
        inputs.amino_acids,
        inputs.metabolic,
        inputs.secretory,
        inputs.combined,
        target_exchange_ratio=1e-6,
        carbon_weights={},
    )
    summary = summarize_mixed_carbon_objective_result(result)

    assert result.success is False
    assert result.result_status == "draft_mixed_carbon_objective_unavailable"
    assert summary["objective_mode"] == "weighted_carbon_uptake"
    assert "positive carbon objective weight" in summary["message"]
    assert any("does not change the corrected default pipeline" in item for item in summary["warnings"])


def test_growth_reaction_context_tracks_single_and_mixed_carbon_sources() -> None:
    glycerol_context = _growth_reaction_context(_inputs_for_carbon_source("glycerol").prepared_model)
    mixed_context = _growth_reaction_context(_inputs_for_carbon_source("glucose_glycerol").prepared_model)

    assert glycerol_context["growth_reaction_id"] == "BIOMASS_glyc"
    assert glycerol_context["open_growth_reaction_ids"] == ("BIOMASS_glyc",)
    assert glycerol_context["growth_reaction_status"] == "single_growth_reaction"

    assert mixed_context["growth_reaction_id"] == "BIOMASS"
    assert mixed_context["open_growth_reaction_ids"] == ("BIOMASS", "BIOMASS_glyc")
    assert mixed_context["growth_reaction_status"] == "multiple_growth_reactions_selected"
    assert _fixed_growth_bounds(mixed_context, 0.1) == {
        "BIOMASS": (0.1, 0.1),
        "BIOMASS_glyc": (0.0, 0.0),
    }
    assert any("draft mixed-carbon convention" in warning for warning in _growth_context_warnings(mixed_context))


def test_uptake_cost_only_accepts_exchange_uptake_direction() -> None:
    assert _uptake_cost(-2.5) == pytest.approx(2.5)
    assert _uptake_cost(0.0) == pytest.approx(0.0)
    assert _uptake_cost(1.0) is None
    assert _uptake_cost_status(-2.5) == "uptake_flux"
    assert _uptake_cost_status(0.0) == "zero_flux"
    assert _uptake_cost_status(1.0) == "non_uptake_flux"


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
    assert result.lp_sensitivity is not None
    assert len(result.lp_sensitivity["eq_marginals"]) == result.constraint_counts["eq_total"]
    assert len(result.lp_sensitivity["ub_marginals"]) == result.constraint_counts["ub_total"]
    assert "eq_marginals" not in summary
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


def test_opn_cost_slope_compatibility_runs_tiny_fixed_exchange_grid() -> None:
    inputs = _inputs()
    opn = _builtin_targets()["OPN_ALPHA_FULL_PROJECT"]

    result = run_protein_cost_slope_compatibility(
        inputs.prepared_model,
        opn,
        inputs.amino_acids,
        inputs.metabolic,
        inputs.secretory,
        inputs.combined,
        growth_rates=(0.10,),
        secretion_ratios=(5e-7, 1e-6),
        write_ribosome_translation_constraint=True,
        write_misfolding_constraints=True,
    )

    assert isinstance(result, ProteinCostSlopeCompatibilityResult)
    assert result.enabled is True
    assert result.result_status == "draft_matlab_compatible_cost_slope"
    assert len(result.rows) == 2
    assert len(result.glucose_cost_slopes) == 1
    assert result.glucose_cost_slopes[0]["success"] is True
    assert result.glucose_cost_slopes[0]["slope"] is not None
    assert result.medium_compatibility_mode == "corrected"
    assert result.medium_bound_overrides == ()
    for row in result.rows:
        assert row["objective_reaction"] == "Ex_glc_D"
        assert row["target_exchange_ratio"] in {5e-7, 1e-6}
        assert row["success"] is True


def test_cost_slope_matlab_legacy_medium_mode_is_explicit() -> None:
    inputs = _inputs()
    opn = _builtin_targets()["OPN_ALPHA_FULL_PROJECT"]

    result = run_protein_cost_slope_compatibility(
        inputs.prepared_model,
        opn,
        inputs.amino_acids,
        inputs.metabolic,
        inputs.secretory,
        inputs.combined,
        growth_rates=(0.10,),
        secretion_ratios=(5e-7, 1e-6),
        medium_compatibility_mode="matlab_legacy_cost",
    )

    assert result.medium_compatibility_mode == "matlab_legacy_cost"
    assert tuple(item["reaction_id"] for item in result.medium_bound_overrides) == MATLAB_LEGACY_COST_MEDIUM_EXCHANGES
    assert all(item["legacy_lower_bound"] == 0.0 for item in result.medium_bound_overrides)
    assert all(row["medium_compatibility_mode"] == "matlab_legacy_cost" for row in result.rows)
