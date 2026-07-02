from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import pytest
from scipy.io import loadmat

from pcsec_pichia.analysis import analyze_target_protein_lp_attribution, summarize_protein_lp_attribution
from pcsec_pichia.loading import PcSecPichiaInputs, load_pcsec_pichia_inputs
from pcsec_pichia.secretion_plan import build_secretion_plan
from pcsec_pichia.simulation import (
    build_supported_target_model,
    build_target_enzymedata,
    solve_pcsec_maximize,
    solve_secretion_capacity,
)
from pcsec_pichia.targets import load_builtin_targets


REPO_ROOT = Path(__file__).resolve().parents[2]


@lru_cache(maxsize=1)
def _inputs() -> PcSecPichiaInputs:
    return load_pcsec_pichia_inputs(REPO_ROOT)


@lru_cache(maxsize=1)
def _opn():
    return {target.target_id: target for target in load_builtin_targets(REPO_ROOT)}["OPN_ALPHA_FULL_PROJECT"]


@lru_cache(maxsize=1)
def _opn_optional_result():
    inputs = _inputs()
    target = _opn()
    simulation = solve_secretion_capacity(
        inputs.prepared_model,
        target,
        inputs.amino_acids,
        inputs.metabolic,
        inputs.secretory,
        inputs.combined,
        write_ribosome_translation_constraint=True,
        write_misfolding_constraints=True,
    )
    plan = build_secretion_plan(target)
    attribution = analyze_target_protein_lp_attribution(
        target,
        plan,
        simulation.constraint_counts,
        simulation,
        reaction_ids=tuple(inputs.prepared_model.rxns),
    )
    return simulation, summarize_protein_lp_attribution(attribution)


def test_opn_lp_attribution_has_consistent_shapes_and_blocks() -> None:
    simulation, attribution = _opn_optional_result()

    assert simulation.success is True
    assert simulation.lp_sensitivity is not None
    assert len(simulation.lp_sensitivity["eq_marginals"]) == simulation.constraint_counts["eq_total"]
    assert len(simulation.lp_sensitivity["ub_marginals"]) == simulation.constraint_counts["ub_total"]
    assert len(simulation.lp_sensitivity["lower_marginals"]) >= len(_inputs().prepared_model.rxns)
    assert attribution["result_status"] == "draft_lp_sensitivity"
    assert "eq_marginals" not in attribution
    assert "lower_marginals" not in attribution

    blocks = {row["block"] for row in attribution["dominant_constraint_blocks"]}
    assert {"stoichiometric", "secretory_coupling", "misfolding"}.issubset(blocks)
    assert any(row["block"] == "misfolding" and row["nonzero_marginal_count"] > 0 for row in attribution["dominant_constraint_blocks"])
    assert attribution["target_related_fluxes"][-1]["is_target_exchange"] is True


def test_bound_marginal_direction_matches_finite_difference_for_opn() -> None:
    simulation, attribution = _opn_optional_result()
    inputs = _inputs()
    target = _opn()
    candidate = next(
        row for row in attribution["top_bound_marginals"]
        if row["bound_type"] == "lower" and abs(float(row["marginal"])) > 1.0
    )
    reaction_id = str(candidate["reaction_id"])
    index = int(candidate["variable_index_0based"])
    epsilon = 1e-7

    prepared = build_supported_target_model(inputs.prepared_model, target, inputs.amino_acids)
    assert prepared.model is not None
    target_enzymedata = build_target_enzymedata(target, prepared.model, inputs.secretory)
    fixed_model = prepared.model.with_bounds({"BIOMASS": (0.10, 0.10)})
    secretory = inputs.secretory.with_reaction_coefficients(target_enzymedata.reaction_coefficients)
    combined = inputs.combined.with_target(target_enzymedata)
    old_lower = float(fixed_model.lb[index])
    perturbed = fixed_model.with_bounds({reaction_id: (old_lower - epsilon, None)})

    resolved, _counts = solve_pcsec_maximize(
        perturbed,
        str(simulation.exchange_reaction_id),
        metabolic=inputs.metabolic,
        secretory=secretory,
        combined=combined,
        mu=0.10,
        write_ribosome_translation_constraint=True,
        write_misfolding_constraints=True,
    )

    assert resolved.success is True
    actual_delta = float(resolved.objective_value or 0.0) - float(simulation.objective_value or 0.0)
    predicted_delta = float(candidate["marginal"]) * epsilon
    assert actual_delta == pytest.approx(predicted_delta, rel=0.25, abs=1e-8)
    assert actual_delta * predicted_delta >= 0


def test_matlab_protein_cost_tp_result_is_slope_based_not_lp_attribution() -> None:
    matlab_path = REPO_ROOT / "Results" / "Protein_cost_TP" / "all_proteincost_gluPP.mat"
    data = loadmat(matlab_path, squeeze_me=True, struct_as_record=False)

    assert {"all_slope_glc", "all_slope_ribo", "all_glcCost", "all_TP", "all_mu"}.issubset(data)
    assert data["all_slope_glc"].shape == (351, 2)
    assert data["all_slope_ribo"].shape == (351, 2)

    _simulation, attribution = _opn_optional_result()
    assert attribution["result_status"] == "draft_lp_sensitivity"
    assert "top_constraint_marginals" in attribution
    assert "all_slope_glc" not in attribution
    assert "all_glcCost" not in attribution
