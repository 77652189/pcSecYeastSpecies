from __future__ import annotations

from pathlib import Path

from pcsec_pichia.analysis import (
    ProteinLpAttributionResult,
    analyze_target_protein_lp_attribution,
    compare_solver_robustness,
    summarize_protein_cost_slope_compatibility,
    summarize_protein_lp_attribution,
    summarize_solver_robustness,
)
from pcsec_pichia.simulation import ProteinCostSlopeCompatibilityResult
from pcsec_pichia.secretion_plan import build_secretion_plan
from pcsec_pichia.simulation import SecretionSimulationResult
from pcsec_pichia.targets import load_builtin_targets


REPO_ROOT = Path(__file__).resolve().parents[2]


def _builtin(target_id: str):
    return {target.target_id: target for target in load_builtin_targets(REPO_ROOT)}[target_id]


def test_lp_attribution_summarizes_top_n_without_full_marginal_arrays() -> None:
    target = _builtin("OPN_ALPHA_FULL_PROJECT")
    plan = build_secretion_plan(target)
    simulation = SecretionSimulationResult(
        success=True,
        target_id=target.target_id,
        objective_value=0.006,
        growth_rate=0.10,
        secretion_flux=0.006,
        status="0",
        message="ok",
        constraint_counts={
            "stoichiometric": 2,
            "secretory_coupling": 2,
            "protein_mass": 1,
            "ribosome_translation": 1,
            "misfolding": 1,
            "mitochondrial": 1,
            "eq_total": 7,
            "ub_total": 1,
        },
        result_status="draft",
        target_parameter_status="draft",
        matlab_alignment_status="pending",
        exchange_reaction_id="r_OPN_ALPHA_FULL_PROJECT_exchange",
        build_status="supported",
        lp_sensitivity={
            "eq_marginals": (0.0, 0.2, -0.5, 0.1, 0.7, -0.9, 0.3),
            "ub_marginals": (0.4,),
            "lower_marginals": (0.0, -0.8, 0.0, 0.0),
            "upper_marginals": (0.0, 0.0, 0.6, 0.0),
        },
        key_fluxes={
            "BIOMASS": 0.10,
            "r_OPN_ALPHA_FULL_PROJECT_exchange": 0.006,
        },
    )

    result = analyze_target_protein_lp_attribution(
        target,
        plan,
        simulation.constraint_counts,
        simulation,
        reaction_ids=("BIOMASS", "r_OPN_ALPHA_FULL_PROJECT_exchange", "Ex_glc_D", "OTHER"),
        top_n=3,
    )
    summary = summarize_protein_lp_attribution(result)

    assert summary["result_status"] == "draft_lp_sensitivity"
    assert len(summary["top_constraint_marginals"]) == 3
    assert summary["top_constraint_marginals"][0]["block"] == "ribosome_translation"
    assert len(summary["top_bound_marginals"]) == 2
    assert summary["top_bound_marginals"][0]["reaction_id"] == "r_OPN_ALPHA_FULL_PROJECT_exchange"
    assert summary["active_bound_counts"]["total_bound_marginal_nonzero"] == 2
    assert "eq_marginals" not in summary
    assert "lower_marginals" not in summary


def test_oe_actionable_bottlenecks_exclude_lower_bound_floor() -> None:
    # R1 (ADR-004): the biggest bound marginal here is a LOWER bound (floor, abs 0.8) that
    # OE cannot relax; the OE-addressable one is a smaller UPPER bound (abs 0.6). The derived
    # oe_actionable list must contain only the upper bound and never the larger lower one -
    # this is the PDI1-alone / ribosome false-positive guard made programmatic.
    target = _builtin("OPN_ALPHA_FULL_PROJECT")
    plan = build_secretion_plan(target)
    simulation = SecretionSimulationResult(
        success=True,
        target_id=target.target_id,
        objective_value=0.006,
        growth_rate=0.10,
        secretion_flux=0.006,
        status="0",
        message="ok",
        constraint_counts={"stoichiometric": 7, "eq_total": 7, "ub_total": 1},
        result_status="draft",
        target_parameter_status="draft",
        matlab_alignment_status="pending",
        exchange_reaction_id="r_OPN_ALPHA_FULL_PROJECT_exchange",
        build_status="supported",
        lp_sensitivity={
            "eq_marginals": (0.0,) * 7,
            "ub_marginals": (0.0,),
            "lower_marginals": (0.0, -0.8, 0.0, 0.0),
            "upper_marginals": (0.0, 0.0, 0.6, 0.0),
        },
        key_fluxes={"BIOMASS": 0.10},
    )
    reaction_ids = ("BIOMASS", "r_OPN_ALPHA_FULL_PROJECT_exchange", "Ex_glc_D", "OTHER")

    result = analyze_target_protein_lp_attribution(
        target, plan, simulation.constraint_counts, simulation, reaction_ids=reaction_ids, top_n=5
    )

    oe_ids = [row["reaction_id"] for row in result.oe_actionable_bottlenecks]
    floor_ids = [row["reaction_id"] for row in result.floor_constraints_not_oe_addressable]
    assert oe_ids == ["Ex_glc_D"]
    assert all(row["bound_type"] == "upper" and row["oe_actionable"] is True for row in result.oe_actionable_bottlenecks)
    # the larger lower-bound floor is segregated, never promoted into the OE-actionable list
    assert "r_OPN_ALPHA_FULL_PROJECT_exchange" in floor_ids
    assert "r_OPN_ALPHA_FULL_PROJECT_exchange" not in oe_ids
    assert all(row["bound_type"] == "lower" for row in result.floor_constraints_not_oe_addressable)
    summary = summarize_protein_lp_attribution(result)
    assert summary["oe_actionable_bottlenecks"] == list(result.oe_actionable_bottlenecks)
    assert summary["floor_constraints_not_oe_addressable"] == list(result.floor_constraints_not_oe_addressable)


def _attr_with_top_bottleneck(reaction_id: str, block: str) -> ProteinLpAttributionResult:
    return ProteinLpAttributionResult(
        target_id="hLF",
        result_status="draft_lp_sensitivity",
        objective_evidence={},
        dominant_constraint_blocks=({"block": block, "sum_abs_marginal": 1.0},),
        top_constraint_marginals=(),
        top_bound_marginals=(),
        target_related_fluxes=(),
        active_bound_counts={},
        warnings=(),
        oe_actionable_bottlenecks=({"reaction_id": reaction_id, "bound_type": "upper", "abs_marginal": 1.0},),
    )


def test_compare_solver_robustness_flags_solver_dependent_bottleneck() -> None:
    stable = compare_solver_robustness(
        "hLF",
        {
            "highs": _attr_with_top_bottleneck("sec_X_complex_formation", "secretory_coupling"),
            "highs-ds": _attr_with_top_bottleneck("sec_X_complex_formation", "secretory_coupling"),
        },
    )
    assert stable.classification == "ranking-insensitive-to-solver"
    assert stable.top_bottleneck_stable is True

    flipped = compare_solver_robustness(
        "hLF",
        {
            "highs": _attr_with_top_bottleneck("sec_X_complex_formation", "secretory_coupling"),
            "highs-ipm": _attr_with_top_bottleneck("sec_Y_complex_formation", "secretory_coupling"),
        },
    )
    assert flipped.classification == "ranking-sensitive-to-solver"
    assert flipped.top_bottleneck_stable is False
    assert summarize_solver_robustness(flipped)["classification"] == "ranking-sensitive-to-solver"

    inconclusive = compare_solver_robustness(
        "hLF", {"highs": _attr_with_top_bottleneck("sec_X_complex_formation", "secretory_coupling")}
    )
    assert inconclusive.classification == "inconclusive"

    # a re-solve that errored is fed in as an unavailable attribution: even with a good
    # method present, the verdict degrades to inconclusive rather than pretending stability.
    errored = ProteinLpAttributionResult(
        target_id="hLF",
        result_status="draft_lp_sensitivity_unavailable",
        objective_evidence={},
        dominant_constraint_blocks=(),
        top_constraint_marginals=(),
        top_bound_marginals=(),
        target_related_fluxes=(),
        active_bound_counts={},
        warnings=(),
    )
    with_error = compare_solver_robustness(
        "hLF",
        {
            "highs": _attr_with_top_bottleneck("sec_X_complex_formation", "secretory_coupling"),
            "highs-ipm": errored,
        },
    )
    assert with_error.classification == "inconclusive"
    assert any(row["result_status"] == "draft_lp_sensitivity_unavailable" for row in with_error.per_method)


def test_lp_attribution_handles_missing_sensitivity_without_crashing() -> None:
    target = _builtin("hLF")
    plan = build_secretion_plan(target)
    simulation = SecretionSimulationResult(
        success=False,
        target_id=target.target_id,
        objective_value=None,
        growth_rate=0.10,
        secretion_flux=None,
        status="2",
        message="infeasible",
        constraint_counts={},
        result_status="draft",
        target_parameter_status="draft_matlab_alignment_pending",
        matlab_alignment_status="pending",
        exchange_reaction_id=None,
        build_status="failed",
    )

    result = analyze_target_protein_lp_attribution(target, plan, {}, simulation)

    assert result.result_status == "draft_lp_sensitivity_unavailable"
    assert result.top_constraint_marginals == ()
    assert any("unavailable" in warning for warning in result.warnings)


def test_cost_slope_compatibility_summary_is_disabled_by_default_payload() -> None:
    summary = summarize_protein_cost_slope_compatibility(None)

    assert summary["enabled"] is False
    assert summary["result_status"] == "disabled"


def test_cost_slope_compatibility_summary_keeps_matlab_style_definition() -> None:
    result = ProteinCostSlopeCompatibilityResult(
        target_id="OPN_ALPHA_FULL_PROJECT",
        enabled=True,
        success=True,
        growth_rates=(0.05,),
        secretion_ratios=(5e-7, 1e-6),
        rows=(
            {
                "mu": 0.05,
                "target_exchange_ratio": 5e-7,
                "success": True,
                "glucose_cost": 1.0,
            },
        ),
        glucose_cost_slopes=(
            {
                "mu": 0.05,
                "cost_key": "glucose_cost",
                "success": True,
                "slope": 10.0,
                "point_count": 2,
                "status": "slope_estimated",
            },
        ),
        ribosome_cost_slopes=(),
        result_status="draft_matlab_compatible_cost_slope",
        warnings=("draft",),
    )

    summary = summarize_protein_cost_slope_compatibility(result)

    assert summary["enabled"] is True
    assert summary["result_status"] == "draft_matlab_compatible_cost_slope"
    assert summary["medium_compatibility_mode"] == "corrected"
    assert summary["medium_bound_overrides"] == []
    assert summary["comparison_scope"]["not_default_pipeline"] is True
    assert summary["comparison_scope"]["medium_compatibility"] == "corrected"
    assert summary["comparison_scope"]["current_default_definition"] == (
        "fixed growth rate, corrected medium, maximize target secretion flux"
    )
    assert summary["comparison_scope"]["ratio_policy"] == "explicit_absolute_ratios"
    assert summary["secretion_ratio_policy"] == "explicit_absolute_ratios"
    assert "fixed target exchange" in summary["comparison_scope"]["definition"]
