from __future__ import annotations

from pathlib import Path

import pytest

from pcsec_pichia.analysis import (
    ProteinCostAnalysisResult,
    analyze_target_protein_lp_attribution,
    analyze_target_protein_cost,
    build_cost_item_table,
    summarize_protein_cost_slope_compatibility,
    summarize_protein_lp_attribution,
    summarize_protein_cost_analysis,
)
from pcsec_pichia.simulation import ProteinCostSlopeCompatibilityResult
from pcsec_pichia.secretion_plan import build_secretion_plan
from pcsec_pichia.simulation import SecretionSimulationResult
from pcsec_pichia.targets import load_builtin_targets, target_spec_from_mapping


REPO_ROOT = Path(__file__).resolve().parents[2]


def _builtin(target_id: str):
    return {target.target_id: target for target in load_builtin_targets(REPO_ROOT)}[target_id]


def _categories(result: ProteinCostAnalysisResult) -> set[str]:
    return {item.category for item in result.cost_items}


def test_opn_protein_cost_analysis_includes_og_burden() -> None:
    target = _builtin("OPN_ALPHA_FULL_PROJECT")
    plan = build_secretion_plan(target)

    result = analyze_target_protein_cost(target, plan)
    summary = summarize_protein_cost_analysis(result)

    assert isinstance(result, ProteinCostAnalysisResult)
    assert result.target_id == "OPN_ALPHA_FULL_PROJECT"
    assert result.route_kind == "opn_like_soluble_secretory"
    assert result.ptm_counts["o_glycosylation_sites"] == 7
    assert "o_glycosylation" in _categories(result)
    og_item = next(item for item in result.cost_items if item.category == "o_glycosylation")
    assert og_item.raw_value > 0
    assert og_item.relative_score > 0
    assert summary["result_status"] == "draft_explanatory"
    assert summary["total_relative_score"] == pytest.approx(100.0, abs=0.01)


def test_hlf_protein_cost_analysis_uses_project_710_sequence_and_ptms() -> None:
    target = _builtin("hLF")

    result = analyze_target_protein_cost(target)
    summary = summarize_protein_cost_analysis(result)

    assert result.route_kind == "soluble_secretory"
    assert result.sequence_features["full_sequence_length"] == 710
    assert result.sequence_features["leader_sequence_length"] == 19
    assert result.sequence_features["signal_peptide_length"] == 19
    assert result.sequence_features["mature_sequence_length"] == 691
    assert result.ptm_counts["disulfide_sites"] == 21
    assert result.ptm_counts["n_glycosylation_sites"] == 4
    assert {"folding_dsb", "n_glycosylation", "misfolding_erad"}.issubset(_categories(result))
    assert any("710aa" in warning for warning in result.warnings)
    assert summary["sequence_features"]["full_sequence_length"] == 710
    assert summary["total_relative_score"] == pytest.approx(100.0, abs=0.01)


def test_custom_target_cost_analysis_is_deterministic_and_warns_for_unset_ptms() -> None:
    target = target_spec_from_mapping(
        {
            "target_id": "CUSTOM_COST",
            "protein_id": "CUSTOM_COST",
            "mature_sequence": "ACDEFGHIKLMNPQRSTVWY",
            "leader_sequence": "MMAA",
            "signal_peptide_sequence": "MM",
            "through_er": True,
            "localization": "e",
            "disulfide_sites": 0,
            "n_glycosylation_sites": 0,
            "o_glycosylation_sites": 0,
        },
        source="request.target_input",
    )

    first = analyze_target_protein_cost(target)
    second = analyze_target_protein_cost(target)

    assert summarize_protein_cost_analysis(first) == summarize_protein_cost_analysis(second)
    assert first.total_relative_score == pytest.approx(100.0, abs=0.01)
    assert any("DSB/NG/OG" in warning for warning in first.warnings)
    assert build_cost_item_table(first)[0]["category"] == "translation"


def test_custom_target_ng_motif_is_not_counted_without_declared_ng_site() -> None:
    target = target_spec_from_mapping(
        {
            "target_id": "CUSTOM_NG_MOTIF_ZERO",
            "protein_id": "CUSTOM_NG_MOTIF_ZERO",
            "mature_sequence": "ANSTA",
            "leader_sequence": "MMAA",
            "signal_peptide_sequence": "MM",
            "through_er": True,
            "localization": "e",
            "disulfide_sites": 0,
            "n_glycosylation_sites": 0,
            "o_glycosylation_sites": 0,
        },
        source="request.target_input",
    )

    result = analyze_target_protein_cost(target)
    ng_item = next(item for item in result.cost_items if item.category == "n_glycosylation")

    assert result.sequence_features["n_glycosylation_motif_positions"] == [2]
    assert ng_item.raw_value == 0.0
    assert ng_item.relative_score == 0.0


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
