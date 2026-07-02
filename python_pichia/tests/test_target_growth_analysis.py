from __future__ import annotations

from pcsec_pichia.analysis import (
    TargetGrowthAnalysisResult,
    analyze_target_growth_impact,
    build_growth_tradeoff_item_table,
    summarize_target_growth_analysis,
)
from pcsec_pichia.simulation import GrowthTradeoffResult


def _tradeoff(target_id: str = "OPN_ALPHA_FULL_PROJECT") -> GrowthTradeoffResult:
    return GrowthTradeoffResult(
        target_id=target_id,
        success=False,
        growth_points=(0.05, 0.10, 0.20),
        tradeoff_rows=(
            {
                "mu": 0.05,
                "success": True,
                "status": "0",
                "secretion_flux": 0.004,
                "secretion_per_biomass": 0.08,
                "message": "ok",
            },
            {
                "mu": 0.10,
                "success": True,
                "status": "0",
                "secretion_flux": 0.006,
                "secretion_per_biomass": 0.06,
                "message": "ok",
            },
            {
                "mu": 0.20,
                "success": False,
                "status": "2",
                "secretion_flux": None,
                "secretion_per_biomass": None,
                "message": "infeasible",
            },
        ),
        result_status="draft",
        target_parameter_status="draft",
        matlab_alignment_status="pending",
    )


def test_growth_analysis_selects_best_points_and_keeps_failure_rows() -> None:
    result = analyze_target_growth_impact(_tradeoff(), baseline_growth_rate=0.10)
    summary = summarize_target_growth_analysis(result)

    assert isinstance(result, TargetGrowthAnalysisResult)
    assert result.target_id == "OPN_ALPHA_FULL_PROJECT"
    assert result.valid_point_count == 2
    assert result.growth_sensitivity_label == "mixed"
    assert result.growth_sensitivity_reason == "contains_failed_or_missing_points"
    assert result.best_secretion_point == {
        "mu": 0.1,
        "success": True,
        "secretion_flux": 0.006,
        "secretion_per_biomass": 0.06,
        "status": "0",
        "interpretation": "最高分泌通量；基准生长点",
    }
    assert result.best_secretion_per_biomass_point == {
        "mu": 0.05,
        "success": True,
        "secretion_flux": 0.004,
        "secretion_per_biomass": 0.08,
        "status": "0",
        "interpretation": "最高单位生物量分泌",
    }
    assert summary["result_status"] == "draft_explanatory"
    assert summary["growth_sensitivity_reason"] == "contains_failed_or_missing_points"
    assert len(summary["tradeoff_points"]) == 3
    assert summary["tradeoff_points"][1]["interpretation"] == "最高分泌通量；基准生长点"
    assert summary["tradeoff_points"][2]["interpretation"] == "求解失败或无可比较分泌值"


def test_growth_analysis_handles_empty_tradeoff_rows() -> None:
    result = analyze_target_growth_impact(
        GrowthTradeoffResult(
            target_id="hLF",
            success=False,
            growth_points=(0.10,),
            tradeoff_rows=(),
            result_status="draft",
            target_parameter_status="draft_matlab_alignment_pending",
            matlab_alignment_status="pending",
        )
    )

    assert result.target_id == "hLF"
    assert result.valid_point_count == 0
    assert result.best_secretion_point is None
    assert result.best_secretion_per_biomass_point is None
    assert result.growth_sensitivity_label == "insufficient_points"
    assert result.growth_sensitivity_reason == "insufficient_tradeoff_rows"
    assert build_growth_tradeoff_item_table(result) == ()


def test_growth_analysis_classifies_mixed_tradeoff() -> None:
    tradeoff = GrowthTradeoffResult(
        target_id="CUSTOM",
        success=True,
        growth_points=(0.05, 0.10, 0.15),
        tradeoff_rows=(
            {"mu": 0.05, "success": True, "status": "0", "secretion_flux": 0.004, "secretion_per_biomass": 0.08},
            {"mu": 0.10, "success": True, "status": "0", "secretion_flux": 0.006, "secretion_per_biomass": 0.06},
            {"mu": 0.15, "success": True, "status": "0", "secretion_flux": 0.005, "secretion_per_biomass": 0.033},
        ),
        result_status="draft",
        target_parameter_status="draft",
        matlab_alignment_status="pending",
    )

    result = analyze_target_growth_impact(tradeoff)

    assert result.growth_sensitivity_label == "mixed"
    assert result.growth_sensitivity_reason == "non_monotonic_successful_grid"
    assert result.best_secretion_point is not None
    assert result.best_secretion_point["mu"] == 0.10


def test_growth_analysis_parses_string_false_as_failure() -> None:
    tradeoff = GrowthTradeoffResult(
        target_id="CSV_STYLE_BOOL",
        success=False,
        growth_points=(0.05, 0.10),
        tradeoff_rows=(
            {
                "mu": 0.05,
                "success": "True",
                "status": "0",
                "secretion_flux": 0.004,
                "secretion_per_biomass": 0.08,
            },
            {
                "mu": 0.10,
                "success": "False",
                "status": "2",
                "secretion_flux": 0.006,
                "secretion_per_biomass": 0.06,
            },
        ),
        result_status="draft",
        target_parameter_status="draft",
        matlab_alignment_status="pending",
    )

    result = analyze_target_growth_impact(tradeoff)

    assert result.valid_point_count == 1
    assert result.tradeoff_points[1].success is False
    assert result.growth_sensitivity_label == "mixed"
    assert result.growth_sensitivity_reason == "contains_failed_or_missing_points"


def test_growth_analysis_does_not_call_discontinuous_grid_increasing() -> None:
    tradeoff = GrowthTradeoffResult(
        target_id="CUSTOM_FAIL_MIDDLE",
        success=False,
        growth_points=(0.05, 0.10, 0.20),
        tradeoff_rows=(
            {"mu": 0.05, "success": True, "status": "0", "secretion_flux": 0.004, "secretion_per_biomass": 0.08},
            {"mu": 0.10, "success": False, "status": "2", "secretion_flux": None, "secretion_per_biomass": None},
            {"mu": 0.20, "success": True, "status": "0", "secretion_flux": 0.006, "secretion_per_biomass": 0.03},
        ),
        result_status="draft",
        target_parameter_status="draft",
        matlab_alignment_status="pending",
    )

    result = analyze_target_growth_impact(tradeoff)

    assert result.growth_sensitivity_label == "mixed"
    assert result.growth_sensitivity_reason == "contains_failed_or_missing_points"
    assert result.valid_point_count == 2
    assert result.tradeoff_points[1].interpretation == "求解失败或无可比较分泌值"
