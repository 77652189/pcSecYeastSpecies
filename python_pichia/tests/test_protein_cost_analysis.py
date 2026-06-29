from __future__ import annotations

from pathlib import Path

import pytest

from pcsec_pichia.analysis import (
    ProteinCostAnalysisResult,
    analyze_target_protein_cost,
    build_cost_item_table,
    summarize_protein_cost_analysis,
)
from pcsec_pichia.secretion_plan import build_secretion_plan
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
