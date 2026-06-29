from __future__ import annotations

import ast
import os
from functools import lru_cache
from pathlib import Path

import pytest

from pcsec_pichia.loading import PcSecPichiaInputs, load_pcsec_pichia_inputs
from pcsec_pichia.screens import (
    ScreenResult,
    build_gene_perturbation_map,
    run_knockout_screen,
    run_overexpression_screen,
    run_reaction_knockout_screen,
    summarize_screen_result,
)
from pcsec_pichia.targets import TargetSpec, load_builtin_targets


REPO_ROOT = Path(__file__).resolve().parents[2]
KO_GENES = ["AT250_GQ_6803479", "AT250_GQ_6803809"]
OE_REACTIONS = ["sec_BIP_NEFS_complex_formation", "sec_Kar2p_complex_formation"]
KO_REACTIONS = ["sec_Och1p_complex_formation"]
slow_screen = pytest.mark.skipif(
    os.environ.get("PCSEC_RUN_SLOW_SCREEN_TESTS") != "1",
    reason="slow pcSec screen solve; set PCSEC_RUN_SLOW_SCREEN_TESTS=1 to run",
)
REQUIRED_ROW_FIELDS = {
    "candidate_id",
    "intervention_type",
    "objective_value",
    "baseline_objective_value",
    "delta_objective",
    "input_gene_id",
    "resolved_reaction_id",
    "effect_label",
    "solver_status_label",
    "failure_reason",
    "secretory_process",
    "mapping_level",
    "mapping_confidence",
    "mapping_interpretation",
    "complex_id",
    "complex_subunit_ids",
    "complex_subunit_stoichiometry",
}


def test_gene_perturbation_map_explains_gene_reaction_process_and_confidence() -> None:
    class TinyModel:
        rxns = [
            "sec_Kar2p_complex_formation",
            "sec_PDI1_ERV2_Ero1p_complex_formation",
            "sec_Kar2p_activity",
            "BIOMASS",
        ]
        rules = ["x(1)", "x(1)", "x(2)", "x(3)"]
        gr_rules = ["G1", "G1", "G2", "G3"]
        gene_index = {"G1": 0, "G2": 1, "G3": 2}
        reaction_index = {
            "sec_Kar2p_complex_formation": 0,
            "sec_PDI1_ERV2_Ero1p_complex_formation": 1,
            "sec_Kar2p_activity": 2,
            "BIOMASS": 3,
        }

    result = build_gene_perturbation_map(
        TinyModel(),
        ("G1", "G2", "G3", "NO_SUCH_GENE"),
        complex_subunits={
            "sec_Kar2p_complex": [
                {"subunit_id": "Kar2p", "stoichiometry": 1.0},
            ],
        },
    )
    rows = result.rows

    g1_rows = [row for row in rows if row["gene_id"] == "G1"]
    assert len(g1_rows) == 2
    assert all(row["reaction_count"] == 2 for row in g1_rows)
    assert any(row["secretory_process"] == "ER 折叠 / 分子伴侣" for row in g1_rows)

    complex_row = next(row for row in rows if row["reaction_id"] == "sec_Kar2p_complex_formation")
    assert complex_row["mapping_level"] == "complex_subunit"
    assert complex_row["mapping_confidence"] == "medium"
    assert complex_row["complex_subunit_ids"] == ["Kar2p"]
    assert complex_row["complex_subunit_stoichiometry"] == [1.0]

    direct_row = next(row for row in rows if row["reaction_id"] == "sec_Kar2p_activity")
    assert direct_row["mapping_level"] == "direct_gpr"
    assert direct_row["mapping_confidence"] == "high"

    metabolic_row = next(row for row in rows if row["reaction_id"] == "BIOMASS")
    assert metabolic_row["mapping_level"] == "metabolic_or_other"
    assert metabolic_row["mapping_confidence"] == "low"

    unresolved_row = next(row for row in rows if row["gene_id"] == "NO_SUCH_GENE")
    assert unresolved_row["mapping_level"] == "unresolved"
    assert unresolved_row["mapping_confidence"] == "unresolved"
    assert unresolved_row["resolved"] is False


def test_formal_screen_row_normalization_adds_mapping_explanation_fields() -> None:
    import pcsec_pichia.screens as screens_module

    row = screens_module._normalize_screen_row(
        {
            "reaction": "sec_Kar2p_complex_formation",
            "success": True,
            "status": "0",
            "objective_value": 0.002,
            "delta_vs_baseline": 0.0001,
        },
        target_id="OPN_ALPHA_FULL_PROJECT",
        screen_type="overexpression",
        intervention_type="OE_reaction",
        baseline_objective_value=0.0019,
        complex_subunits={
            "sec_Kar2p_complex": [
                {"subunit_id": "Kar2p", "stoichiometry": 1.0},
                {"subunit_id": "Sil1p", "stoichiometry": 1.0},
            ],
        },
    )

    assert row["mapping_level"] == "complex_subunit"
    assert row["mapping_confidence"] == "medium"
    assert row["mapping_interpretation"]
    assert row["complex_id"] == "sec_Kar2p_complex"
    assert row["complex_subunit_ids"] == ["Kar2p", "Sil1p"]
    assert row["secretory_process"] == "ER 折叠 / 分子伴侣"


def test_screen_solve_tests_are_slow_gated() -> None:
    module_ast = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    solve_calls = {"run_knockout_screen", "run_overexpression_screen", "run_reaction_knockout_screen"}
    ungated: list[str] = []
    for node in module_ast.body:
        if not isinstance(node, ast.FunctionDef) or not node.name.startswith("test_"):
            continue
        calls_screen_solve = any(
            isinstance(child, ast.Call)
            and isinstance(child.func, ast.Name)
            and child.func.id in solve_calls
            for child in ast.walk(node)
        )
        has_slow_marker = any(
            isinstance(decorator, ast.Name) and decorator.id == "slow_screen"
            for decorator in node.decorator_list
        )
        if calls_screen_solve and not has_slow_marker:
            ungated.append(node.name)

    assert ungated == []


@lru_cache(maxsize=1)
def _inputs() -> PcSecPichiaInputs:
    return load_pcsec_pichia_inputs(REPO_ROOT)


@lru_cache(maxsize=1)
def _builtin_targets() -> dict[str, TargetSpec]:
    return {target.target_id: target for target in load_builtin_targets(REPO_ROOT)}


def _assert_common_screen_result(result: ScreenResult, target_id: str, screen_type: str) -> None:
    summary = summarize_screen_result(result)

    assert isinstance(result, ScreenResult)
    assert result.target_id == target_id
    assert result.screen_type == screen_type
    assert result.success is True
    assert result.candidate_count == 2
    assert len(result.rows) == 2
    assert result.constraint_counts["eq_total"] > 0
    assert result.constraint_counts["ub_total"] == 1
    assert result.baseline_objective_value is not None
    assert result.result_status == "draft"
    assert result.matlab_alignment_status == "pending"
    assert summary["candidate_count"] == 2

    for row in result.rows:
        assert REQUIRED_ROW_FIELDS.issubset(row)
        assert row["objective_value"] is not None
        assert row["baseline_objective_value"] == pytest.approx(result.baseline_objective_value)


@slow_screen
@pytest.mark.parametrize("target_id", ("OPN_ALPHA_FULL_PROJECT", "hLF"))
def test_builtin_targets_run_knockout_screen_smoke(target_id: str) -> None:
    inputs = _inputs()
    target = _builtin_targets()[target_id]

    result = run_knockout_screen(
        inputs.prepared_model,
        target,
        inputs.amino_acids,
        inputs.metabolic,
        inputs.secretory,
        inputs.combined,
        genes=list(KO_GENES),
    )

    _assert_common_screen_result(result, target_id, "knockout")
    for row in result.rows:
        assert row["intervention_type"] == "KO"
        assert row["gene_id"] in KO_GENES
        assert row["input_gene_id"] == row["gene_id"]
        assert row["resolved_reaction_id"] == row["reaction_id"]
        assert row["effect_label"] in {"提升分泌", "降低分泌", "无明显变化", "求解失败", "未解析"}
        assert row["secretory_process"]
        assert row["reaction_id"]


@slow_screen
@pytest.mark.parametrize("target_id", ("OPN_ALPHA_FULL_PROJECT", "hLF"))
def test_builtin_targets_run_overexpression_screen_smoke(target_id: str) -> None:
    inputs = _inputs()
    target = _builtin_targets()[target_id]

    result = run_overexpression_screen(
        inputs.prepared_model,
        target,
        inputs.amino_acids,
        inputs.metabolic,
        inputs.secretory,
        inputs.combined,
        reactions=list(OE_REACTIONS),
    )

    _assert_common_screen_result(result, target_id, "overexpression")
    for row in result.rows:
        assert row["intervention_type"] == "OE_reaction"
        assert row["gene_id"] is None
        assert row["input_gene_id"] is None
        assert row["resolved_reaction_id"] == row["reaction_id"]
        assert row["effect_label"] in {"提升分泌", "降低分泌", "无明显变化", "求解失败", "未解析"}
        assert row["secretory_process"]
        assert row["reaction_id"] in OE_REACTIONS
        assert row["complex_subunit_ids"]
        assert row["complex_subunit_stoichiometry"]


@slow_screen
def test_builtin_opn_runs_reaction_knockout_screen_smoke() -> None:
    inputs = _inputs()
    target = _builtin_targets()["OPN_ALPHA_FULL_PROJECT"]

    result = run_reaction_knockout_screen(
        inputs.prepared_model,
        target,
        inputs.amino_acids,
        inputs.metabolic,
        inputs.secretory,
        inputs.combined,
        reactions=list(KO_REACTIONS),
        write_ribosome_translation_constraint=True,
        write_misfolding_constraints=True,
    )

    summary = summarize_screen_result(result)

    assert isinstance(result, ScreenResult)
    assert result.target_id == "OPN_ALPHA_FULL_PROJECT"
    assert result.screen_type == "knockout"
    assert result.candidate_count == 1
    assert len(result.rows) == 1
    assert summary["candidate_count"] == 1
    assert result.constraint_counts["ribosome_translation"] == 1
    assert result.constraint_counts["misfolding"] == 1418
    row = result.rows[0]
    assert REQUIRED_ROW_FIELDS.issubset(row)
    assert row["intervention_type"] == "KO_reaction"
    assert row["gene_id"] is None
    assert row["reaction_id"] == KO_REACTIONS[0]
    assert row["resolved_reaction_id"] == KO_REACTIONS[0]
    assert row["baseline_objective_value"] == pytest.approx(result.baseline_objective_value)
    assert row["complex_subunit_ids"]
    assert row["status"] == "2"
    assert row["success"] is False
    assert row["objective_value"] is None
    assert row["delta_objective"] is None
    assert row["effect_label"] == "约束不可行"
    assert row["solver_status_label"] == "约束不可行"
    assert row["failure_reason"] == "约束不可行"
