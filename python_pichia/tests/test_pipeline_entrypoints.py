from __future__ import annotations

import ast
import csv
import json
import os
from pathlib import Path

import pytest

from pcsec_pichia.engines.base import PichiaSimulationRequest, PichiaSimulationRunResult
from pcsec_pichia.pipeline import (
    _alignment_request_for_target,
    _build_screen_plan,
    _build_pipeline_report,
    _target_metadata,
    run_pichia_secretion_simulation,
)
from pcsec_pichia.screens import resolve_oe_gene_reactions, split_existing_genes
from pcsec_pichia.targets import target_spec_from_mapping


REPO_ROOT = Path(__file__).resolve().parents[2]
CUSTOM_TARGETS = REPO_ROOT / "local_runs" / "pichia_hlf_opn_probe" / "targets.example.json"
slow_pipeline = pytest.mark.skipif(
    os.environ.get("PCSEC_RUN_SLOW_PIPELINE_TESTS") != "1",
    reason="slow pcSec pipeline solve; set PCSEC_RUN_SLOW_PIPELINE_TESTS=1 to run",
)


def test_pipeline_solve_tests_are_slow_gated() -> None:
    module_ast = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    ungated: list[str] = []
    for node in module_ast.body:
        if not isinstance(node, ast.FunctionDef) or not node.name.startswith("test_"):
            continue
        calls_pipeline_solve = any(
            isinstance(child, ast.Call)
            and isinstance(child.func, ast.Name)
            and child.func.id == "run_pichia_secretion_simulation"
            for child in ast.walk(node)
        )
        has_slow_marker = any(
            isinstance(decorator, ast.Name) and decorator.id == "slow_pipeline"
            for decorator in node.decorator_list
        )
        if calls_pipeline_solve and not has_slow_marker:
            ungated.append(node.name)
    assert ungated == []


def _assert_common_pipeline_outputs(result: PichiaSimulationRunResult, output_dir: Path) -> dict[str, object]:
    assert isinstance(result, PichiaSimulationRunResult)
    assert result.success is True
    assert result.summary_path is not None
    assert result.report_path is not None
    assert result.candidate_table_path is not None
    assert result.tradeoff_path is not None
    for path in (result.summary_path, result.report_path, result.candidate_table_path, result.tradeoff_path):
        assert path.exists()
        assert path.parent == output_dir
        assert path.stat().st_size > 0

    summary = json.loads(result.summary_path.read_text(encoding="utf-8"))
    assert summary["target_id"] == result.target_id
    assert summary["success"] is True
    assert summary["objective_value"] is not None
    assert summary["constraint_counts"] == result.constraint_counts
    assert summary["result_status"] == "draft"
    assert summary["matlab_alignment_status"] == "pending"
    assert summary["alignment_summary"] == result.alignment_summary
    assert summary["compatibility_mode"] == "corrected"
    assert result.result_status == "corrected_condition"
    assert result.matlab_alignment_status != "aligned"

    with result.candidate_table_path.open(newline="", encoding="utf-8") as handle:
        candidate_rows = list(csv.DictReader(handle))
    assert candidate_rows
    assert {"gene_id", "reaction_id", "input_gene_id", "resolved_reaction_id", "intervention_type", "effect_label", "solver_status_label", "failure_reason", "secretory_process", "objective_value", "baseline_objective_value", "delta_objective", "complex_subunit_ids", "complex_subunit_stoichiometry"}.issubset(candidate_rows[0])

    with result.tradeoff_path.open(newline="", encoding="utf-8") as handle:
        tradeoff_rows = list(csv.DictReader(handle))
    assert len(tradeoff_rows) == 1
    return summary


def test_oe_gene_resolution_expands_model_rules() -> None:
    class TinyModel:
        rxns = ["R1", "R2", "R3"]
        rules = ["x(1)", "x(1) | x(2)", "[]"]
        gr_rules = ["", "G1 or G2", ""]
        gene_index = {"G1": 0, "G2": 1}
        reaction_index = {"R1": 0, "R2": 1, "R3": 2}

    reactions, mapping, unresolved, warnings = resolve_oe_gene_reactions(TinyModel(), ("G1", "G2", "NO_SUCH_GENE"), 10)

    assert reactions == ["R1", "R2"]
    assert mapping == {"R1": "G1", "R2": "G1,G2"}
    assert unresolved == ("NO_SUCH_GENE",)
    assert warnings == []


def test_gene_resolution_splits_existing_and_unresolved_ids() -> None:
    class TinyModel:
        gene_index = {"G1": 0, "G2": 1}

    existing, unresolved = split_existing_genes(TinyModel(), ("G1", "NO_SUCH_GENE", "G2"))

    assert existing == ["G1", "G2"]
    assert unresolved == ("NO_SUCH_GENE",)


def test_screen_plan_uses_manual_candidates_without_solving() -> None:
    class TinyModel:
        rxns = ["R1", "R2", "R3"]
        rules = ["x(1)", "x(1) | x(2)", "[]"]
        gr_rules = ["G1", "G1 or G2", ""]
        gene_index = {"G1": 0, "G2": 1}
        reaction_index = {"R1": 0, "R2": 1, "R3": 2}

    request = PichiaSimulationRequest(
        target_id="OPN_ALPHA_FULL_PROJECT",
        candidate_id="OPN_ALPHA_FULL_PROJECT",
        ko_gene_ids=("G1", "NO_SUCH_KO_GENE"),
        ko_reaction_ids=("R1", "NO_SUCH_KO_RXN"),
        oe_gene_ids=("G1", "NO_SUCH_OE_GENE"),
        oe_reaction_ids=("R2", "NO_SUCH_OE_RXN"),
        screen_candidate_limit=2,
    )

    plan = _build_screen_plan(TinyModel(), request)

    assert plan["ko_gene_ids"] == ["G1"]
    assert plan["unresolved_ko_gene_ids"] == ("NO_SUCH_KO_GENE",)
    assert plan["ko_reaction_ids"] == ["R1"]
    assert plan["unresolved_ko_reaction_ids"] == ("NO_SUCH_KO_RXN",)
    assert plan["oe_gene_reaction_ids"] == ["R1", "R2"]
    assert plan["oe_gene_by_reaction"] == {"R1": "G1", "R2": "G1"}
    assert plan["unresolved_oe_gene_ids"] == ("NO_SUCH_OE_GENE",)
    assert plan["oe_reaction_ids"] == ["R2"]
    assert plan["unresolved_oe_reaction_ids"] == ("NO_SUCH_OE_RXN",)
    assert plan["candidate_limit"] == 2
    assert plan["warnings"]


def test_alignment_request_maps_project_hlf_to_harness_artifact() -> None:
    request = _alignment_request_for_target("hLF", "corrected", REPO_ROOT)

    assert request["target_id"] == "hLF_PROJECT_710"
    assert request["artifact_path"].name == "hlf_project_sequence_matlab_harness_summary.json"
    assert request["constraint_diff_status"] == "known_matlab_compatibility_differences"
    assert {item["id"] for item in request["compatibility_exceptions"]} == {
        "corrected_medium_exchange_bounds",
        "misfolding_dilution_bounds",
        "ribosome_optional_row_mapping",
    }


@slow_pipeline
def test_pipeline_runs_builtin_opn_with_optional_constraints(tmp_path: Path) -> None:
    output_dir = tmp_path / "opn"
    result = run_pichia_secretion_simulation(
        PichiaSimulationRequest(
            target_id="OPN_ALPHA_FULL_PROJECT",
            candidate_id="OPN_ALPHA_FULL_PROJECT",
            enable_ribosome_translation_constraint=True,
            enable_misfolding_constraint=True,
        ),
        output_dir=output_dir,
    )

    summary = _assert_common_pipeline_outputs(result, output_dir)
    assert result.target_id == "OPN_ALPHA_FULL_PROJECT"
    assert result.constraint_counts["ribosome_translation"] == 1
    assert result.constraint_counts["misfolding"] == 1418
    assert result.alignment_summary["matlab_alignment_status"] == "aligned_except_known_matlab_compatibility_differences"
    assert result.alignment_summary["is_fully_aligned"] is False
    assert len(result.alignment_summary["compatibility_exceptions"]) == 4
    assert summary["secretion_plan"]["route_kind"] == "opn_like_soluble_secretory"
    report = result.report_path.read_text(encoding="utf-8")
    assert "Known MATLAB compatibility exceptions" in report or "已知 MATLAB 兼容差异" in report
    assert "corrected_condition" in report


@slow_pipeline
def test_pipeline_runs_builtin_hlf_with_project_710_alignment_artifact(tmp_path: Path) -> None:
    output_dir = tmp_path / "hlf"
    result = run_pichia_secretion_simulation(
        PichiaSimulationRequest(target_id="hLF", candidate_id="hLF"),
        output_dir=output_dir,
    )

    summary = _assert_common_pipeline_outputs(result, output_dir)
    assert result.target_id == "hLF"
    assert result.alignment_summary["target_id"] == "hLF_PROJECT_710"
    assert result.alignment_summary["python_target_id"] == "hLF"
    assert result.alignment_summary["alignment_artifact_target_id"] == "hLF_PROJECT_710"
    assert result.alignment_summary["matlab_alignment_status"] == "aligned_except_known_matlab_compatibility_differences"
    assert result.alignment_summary["is_fully_aligned"] is False
    assert summary["target_parameter_status"] == "draft_matlab_alignment_pending"
    assert summary["target_metadata"]["alignment_target_kind"] == "project_defined_hLF"
    assert summary["target_metadata"]["sequence_role"] == "native_signal_plus_mature_hLF"
    assert summary["target_metadata"]["normalization_mode"] == "user_provided_as_provided"
    assert summary["target_metadata"]["full_sequence_length"] == 710
    assert summary["target_metadata"]["mature_sequence_length"] == 691
    assert summary["target_metadata"]["leader_sequence_length"] == 19
    assert summary["target_metadata"]["signal_peptide_length"] == 19
    assert any("用户提供" in item for item in summary["target_warnings"])
    report = result.report_path.read_text(encoding="utf-8")
    assert "project_defined_hLF" in report
    assert "hLF_PROJECT_710" in report


def test_pipeline_metadata_records_custom_sequence_contract_without_solving() -> None:
    target = target_spec_from_mapping(
        {
            "target_id": "CUSTOM_STOP_STRIPPED",
            "protein_id": "CUSTOM_STOP_STRIPPED",
            "mature_sequence": "ACD",
            "leader_sequence": "MMAA",
            "signal_peptide_sequence": "MM",
            "through_er": True,
            "localization": "e",
        },
        source="request.target_input",
    )
    request = PichiaSimulationRequest(
        target_id="CUSTOM_STOP_STRIPPED",
        candidate_id="CUSTOM_STOP_STRIPPED",
        target_input={"target_id": "CUSTOM_STOP_STRIPPED"},
        sequence_role="mature_secreted",
        normalization_mode="remove_terminal_stop",
        contains_signal_peptide=True,
        contains_leader=True,
        terminal_stop_policy="strip",
        original_sequence_length=4,
        normalized_sequence_length=3,
        original_full_sequence_length=8,
        normalized_full_sequence_length=7,
        original_leader_sequence_length=4,
        normalized_leader_sequence_length=4,
        original_signal_peptide_length=2,
        normalized_signal_peptide_length=2,
        terminal_stop_present=True,
        terminal_stop_removed=True,
    )

    metadata = _target_metadata(target, request)

    assert metadata["sequence_role"] == "mature_secreted"
    assert metadata["normalization_mode"] == "remove_terminal_stop"
    assert metadata["contains_signal_peptide"] is True
    assert metadata["contains_leader"] is True
    assert metadata["terminal_stop_policy"] == "strip"
    assert metadata["original_sequence_length"] == 4
    assert metadata["normalized_sequence_length"] == 3
    assert metadata["original_full_sequence_length"] == 8
    assert metadata["normalized_full_sequence_length"] == 7
    assert metadata["terminal_stop_present"] is True
    assert metadata["terminal_stop_removed"] is True


def test_pipeline_report_includes_draft_and_candidate_interpretation() -> None:
    report = _build_pipeline_report(
        {
            "target_id": "OPN_ALPHA_FULL_PROJECT",
            "result_status": "draft",
            "matlab_alignment_status": "pending",
            "target_parameter_status": "draft",
            "success": True,
            "objective_value": 0.1,
            "growth_rate": 0.1,
            "constraint_counts": {"misfolding": 1418},
            "candidate_count": 1,
            "tradeoff": {"tradeoff_rows": [{"mu": 0.1}]},
            "candidate_interpretation": {
                "effect_counts": {"提升分泌": 1},
                "top_explanations": [
                    {
                        "summary": "OE_reaction `sec_TEST_complex_formation`：提升分泌；关联环节：ER；Δobjective=1e-06。",
                    }
                ],
            },
        }
    )

    assert "Python 草稿结果" in report
    assert "MATLAB 对齐状态" in report
    assert "## 候选解释" in report
    assert "提升分泌" in report


def test_pipeline_report_keeps_hlf_project_710_alignment_semantics() -> None:
    report = _build_pipeline_report(
        {
            "target_id": "hLF",
            "result_status": "corrected_condition",
            "matlab_alignment_status": "aligned_except_known_matlab_compatibility_differences",
            "target_parameter_status": "draft_matlab_alignment_pending",
            "success": True,
            "objective_value": 0.1,
            "growth_rate": 0.1,
            "constraint_counts": {},
            "candidate_count": 0,
            "tradeoff": {"tradeoff_rows": []},
            "alignment_summary": {
                "target_id": "hLF_PROJECT_710",
                "python_target_id": "hLF",
                "alignment_artifact_target_id": "hLF_PROJECT_710",
                "matlab_alignment_status": "aligned_except_known_matlab_compatibility_differences",
                "is_fully_aligned": False,
                "compatibility_exceptions": [
                    {"id": "corrected_medium_exchange_bounds", "count": 9, "category": "bound_difference"}
                ],
            },
            "target_metadata": {
                "alignment_target_kind": "project_defined_hLF",
                "sequence_role": "native_signal_plus_mature_hLF",
                "normalization_mode": "user_provided_as_provided",
                "full_sequence_length": 710,
                "mature_sequence_length": 691,
                "leader_sequence_length": 19,
                "signal_peptide_length": 19,
                "disulfide_sites": 21,
                "n_glycosylation_sites": 4,
                "o_glycosylation_sites": 0,
            },
            "target_warnings": [
                "hLF 使用用户提供的 710aa 目标序列：人源天然信号肽 19aa + mature hLF 691aa。",
                "当前 Python target `hLF` 对应 MATLAB artifact target `hLF_PROJECT_710`，可报告为 `aligned_except_known_matlab_compatibility_differences`，但不是 fully aligned。",
                "旧 MATLAB hLF baseline 仍保持 matlab_failed；Python corrected 结果不能声明为旧 MATLAB fully aligned。",
            ],
        }
    )

    assert "MATLAB artifact target: `hLF_PROJECT_710`" in report
    assert "Python target `hLF` 对应 MATLAB artifact target `hLF_PROJECT_710`" in report
    assert "aligned_except_known_matlab_compatibility_differences" in report
    assert "fully aligned" in report
    assert "matlab_failed" in report


@slow_pipeline
def test_pipeline_uses_manual_ko_oe_candidates(tmp_path: Path) -> None:
    output_dir = tmp_path / "manual_candidates"
    result = run_pichia_secretion_simulation(
        PichiaSimulationRequest(
            target_id="OPN_ALPHA_FULL_PROJECT",
            candidate_id="OPN_ALPHA_FULL_PROJECT",
            ko_gene_ids=("PAS_chr2-2_0107", "NO_SUCH_KO_GENE"),
            ko_reaction_ids=("sec_Och1p_complex_formation", "NO_SUCH_KO_REACTION"),
            oe_gene_ids=("NO_SUCH_OE_GENE",),
            oe_reaction_ids=("sec_BIP_NEFS_complex_formation", "NO_SUCH_OE_REACTION"),
            screen_candidate_limit=2,
        ),
        output_dir=output_dir,
    )

    summary = _assert_common_pipeline_outputs(result, output_dir)
    with result.candidate_table_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    assert summary["screen_request"]["ko_gene_ids"] == ["PAS_chr2-2_0107", "NO_SUCH_KO_GENE"]
    assert summary["screen_request"]["ko_reaction_ids"] == ["sec_Och1p_complex_formation", "NO_SUCH_KO_REACTION"]
    assert summary["screen_request"]["oe_gene_ids"] == ["NO_SUCH_OE_GENE"]
    assert summary["screen_request"]["oe_reaction_ids"] == ["sec_BIP_NEFS_complex_formation", "NO_SUCH_OE_REACTION"]
    assert any(row["intervention_type"] == "KO" and row["input_gene_id"] == "PAS_chr2-2_0107" for row in rows)
    assert any(row["intervention_type"] == "KO" and row["status"] == "unresolved_gene" for row in rows)
    assert any(
        row["intervention_type"] == "KO_reaction"
        and row["reaction_id"] == "sec_Och1p_complex_formation"
        and row["effect_label"] in {"约束不可行", "求解失败"}
        for row in rows
    )
    assert any(row["intervention_type"] == "KO_reaction" and row["status"] == "unresolved_reaction" for row in rows)
    assert any(row["intervention_type"] == "OE_gene_proxy" and row["status"] == "unresolved_gene" for row in rows)
    assert any(row["intervention_type"] == "OE_reaction" and row["reaction_id"] == "sec_BIP_NEFS_complex_formation" for row in rows)
    assert any(row["intervention_type"] == "OE_reaction" and row["status"] == "unresolved_reaction" for row in rows)
    assert any("过表达基因" in item for item in summary["screen_warnings"])


@slow_pipeline
def test_pipeline_runs_custom_json_target(tmp_path: Path) -> None:
    output_dir = tmp_path / "custom"
    result = run_pichia_secretion_simulation(
        PichiaSimulationRequest(
            target_id="OPN_CUSTOM",
            candidate_id="OPN_CUSTOM",
            target_input=CUSTOM_TARGETS,
        ),
        output_dir=output_dir,
    )

    summary = _assert_common_pipeline_outputs(result, output_dir)
    assert result.target_id == "OPN_CUSTOM"
    assert summary["secretion_plan"]["supported"] is True
    assert result.alignment_summary["matlab_alignment_status"] in {"baseline_missing", "pending"}


@slow_pipeline
def test_pipeline_raises_clear_error_for_unknown_target(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="Available targets"):
        run_pichia_secretion_simulation(
            PichiaSimulationRequest(target_id="NO_SUCH_TARGET", candidate_id="NO_SUCH_TARGET"),
            output_dir=tmp_path,
        )
