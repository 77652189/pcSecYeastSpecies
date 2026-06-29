from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from pcsec_pichia.reports import (
    CANDIDATE_COLUMNS,
    ReportOutputResult,
    build_markdown_report,
    build_candidate_interpretation,
    build_summary_payload,
    audit_candidate_report_consistency,
    normalize_candidate_explanation_row,
    summarize_report_outputs,
    write_candidate_table,
    write_simulation_outputs,
    write_tradeoff_table,
)
from pcsec_pichia.screens import ScreenResult
from pcsec_pichia.simulation import GrowthTradeoffResult, SecretionSimulationResult


def _simulation_result(target_id: str) -> SecretionSimulationResult:
    return SecretionSimulationResult(
        success=True,
        target_id=target_id,
        objective_value=0.0065 if target_id != "hLF" else 0.0031,
        growth_rate=0.10,
        secretion_flux=0.0065 if target_id != "hLF" else 0.0031,
        status="0",
        message="ok",
        constraint_counts={"eq_total": 23015, "ub_total": 1},
        result_status="draft",
        target_parameter_status="draft" if target_id != "hLF" else "draft_matlab_alignment_pending",
        matlab_alignment_status="pending",
        exchange_reaction_id=f"{target_id} exchange",
        build_status="stoichiometric_target_model_built",
    )


def _tradeoff_result(target_id: str) -> GrowthTradeoffResult:
    return GrowthTradeoffResult(
        target_id=target_id,
        success=True,
        growth_points=(0.05, 0.10),
        tradeoff_rows=(
            {
                "mu": 0.05,
                "success": True,
                "status": "0",
                "secretion_flux": 0.004,
                "secretion_per_biomass": 0.08,
                "message": "ok",
                "constraint_counts": {"eq_total": 23015, "ub_total": 1},
            },
            {
                "mu": 0.10,
                "success": True,
                "status": "0",
                "secretion_flux": 0.0065,
                "secretion_per_biomass": 0.065,
                "message": "ok",
                "constraint_counts": {"eq_total": 23015, "ub_total": 1},
            },
        ),
        result_status="draft",
        target_parameter_status="draft" if target_id != "hLF" else "draft_matlab_alignment_pending",
        matlab_alignment_status="pending",
    )


def _screen_results(target_id: str) -> tuple[ScreenResult, ScreenResult]:
    ko = ScreenResult(
        target_id=target_id,
        screen_type="knockout",
        success=True,
        candidate_count=1,
        rows=(
            {
                "target_id": target_id,
                "screen_type": "knockout",
                "candidate_id": "AT250_GQ_6803479",
                "gene_id": "AT250_GQ_6803479",
                "reaction_id": "RPE_no_1_fwd",
                "input_gene_id": "AT250_GQ_6803479",
                "resolved_reaction_id": "RPE_no_1_fwd",
                "effect_label": "降低分泌",
                "secretory_process": "代谢或其它反应",
                "mapping_level": "metabolic_or_other",
                "mapping_confidence": "low",
                "mapping_interpretation": "该基因关联代谢或其它反应；可能间接影响分泌，解释置信度较低。",
                "complex_id": None,
                "intervention_type": "KO",
                "success": True,
                "status": "0",
                "objective_value": 0.0064,
                "baseline_objective_value": 0.0065,
                "delta_objective": -0.0001,
                "complex_subunit_ids": [],
                "complex_subunit_stoichiometry": [],
            },
        ),
        constraint_counts={"eq_total": 23015, "ub_total": 1},
        baseline_objective_value=0.0065,
        result_status="draft",
        matlab_alignment_status="pending",
    )
    oe = ScreenResult(
        target_id=target_id,
        screen_type="overexpression",
        success=True,
        candidate_count=1,
        rows=(
            {
                "target_id": target_id,
                "screen_type": "overexpression",
                "candidate_id": "sec_BIP_NEFS_complex_formation",
                "gene_id": None,
                "reaction_id": "sec_BIP_NEFS_complex_formation",
                "input_gene_id": None,
                "resolved_reaction_id": "sec_BIP_NEFS_complex_formation",
                "effect_label": "提升分泌",
                "secretory_process": "ER 折叠 / 分子伴侣",
                "mapping_level": "complex_subunit",
                "mapping_confidence": "medium",
                "mapping_interpretation": "该基因关联分泌复合体反应（ER 折叠 / 分子伴侣）；OE gene 应解释为 reaction-level OE proxy。",
                "complex_id": "sec_BIP_NEFS_complex",
                "intervention_type": "OE_reaction",
                "success": True,
                "status": "0",
                "objective_value": 0.0066,
                "baseline_objective_value": 0.0065,
                "delta_objective": 0.0001,
                "complex_subunit_ids": ["Kar2p", "Sil1p"],
                "complex_subunit_stoichiometry": [1.0, 1.0],
            },
        ),
        constraint_counts={"eq_total": 23015, "ub_total": 1},
        baseline_objective_value=0.0065,
        result_status="draft",
        matlab_alignment_status="pending",
    )
    return ko, oe


@pytest.mark.parametrize("target_id", ("OPN_ALPHA_FULL_PROJECT", "hLF"))
def test_write_simulation_outputs_creates_stable_report_bundle(tmp_path: Path, target_id: str) -> None:
    simulation = _simulation_result(target_id)
    tradeoff = _tradeoff_result(target_id)
    screens = _screen_results(target_id)

    result = write_simulation_outputs(
        simulation,
        tradeoff,
        screens,
        output_dir=tmp_path / target_id,
        output_prefix=f"{target_id}_smoke",
    )
    summary = summarize_report_outputs(result)

    assert isinstance(result, ReportOutputResult)
    assert result.output_dir == tmp_path / target_id
    assert len(result.written_files) == 4
    for path in result.written_files:
        assert path.exists()
        assert path.parent == tmp_path / target_id
        assert path.stat().st_size > 0

    payload = json.loads(result.summary_path.read_text(encoding="utf-8"))
    assert payload["target_id"] == target_id
    assert payload["success"] is True
    assert payload["objective_value"] == simulation.objective_value
    assert payload["constraint_counts"] == simulation.constraint_counts
    assert payload["result_status"] == "draft"
    assert payload["matlab_alignment_status"] == "pending"
    assert payload["candidate_count"] == 2
    assert payload["candidate_interpretation"]["effect_counts"]["提升分泌"] == 1
    assert payload["candidate_interpretation"]["effect_counts"]["降低分泌"] == 1
    if target_id == "hLF":
        assert payload["target_parameter_status"] == "draft_matlab_alignment_pending"

    report = result.report_path.read_text(encoding="utf-8")
    assert "Python 草稿" in report
    assert "候选解释" in report
    assert "Optional constraints / MATLAB LP alignment still require separate alignment" in report or "可选约束" in report
    assert "hLF" in report and "MATLAB" in report

    with result.candidate_table_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        candidate_rows = list(reader)
    assert len(candidate_rows) == 2
    assert set(CANDIDATE_COLUMNS).issubset(reader.fieldnames or [])
    for field in (
        "gene_id",
        "reaction_id",
        "intervention_type",
        "objective_value",
        "baseline_objective_value",
        "delta_objective",
        "input_gene_id",
        "resolved_reaction_id",
        "effect_label",
        "secretory_process",
        "complex_subunit_ids",
        "complex_subunit_stoichiometry",
        "mapping_level",
        "mapping_confidence",
        "mapping_interpretation",
        "complex_id",
    ):
        assert field in (reader.fieldnames or [])

    with result.tradeoff_path.open(newline="", encoding="utf-8") as handle:
        tradeoff_rows = list(csv.DictReader(handle))
    assert len(tradeoff_rows) == 2
    assert tradeoff_rows[0]["target_id"] == target_id
    assert summary["summary_path"] == str(result.summary_path)


def test_report_helpers_can_write_tables_independently(tmp_path: Path) -> None:
    target_id = "OPN_ALPHA_FULL_PROJECT"
    simulation = _simulation_result(target_id)
    tradeoff = _tradeoff_result(target_id)
    screens = _screen_results(target_id)

    summary = build_summary_payload(simulation, tradeoff, screens)
    markdown = build_markdown_report(summary)
    candidates_path = write_candidate_table(screens, tmp_path / "candidates.csv")
    tradeoff_path = write_tradeoff_table(tradeoff, tmp_path / "tradeoff.csv", target_id=target_id)

    assert summary["target_id"] == target_id
    assert summary["candidate_count"] == 2
    assert "Python 草稿" in markdown
    assert candidates_path.exists()
    assert tradeoff_path.exists()


def test_candidate_interpretation_summarizes_effects_and_infeasible_rows() -> None:
    rows = [
        {
            "candidate_id": "sec_BIP_NEFS_complex_formation",
            "reaction_id": "sec_BIP_NEFS_complex_formation",
            "intervention_type": "OE_reaction",
            "effect_label": "提升分泌",
            "secretory_process": "ER 折叠 / 分子伴侣",
            "mapping_level": "complex_subunit",
            "mapping_confidence": "medium",
            "mapping_interpretation": "该反应映射到分子伴侣复合体。",
            "success": True,
            "status": "0",
            "delta_objective": 0.0001,
        },
        {
            "candidate_id": "sec_Och1p_complex_formation",
            "reaction_id": "sec_Och1p_complex_formation",
            "intervention_type": "KO_reaction",
            "effect_label": "求解失败",
            "solver_status_label": "约束不可行",
            "secretory_process": "N-糖基化 NG",
            "mapping_level": "complex_subunit",
            "mapping_confidence": "medium",
            "mapping_interpretation": "该反应映射到 N-糖基化复合体。",
            "success": False,
            "status": "2",
            "delta_objective": None,
        },
        {
            "candidate_id": "NO_SUCH_GENE",
            "input_gene_id": "NO_SUCH_GENE",
            "intervention_type": "KO",
            "effect_label": "未解析",
            "secretory_process": "未解析",
            "mapping_level": "unresolved",
            "mapping_confidence": "unresolved",
            "mapping_interpretation": "未解析到可解释的模型反应。",
            "success": False,
            "status": "unresolved_gene",
            "delta_objective": None,
        },
    ]

    interpretation = build_candidate_interpretation(rows)
    infeasible = normalize_candidate_explanation_row(rows[1])

    assert interpretation["effect_counts"]["提升分泌"] == 1
    assert interpretation["effect_counts"]["约束不可行"] == 1
    assert interpretation["effect_counts"]["未解析"] == 1
    assert "固定生长下约束不可行" in infeasible["summary"]
    assert "N-糖基化 NG" in infeasible["summary"]
    assert "映射置信度：中" in infeasible["summary"]


def test_candidate_report_consistency_audit_matches_summary_and_markdown() -> None:
    target_id = "OPN_ALPHA_FULL_PROJECT"
    simulation = _simulation_result(target_id)
    tradeoff = _tradeoff_result(target_id)
    screens = _screen_results(target_id)
    rows = [row for screen in screens for row in screen.rows]
    summary = build_summary_payload(simulation, tradeoff, screens)
    markdown = build_markdown_report(summary)

    audit = audit_candidate_report_consistency(rows, summary, markdown)

    assert audit["passed"] is True
    assert audit["issue_count"] == 0
    assert audit["candidate_count"] == 2
    assert audit["effect_counts"] == {"降低分泌": 1, "提升分泌": 1}
    assert "complex_subunit_ids" in audit["required_columns"]
    assert "complex_subunit_stoichiometry" in audit["required_columns"]


def test_candidate_report_consistency_audit_flags_mismatched_counts() -> None:
    rows = [
        {
            column: None
            for column in CANDIDATE_COLUMNS
        }
    ]
    rows[0].update(
        {
            "candidate_id": "NO_SUCH_GENE",
            "input_gene_id": "NO_SUCH_GENE",
            "intervention_type": "KO",
            "effect_label": "未解析",
            "secretory_process": "未解析",
            "mapping_level": "unresolved",
            "mapping_confidence": "unresolved",
            "mapping_interpretation": "未解析到可解释的模型反应。",
            "success": False,
            "status": "unresolved_gene",
        }
    )
    summary = {
        "candidate_count": 2,
        "candidate_interpretation": {"effect_counts": {"提升分泌": 1}},
    }

    audit = audit_candidate_report_consistency(rows, summary, "no useful report text")

    assert audit["passed"] is False
    assert any("candidate_count mismatch" in issue for issue in audit["issues"])
    assert any("effect_counts mismatch" in issue for issue in audit["issues"])
    assert any("Python draft" in issue or "Python 草稿" in issue for issue in audit["issues"])


def test_markdown_report_renders_protein_cost_analysis_section() -> None:
    summary = build_summary_payload(_simulation_result("OPN_ALPHA_FULL_PROJECT"), None, ())
    summary["protein_cost_analysis"] = {
        "result_status": "draft_explanatory",
        "total_relative_score": 100.0,
        "dominant_cost_categories": ["translation", "o_glycosylation"],
        "cost_items": [
            {
                "category": "translation",
                "label": "翻译负担",
                "basis": "full sequence length",
                "raw_value": 383.0,
                "relative_score": 55.0,
                "interpretation": "sequence length burden",
            }
        ],
        "warnings": ["explanatory score only"],
    }

    markdown = build_markdown_report(summary)

    assert "## 目标蛋白成本分析" in markdown
    assert "draft_explanatory" in markdown
    assert "translation" in markdown
    assert "不代表真实发酵产量" in markdown
