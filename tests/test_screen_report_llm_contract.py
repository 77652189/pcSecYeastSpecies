from __future__ import annotations

import csv
import json
from pathlib import Path

from pcsec_pichia.core.paths import ProjectPaths

from app.services.screen_report_service import (
    build_fact_pack_for_runs,
    generate_judged_screen_report,
    latest_report_runs,
    render_screen_report_markdown,
)


class FakeWriter:
    def __init__(self, drafts: list[dict[str, object]]) -> None:
        self.drafts = drafts
        self.calls: list[dict[str, object]] = []

    def complete_json(self, *, system_prompt: str, payload: dict[str, object]) -> dict[str, object]:
        self.calls.append(payload)
        return self.drafts[min(len(self.calls) - 1, len(self.drafts) - 1)]


class FakeJudge:
    def __init__(self, verdicts: list[str]) -> None:
        self.verdicts = verdicts
        self.calls: list[dict[str, object]] = []

    def complete_json(self, *, system_prompt: str, payload: dict[str, object]) -> dict[str, object]:
        verdict = self.verdicts[min(len(self.calls), len(self.verdicts) - 1)]
        self.calls.append(payload)
        if verdict == "pass":
            return {"verdict": "pass", "blocking_issues": [], "required_fixes": []}
        return {"verdict": "fail", "blocking_issues": [{"type": "readability", "message": "revise"}], "required_fixes": ["revise"]}


def _paths(tmp_path: Path) -> ProjectPaths:
    for directory in ("local_runs", "Results", "Data", "Model", "Enzymedata"):
        (tmp_path / directory).mkdir(parents=True, exist_ok=True)
    return ProjectPaths(tmp_path)


def _csv(tmp_path: Path, run_name: str = "fixture_run") -> Path:
    run_dir = tmp_path / "local_runs" / run_name
    run_dir.mkdir(parents=True)
    path = run_dir / "gene_tradeoff_rows.csv"
    path.write_text(
        "\n".join(
            [
                "target_id,gene_id,candidate_kind,intervention_type,support_status,secretion_ratio_vs_wildtype,growth_retention_ratio,max_feasible_mu",
                "hLF,PAS_HLF_KO,gene,KO,ko_runnable_gpr_gene_deletion,1.2,1.0,0.1",
            ]
        ),
        encoding="utf-8",
    )
    return path


def _valid_report() -> dict[str, object]:
    return {
        "schema_version": 1,
        "targets": {
            "hLF": {
                "executive_summary": "hLF summary",
                "recommended_ko": [
                    {"evidence_id": "hLF-KO-0001", "claim": "PAS_HLF_KO 可优先复核", "rationale": "分泌相对提升", "risk": "模型内结论", "next_step": "小规模验证"}
                ],
                "recommended_oe": [],
                "manual_review": [],
                "not_recommended_or_risky": [],
                "evidence_boundaries": ["不做真实产量承诺"],
            },
            "OPN": {"executive_summary": "", "recommended_ko": [], "recommended_oe": [], "manual_review": [], "not_recommended_or_risky": [], "evidence_boundaries": []},
        },
        "global_warnings": ["不预测真实发酵结果。"],
    }


def test_judge_fail_feedback_reaches_writer_then_pass_generates_final_report(tmp_path: Path) -> None:
    writer = FakeWriter([_valid_report(), _valid_report()])
    judge = FakeJudge(["fail", "pass"])

    result = generate_judged_screen_report(
        _paths(tmp_path),
        csv_paths=(_csv(tmp_path),),
        writer_client=writer,
        judge_client=judge,
        max_rounds=3,
    )

    assert result.success is True
    assert result.final_report_path is not None
    assert result.final_report_path.exists()
    assert len(writer.calls) == 2
    assert writer.calls[1]["feedback_to_fix"]
    assert "evidence_items" not in writer.calls[-1]["fact_pack"]
    assert judge.calls[-1]["fact_pack"]["cited_evidence_items"][0]["gene_id"] == "PAS_HLF_KO"
    assert "fact_pack_summary" in judge.calls[-1]


def test_generation_does_not_write_final_report_when_validator_never_passes(tmp_path: Path) -> None:
    writer = FakeWriter([_valid_report() | {"targets": {"hLF": {"recommended_ko": [{"evidence_id": "hLF-KO-9999", "claim": "PAS_FAKE", "rationale": "x", "risk": "x", "next_step": "x"}]}}}])
    judge = FakeJudge(["pass"])

    result = generate_judged_screen_report(
        _paths(tmp_path),
        csv_paths=(_csv(tmp_path),),
        writer_client=writer,
        judge_client=judge,
        max_rounds=1,
    )

    assert result.success is False
    assert result.final_report_path is None
    assert not (result.output_dir / "final_report.md").exists()


def test_latest_report_runs_can_filter_by_source_run(tmp_path: Path) -> None:
    paths = _paths(tmp_path)
    writer = FakeWriter([_valid_report()])
    judge = FakeJudge(["pass"])

    first = generate_judged_screen_report(
        paths,
        csv_paths=(_csv(tmp_path, "first_run"),),
        writer_client=writer,
        judge_client=judge,
        max_rounds=1,
    )
    second = generate_judged_screen_report(
        paths,
        csv_paths=(_csv(tmp_path, "second_run"),),
        writer_client=FakeWriter([_valid_report()]),
        judge_client=FakeJudge(["pass"]),
        max_rounds=1,
    )

    assert first.output_dir != second.output_dir
    first_runs = latest_report_runs(paths, run_name="first_run")
    second_runs = latest_report_runs(paths, run_name="second_run")
    assert first_runs == [first.output_dir]
    assert second_runs == [second.output_dir]


def test_fact_pack_preserves_external_evidence_fields_without_tier_upgrade(tmp_path: Path) -> None:
    paths = _paths(tmp_path)
    csv_path = _external_evidence_csv(tmp_path)

    fact_pack = build_fact_pack_for_runs(paths, csv_paths=(csv_path,))

    item = fact_pack["evidence_items"][0]
    assert item["recommendation_tier"] == "model_executable"
    assert item["database_annotation_sources"] == ["UniProt", "KEGG"]
    assert item["database_annotation_confidence"] == "reviewed_structured_annotation"
    assert item["model_gpr_executable"] is True
    assert item["oe_reaction_proxy"] is False
    assert item["phenotype_evidence"] == {}
    assert item["external_gene_function_confidence"] == ["reviewed_structured_annotation"]
    assert item["external_gene_function_evidence"][0]["function_description"] == "Annotation-only secretion function"
    assert item["external_gpr_candidate_evidence"][0]["gpr_transfer_status"] == "gene_mapping_required"
    assert item["ko_oe_external_gene_evidence"]["external_name_status"] == "external_match_confirmed"
    assert item["recommendation_tier"] != "experiment_calibrated"

    brief = fact_pack["targets"]["hLF"]["top_candidates"][0]
    assert brief["external_gene_function_sources"] == ["UniProt"]
    assert brief["external_gpr_candidate_evidence"][0]["gpr_transfer_status"] == "gene_mapping_required"


def test_markdown_reference_table_is_program_generated() -> None:
    fact_pack = {
        "evidence_items": [
            {"evidence_id": "hLF-KO-0001", "target_id": "hLF", "gene_id": "PAS_HLF_KO", "intervention_type": "KO", "source_run": "fixture", "secretion_ratio_vs_wildtype": 1.2}
        ],
        "warnings": [],
    }
    markdown = render_screen_report_markdown(fact_pack, _valid_report())

    assert "## hLF 总结" in markdown
    assert "## OPN 总结" in markdown
    assert "| evidence_id | target_id | gene_id/reaction_id |" in markdown
    assert "hLF-KO-0001" in markdown


def _external_evidence_csv(tmp_path: Path) -> Path:
    run_dir = tmp_path / "local_runs" / "external_fixture_run"
    run_dir.mkdir(parents=True)
    path = run_dir / "gene_tradeoff_rows.csv"
    row = {
        "target_id": "hLF",
        "gene_id": "PAS_EXTERNAL_KO",
        "candidate_kind": "gene",
        "intervention_type": "KO",
        "support_status": "ko_runnable_gpr_gene_deletion",
        "recommendation_tier": "model_executable",
        "database_annotation_sources": json.dumps(["UniProt", "KEGG"]),
        "database_annotation_confidence": "reviewed_structured_annotation",
        "model_gpr_executable": "true",
        "oe_reaction_proxy": "false",
        "phenotype_evidence": json.dumps({}),
        "external_gene_function_sources": json.dumps(["UniProt"]),
        "external_gene_function_confidence": json.dumps(["reviewed_structured_annotation"]),
        "external_gene_function_evidence": json.dumps(
            [
                {
                    "source_database": "uniprot",
                    "function_description": "Annotation-only secretion function",
                    "evidence_scope": "reviewed_structured_annotation",
                }
            ]
        ),
        "external_gpr_candidate_evidence": json.dumps(
            [
                {
                    "source_database": "yeast-gem",
                    "gpr_transfer_status": "gene_mapping_required",
                    "blocking_reasons": ["external gene rule is not mapped"],
                }
            ]
        ),
        "ko_oe_external_gene_evidence": json.dumps(
            {
                "pichia_gene_id": "PAS_EXTERNAL_KO",
                "external_name_status": "external_match_confirmed",
                "manual_review_reasons": [],
            }
        ),
        "secretion_ratio_vs_wildtype": "1.2",
        "growth_retention_ratio": "1.0",
        "max_feasible_mu": "0.1",
    }
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row))
        writer.writeheader()
        writer.writerow(row)
    return path
