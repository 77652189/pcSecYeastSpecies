from __future__ import annotations

from pathlib import Path

from pcsec_pichia.core.paths import ProjectPaths

from app.services.screen_report_service import generate_judged_screen_report, render_screen_report_markdown


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
        self.calls = 0

    def complete_json(self, *, system_prompt: str, payload: dict[str, object]) -> dict[str, object]:
        verdict = self.verdicts[min(self.calls, len(self.verdicts) - 1)]
        self.calls += 1
        if verdict == "pass":
            return {"verdict": "pass", "blocking_issues": [], "required_fixes": []}
        return {"verdict": "fail", "blocking_issues": [{"type": "readability", "message": "revise"}], "required_fixes": ["revise"]}


def _paths(tmp_path: Path) -> ProjectPaths:
    for directory in ("local_runs", "Results", "Data", "Model", "Enzymedata"):
        (tmp_path / directory).mkdir(parents=True, exist_ok=True)
    return ProjectPaths(tmp_path)


def _csv(tmp_path: Path) -> Path:
    run_dir = tmp_path / "local_runs" / "fixture_run"
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
