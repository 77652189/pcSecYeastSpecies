from __future__ import annotations

import json
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

from pcsec_pichia.core.paths import ProjectPaths

from app.services.screen_report_fact_pack import build_screen_report_fact_pack, summarize_fact_pack
from app.services.screen_report_llm import (
    JsonLlmClient,
    get_default_screen_report_llm_client,
    judge_screen_report,
    write_screen_report_draft,
)
from app.services.screen_report_validator import validate_screen_report_json


REPORT_RUN_ROOT = Path("screen_report_runs")


@dataclass(frozen=True)
class ScreenReportRunResult:
    success: bool
    output_dir: Path
    fact_pack_path: Path
    draft_report_path: Path | None
    judge_report_path: Path | None
    final_report_path: Path | None
    manifest_path: Path
    validator_result: dict[str, Any]
    judge_result: dict[str, Any] | None
    message: str


def build_fact_pack_for_runs(
    paths: ProjectPaths,
    *,
    run_names: tuple[str, ...] | None = None,
    csv_paths: tuple[Path | str, ...] | None = None,
) -> dict[str, Any]:
    return build_screen_report_fact_pack(paths, run_names=run_names, csv_paths=csv_paths)


def generate_judged_screen_report(
    paths: ProjectPaths,
    *,
    run_names: tuple[str, ...] | None = None,
    csv_paths: tuple[Path | str, ...] | None = None,
    writer_client: JsonLlmClient | None = None,
    judge_client: JsonLlmClient | None = None,
    max_rounds: int = 3,
) -> ScreenReportRunResult:
    output_dir = _new_output_dir(paths)
    output_dir.mkdir(parents=True, exist_ok=True)
    fact_pack = build_fact_pack_for_runs(paths, run_names=run_names, csv_paths=csv_paths)
    fact_pack_path = output_dir / "fact_pack.json"
    _write_json(fact_pack_path, fact_pack)

    writer = writer_client or get_default_screen_report_llm_client()
    judge = judge_client or writer
    feedback: list[dict[str, Any]] = []
    latest_draft: dict[str, Any] | None = None
    latest_validator: dict[str, Any] = {"verdict": "fail", "blocking_issues": [{"type": "not_started"}], "warnings": []}
    latest_judge: dict[str, Any] | None = None
    draft_path: Path | None = None
    judge_path: Path | None = None
    final_path: Path | None = None

    for round_index in range(1, max(1, max_rounds) + 1):
        latest_draft = write_screen_report_draft(writer, fact_pack, feedback=feedback)
        draft_path = output_dir / f"draft_report_round{round_index}.json"
        _write_json(draft_path, latest_draft)
        latest_validator = validate_screen_report_json(fact_pack, latest_draft)
        _write_json(output_dir / f"validator_round{round_index}.json", latest_validator)
        if latest_validator.get("verdict") != "pass":
            feedback = list(latest_validator.get("blocking_issues") or [])
            continue
        latest_judge = judge_screen_report(judge, fact_pack, latest_draft, latest_validator)
        judge_path = output_dir / f"judge_report_round{round_index}.json"
        _write_json(judge_path, latest_judge)
        if latest_judge.get("verdict") == "pass":
            final_path = output_dir / "final_report.md"
            final_path.write_text(render_screen_report_markdown(fact_pack, latest_draft), encoding="utf-8")
            manifest_path = _write_manifest(
                output_dir,
                success=True,
                fact_pack_path=fact_pack_path,
                draft_path=draft_path,
                judge_path=judge_path,
                final_path=final_path,
                fact_pack=fact_pack,
                validator_result=latest_validator,
                judge_result=latest_judge,
            )
            return ScreenReportRunResult(
                success=True,
                output_dir=output_dir,
                fact_pack_path=fact_pack_path,
                draft_report_path=draft_path,
                judge_report_path=judge_path,
                final_report_path=final_path,
                manifest_path=manifest_path,
                validator_result=latest_validator,
                judge_result=latest_judge,
                message="final_report.md generated after validator and judge passed.",
            )
        feedback = list(latest_judge.get("blocking_issues") or latest_judge.get("required_fixes") or [])

    failed_path = output_dir / "failed_draft.json"
    if latest_draft is not None:
        _write_json(failed_path, latest_draft)
    manifest_path = _write_manifest(
        output_dir,
        success=False,
        fact_pack_path=fact_pack_path,
        draft_path=failed_path if latest_draft is not None else None,
        judge_path=judge_path,
        final_path=None,
        fact_pack=fact_pack,
        validator_result=latest_validator,
        judge_result=latest_judge,
    )
    return ScreenReportRunResult(
        success=False,
        output_dir=output_dir,
        fact_pack_path=fact_pack_path,
        draft_report_path=failed_path if latest_draft is not None else None,
        judge_report_path=judge_path,
        final_report_path=None,
        manifest_path=manifest_path,
        validator_result=latest_validator,
        judge_result=latest_judge,
        message="Report generation failed validator or judge; final_report.md was not written.",
    )


def render_screen_report_markdown(fact_pack: dict[str, Any], report_json: dict[str, Any]) -> str:
    lines = ["# KO/OE 筛查结果研发建议报告", ""]
    lines.append("本报告由程序 fact pack 限定事实范围，LLM 只负责组织语言；引用表由程序生成。")
    lines.append("")
    for target_key in ("hLF", "OPN"):
        target_report = ((report_json.get("targets") or {}).get(target_key) or {})
        lines.extend([f"## {target_key} 总结", "", str(target_report.get("executive_summary") or "无总结。"), ""])
        for title, bucket in (
            ("推荐 KO", "recommended_ko"),
            ("推荐 OE proxy", "recommended_oe"),
            ("需人工复核", "manual_review"),
            ("不建议/高风险", "not_recommended_or_risky"),
        ):
            lines.extend([f"### {target_key} {title}", ""])
            rows = target_report.get(bucket) or []
            if not rows:
                lines.extend(["- 暂无。", ""])
                continue
            for row in rows:
                lines.append(f"- `{row.get('evidence_id')}` {row.get('claim', '')}")
                lines.append(f"  - 理由：{row.get('rationale', '')}")
                lines.append(f"  - 风险：{row.get('risk', '')}")
                lines.append(f"  - 下一步：{row.get('next_step', '')}")
            lines.append("")
        boundaries = target_report.get("evidence_boundaries") or []
        if boundaries:
            lines.extend([f"### {target_key} 证据边界", ""])
            lines.extend(f"- {item}" for item in boundaries)
            lines.append("")
    global_warnings = report_json.get("global_warnings") or fact_pack.get("warnings") or []
    lines.extend(["## 全局证据边界", ""])
    if global_warnings:
        lines.extend(f"- {warning}" for warning in global_warnings)
    else:
        lines.append("- 当前报告不包含绝对产量或实验成功率承诺。")
    lines.extend(["", "## 引用表", ""])
    lines.extend(_reference_table_lines(fact_pack))
    lines.append("")
    return "\n".join(lines)


def latest_report_runs(paths: ProjectPaths, *, run_name: str | None = None) -> list[Path]:
    root = paths.local_runs_dir / REPORT_RUN_ROOT
    if not root.exists():
        return []
    runs = sorted((path for path in root.iterdir() if path.is_dir()), key=lambda path: path.stat().st_mtime, reverse=True)
    if run_name is None:
        return runs
    return [path for path in runs if _report_run_matches_source(path, run_name)]


def _reference_table_lines(fact_pack: dict[str, Any]) -> list[str]:
    lines = [
        "| evidence_id | target_id | gene_id/reaction_id | standard_symbol | protein_name | intervention_type | source_run | key metrics | warnings |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for item in fact_pack.get("evidence_items") or []:
        gene_or_reaction = item.get("gene_id") or item.get("reaction_id") or ""
        metrics = (
            f"secretion_ratio={item.get('secretion_ratio_vs_wildtype')}; "
            f"growth_retention={item.get('growth_retention_ratio')}; "
            f"max_mu={item.get('max_feasible_mu')}"
        )
        warnings = "; ".join(str(value) for value in item.get("warnings") or [])
        lines.append(
            "| "
            + " | ".join(
                _md_cell(value)
                for value in (
                    item.get("evidence_id"),
                    item.get("target_id"),
                    gene_or_reaction,
                    item.get("standard_symbol"),
                    item.get("protein_name"),
                    item.get("intervention_type"),
                    item.get("source_run"),
                    metrics,
                    warnings,
                )
            )
            + " |"
        )
    return lines


def _md_cell(value: Any) -> str:
    return str(value or "").replace("|", "\\|").replace("\n", " ")


def _new_output_dir(paths: ProjectPaths) -> Path:
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    return paths.local_runs_dir / REPORT_RUN_ROOT / stamp


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")


def _write_manifest(
    output_dir: Path,
    *,
    success: bool,
    fact_pack_path: Path,
    draft_path: Path | None,
    judge_path: Path | None,
    final_path: Path | None,
    fact_pack: dict[str, Any],
    validator_result: dict[str, Any],
    judge_result: dict[str, Any] | None,
) -> Path:
    manifest_path = output_dir / "report_manifest.json"
    _write_json(
        manifest_path,
        {
            "success": success,
            "fact_pack_path": str(fact_pack_path),
            "draft_report_path": str(draft_path) if draft_path else None,
            "judge_report_path": str(judge_path) if judge_path else None,
            "final_report_path": str(final_path) if final_path else None,
            "source_run_names": _fact_pack_source_run_names(fact_pack),
            "source_files": _fact_pack_source_files(fact_pack),
            "validator_result": validator_result,
            "judge_result": judge_result,
        },
    )
    return manifest_path


def _report_run_matches_source(report_dir: Path, run_name: str) -> bool:
    manifest_path = report_dir / "report_manifest.json"
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (FileNotFoundError, OSError, json.JSONDecodeError, UnicodeDecodeError):
        return False
    source_run_names = {str(value) for value in manifest.get("source_run_names") or []}
    return run_name in source_run_names


def _fact_pack_source_run_names(fact_pack: dict[str, Any]) -> list[str]:
    names = [
        str(source.get("run_name"))
        for source in fact_pack.get("source_runs") or []
        if source.get("run_name")
    ]
    return sorted(dict.fromkeys(names))


def _fact_pack_source_files(fact_pack: dict[str, Any]) -> list[str]:
    files = [
        str(source.get("source_file"))
        for source in fact_pack.get("source_runs") or []
        if source.get("source_file")
    ]
    return sorted(dict.fromkeys(files))


__all__ = [
    "ScreenReportRunResult",
    "build_fact_pack_for_runs",
    "generate_judged_screen_report",
    "latest_report_runs",
    "render_screen_report_markdown",
    "summarize_fact_pack",
]
