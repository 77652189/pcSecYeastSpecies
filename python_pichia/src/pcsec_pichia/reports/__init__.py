from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable

from pcsec_pichia.reports._prototype_adapter import (
    build_candidate_table,
    candidate_table_row,
    classify_candidate_effect,
    classify_secretory_process,
    format_candidate_rows,
    format_tradeoff_rows,
    write_outputs,
)
from pcsec_pichia.screens import ScreenResult
from pcsec_pichia.simulation import GrowthTradeoffResult, SecretionSimulationResult


CANDIDATE_COLUMNS: tuple[str, ...] = (
    "target_id",
    "screen_type",
    "candidate_id",
    "gene_id",
    "reaction_id",
    "input_gene_id",
    "resolved_reaction_id",
    "intervention_type",
    "effect_label",
    "solver_status_label",
    "failure_reason",
    "secretory_process",
    "success",
    "status",
    "objective_value",
    "baseline_objective_value",
    "delta_objective",
    "complex_subunit_ids",
    "complex_subunit_stoichiometry",
)

CORE_CANDIDATE_EXPLANATION_COLUMNS: tuple[str, ...] = (
    "target_id",
    "screen_type",
    "candidate_id",
    "gene_id",
    "reaction_id",
    "input_gene_id",
    "resolved_reaction_id",
    "intervention_type",
    "effect_label",
    "secretory_process",
    "success",
    "status",
    "objective_value",
    "baseline_objective_value",
    "delta_objective",
    "complex_subunit_ids",
    "complex_subunit_stoichiometry",
)

TRADEOFF_COLUMNS: tuple[str, ...] = (
    "target_id",
    "mu",
    "success",
    "status",
    "secretion_flux",
    "secretion_per_biomass",
    "message",
    "constraint_counts",
)


@dataclass(frozen=True)
class ReportOutputResult:
    output_dir: Path
    summary_path: Path
    report_path: Path
    candidate_table_path: Path
    tradeoff_path: Path
    result_status: str
    matlab_alignment_status: str
    written_files: tuple[Path, ...]


def write_simulation_outputs(
    simulation_result: SecretionSimulationResult,
    tradeoff_result: GrowthTradeoffResult | None,
    screen_results: Iterable[ScreenResult],
    output_dir: Path,
    output_prefix: str,
) -> ReportOutputResult:
    output_dir.mkdir(parents=True, exist_ok=True)
    safe_prefix = _safe_prefix(output_prefix)
    summary_path = output_dir / f"{safe_prefix}_summary.json"
    report_path = output_dir / f"{safe_prefix}_REPORT.md"
    candidate_table_path = output_dir / f"{safe_prefix}_candidates.csv"
    tradeoff_path = output_dir / f"{safe_prefix}_tradeoff.csv"

    screens = tuple(screen_results)
    summary = build_summary_payload(simulation_result, tradeoff_result, screens)
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    report_path.write_text(build_markdown_report(summary), encoding="utf-8")
    write_candidate_table(screens, candidate_table_path)
    write_tradeoff_table(tradeoff_result, tradeoff_path, simulation_result.target_id)

    return ReportOutputResult(
        output_dir=output_dir,
        summary_path=summary_path,
        report_path=report_path,
        candidate_table_path=candidate_table_path,
        tradeoff_path=tradeoff_path,
        result_status=simulation_result.result_status,
        matlab_alignment_status=simulation_result.matlab_alignment_status,
        written_files=(summary_path, report_path, candidate_table_path, tradeoff_path),
    )


def write_candidate_table(screen_results: Iterable[ScreenResult], path: Path) -> Path:
    rows = [row for result in screen_results for row in result.rows]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CANDIDATE_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: _csv_value(row.get(key)) for key in CANDIDATE_COLUMNS})
    return path


def write_tradeoff_table(tradeoff_result: GrowthTradeoffResult | None, path: Path, target_id: str | None = None) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=TRADEOFF_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        for row in tuple(tradeoff_result.tradeoff_rows if tradeoff_result else ()):
            record = {"target_id": target_id or getattr(tradeoff_result, "target_id", None), **row}
            writer.writerow({key: _csv_value(record.get(key)) for key in TRADEOFF_COLUMNS})
    return path


def build_summary_payload(
    simulation_result: SecretionSimulationResult,
    tradeoff_result: GrowthTradeoffResult | None,
    screen_results: Iterable[ScreenResult],
) -> dict[str, Any]:
    screens = tuple(screen_results)
    candidate_rows = [row for result in screens for row in result.rows]
    candidate_interpretation = build_candidate_interpretation(candidate_rows)
    payload = {
        "target_id": simulation_result.target_id,
        "success": simulation_result.success,
        "objective_value": simulation_result.objective_value,
        "growth_rate": simulation_result.growth_rate,
        "secretion_flux": simulation_result.secretion_flux,
        "status": simulation_result.status,
        "constraint_counts": simulation_result.constraint_counts,
        "result_status": simulation_result.result_status,
        "target_parameter_status": simulation_result.target_parameter_status,
        "matlab_alignment_status": simulation_result.matlab_alignment_status,
        "candidate_count": len(candidate_rows),
        "candidate_interpretation": candidate_interpretation,
        "screen_results": [asdict(result) for result in screens],
        "tradeoff": asdict(tradeoff_result) if tradeoff_result else None,
        "alignment_notes": [
        ],
    }
    return payload


def build_markdown_report(summary: dict[str, Any]) -> str:
    target_id = summary.get("target_id")
    status = summary.get("target_parameter_status")
    constraint_counts = summary.get("constraint_counts") or {}
    candidate_interpretation = summary.get("candidate_interpretation") or {}
    protein_cost = summary.get("protein_cost_analysis") or {}
    lines = [
            f"# pcSecPichia Python 分泌仿真报告: {target_id}",
            "",
            "## 状态",
            "",
            "- 当前结果是 Python 草稿结果。",
            f"- 结果状态: `{summary.get('result_status')}`.",
            f"- MATLAB 对齐状态: `{summary.get('matlab_alignment_status')}`.",
            f"- 目标参数状态: `{status}`.",
            "- Optional constraints / MATLAB LP alignment still require separate alignment.",
            "- hLF 的 MATLAB 基线目前尚未对齐。",
            "",
            "## 仿真结果",
            "",
            f"- 是否成功: `{summary.get('success')}`.",
            f"- 目标函数值: `{summary.get('objective_value')}`.",
            f"- 生长速率: `{summary.get('growth_rate')}`.",
            f"- 约束计数: `{constraint_counts}`.",
            "",
            "## 输出",
            "",
            f"- 候选行数: `{summary.get('candidate_count')}`.",
            f"- 生长权衡行数: `{len((summary.get('tradeoff') or {}).get('tradeoff_rows') or [])}`.",
            "",
    ]
    if candidate_interpretation:
        lines.extend(
            [
                "## 候选解释",
                "",
                f"- 分类汇总: `{candidate_interpretation.get('effect_counts') or {}}`.",
                "- `约束不可行` 表示在当前固定生长速率和约束组合下，该扰动没有可行解；这不等同于真实生物系统必然不可行。",
                "",
            ]
        )
        for row in candidate_interpretation.get("top_explanations") or []:
            lines.append(f"- {row.get('summary')}")
        lines.append("")
    if protein_cost:
        lines.extend(_protein_cost_markdown_lines(protein_cost))
    return "\n".join(lines)


def _protein_cost_markdown_lines(protein_cost: dict[str, Any]) -> list[str]:
    items = protein_cost.get("cost_items") or []
    dominant = protein_cost.get("dominant_cost_categories") or []
    lines = [
        "## 目标蛋白成本分析",
        "",
        "- 当前结果是 Python draft explanatory score，不代表真实发酵产量或湿实验成本。",
        f"- 成本分析状态: `{protein_cost.get('result_status')}`.",
        f"- 总相对成本分: `{protein_cost.get('total_relative_score')}`.",
        f"- 主要成本类别: `{', '.join(str(item) for item in dominant)}`.",
        "",
        "| 类别 | 成本项 | 相对分 | 依据 |",
        "| --- | --- | ---: | --- |",
    ]
    for item in items:
        lines.append(
            f"| `{item.get('category')}` | {item.get('label')} | "
            f"{item.get('relative_score')} | {item.get('basis')} |"
        )
    warnings = protein_cost.get("warnings") or []
    if warnings:
        lines.extend(["", "提示:"])
        lines.extend(f"- {warning}" for warning in warnings)
    lines.append("")
    return lines


def build_candidate_interpretation(candidate_rows: Iterable[dict[str, Any]], limit: int = 5) -> dict[str, Any]:
    rows = [normalize_candidate_explanation_row(row) for row in candidate_rows]
    counts: dict[str, int] = {}
    for row in rows:
        bucket = str(row["effect_bucket"])
        counts[bucket] = counts.get(bucket, 0) + 1
    top_rows = sorted(rows, key=_candidate_explanation_sort_key)[: max(0, int(limit))]
    return {
        "effect_counts": counts,
        "top_explanations": top_rows,
    }


def audit_candidate_report_consistency(
    candidate_rows: Iterable[dict[str, Any]],
    summary: dict[str, Any],
    markdown_report: str = "",
) -> dict[str, Any]:
    rows = list(candidate_rows)
    issues: list[str] = []
    warnings: list[str] = []

    missing_core_columns = [
        column
        for column in CORE_CANDIDATE_EXPLANATION_COLUMNS
        if any(column not in row for row in rows)
    ]
    if missing_core_columns:
        issues.append(f"candidate rows missing core columns: {', '.join(missing_core_columns)}")
    missing_display_columns = [
        column
        for column in ("solver_status_label", "failure_reason")
        if any(column not in row for row in rows)
    ]
    if missing_display_columns:
        warnings.append(f"candidate rows missing display columns: {', '.join(missing_display_columns)}")

    expected_count = int(summary.get("candidate_count") or 0)
    if expected_count != len(rows):
        issues.append(f"candidate_count mismatch: summary={expected_count}, rows={len(rows)}")

    row_counts = build_candidate_interpretation(rows).get("effect_counts") or {}
    summary_counts = ((summary.get("candidate_interpretation") or {}).get("effect_counts") or {})
    if summary_counts and dict(summary_counts) != dict(row_counts):
        issues.append(f"effect_counts mismatch: summary={summary_counts}, rows={row_counts}")
    if not summary_counts and rows:
        warnings.append("summary missing candidate_interpretation.effect_counts")

    for index, row in enumerate(rows, start=1):
        normalized = normalize_candidate_explanation_row(row)
        if normalized["effect_bucket"] == "未解析" and row.get("status") not in {"unresolved_gene", "unresolved_reaction", "missing_reaction"}:
            warnings.append(f"row {index} classified as unresolved without unresolved status")
        if row.get("intervention_type") in {"OE_gene_proxy", "OE_reaction"} and not row.get("resolved_reaction_id") and row.get("success") not in {False, "False", "false"}:
            warnings.append(f"row {index} OE row has no resolved_reaction_id")
        if row.get("status") == "2" and normalized["effect_bucket"] != "约束不可行":
            issues.append(f"row {index} status=2 is not classified as infeasible")

    if markdown_report:
        if "Python 草稿" not in markdown_report and "Python draft" not in markdown_report:
            issues.append("markdown report does not state Python draft status")
        if "MATLAB" not in markdown_report:
            issues.append("markdown report does not mention MATLAB alignment boundary")
        if "候选解释" not in markdown_report:
            issues.append("markdown report does not include candidate interpretation section")

    return {
        "passed": not issues,
        "issue_count": len(issues),
        "warning_count": len(warnings),
        "issues": issues,
        "warnings": warnings,
        "candidate_count": len(rows),
        "effect_counts": row_counts,
        "required_columns": list(CANDIDATE_COLUMNS),
    }


def normalize_candidate_explanation_row(row: dict[str, Any]) -> dict[str, Any]:
    intervention_type = _text(row.get("intervention_type"), "扰动")
    candidate_id = _text(
        row.get("input_gene_id")
        or row.get("gene_id")
        or row.get("resolved_reaction_id")
        or row.get("reaction_id")
        or row.get("candidate_id"),
        "未解析候选",
    )
    process = _text(row.get("secretory_process"), "未解析")
    effect_label = _effect_bucket(row)
    delta = row.get("delta_objective")
    status = _text(row.get("status"))
    solver_status = _text(row.get("solver_status_label") or row.get("failure_reason"))
    delta_text = _delta_text(delta)
    status_note = _status_note(status, solver_status, effect_label)
    summary = (
        f"{intervention_type} `{candidate_id}`：{effect_label}；"
        f"关联环节：{process}；Δobjective={delta_text}。{status_note}"
    )
    return {
        "candidate_id": candidate_id,
        "intervention_type": intervention_type,
        "effect_bucket": effect_label,
        "secretory_process": process,
        "delta_objective": delta,
        "status": status,
        "solver_status_label": solver_status,
        "summary": summary,
    }


def _effect_bucket(row: dict[str, Any]) -> str:
    status = _text(row.get("status"))
    effect = _text(row.get("effect_label"))
    solver_status = _text(row.get("solver_status_label") or row.get("failure_reason"))
    if status == "2" or "不可行" in effect or "不可行" in solver_status:
        return "约束不可行"
    if "提升" in effect:
        return "提升分泌"
    if "降低" in effect:
        return "降低分泌"
    if "无明显" in effect:
        return "无明显变化"
    if "未解析" in effect or status in {"unresolved_gene", "unresolved_reaction", "missing_reaction"}:
        return "未解析"
    if str(row.get("success")).lower() not in {"true", "1"}:
        return "求解失败"
    return effect or "未解析"


def _candidate_explanation_sort_key(row: dict[str, Any]) -> tuple[int, float, str]:
    priority = {
        "提升分泌": 0,
        "降低分泌": 1,
        "约束不可行": 2,
        "求解失败": 3,
        "未解析": 4,
        "无明显变化": 5,
    }.get(str(row.get("effect_bucket")), 9)
    delta = row.get("delta_objective")
    try:
        magnitude = abs(float(delta)) if delta is not None and str(delta) else 0.0
    except (TypeError, ValueError):
        magnitude = 0.0
    return (priority, -magnitude, str(row.get("candidate_id") or ""))


def _status_note(status: str, solver_status: str, effect_label: str) -> str:
    if status == "2" or effect_label == "约束不可行":
        return "固定生长下约束不可行，需要降低生长点、放宽扰动或单独诊断。"
    if effect_label == "未解析":
        return "候选 ID 未解析到当前模型对象。"
    if solver_status and solver_status != "求解成功":
        return f"求解状态：{solver_status}。"
    return "该解释仅表示当前模型约束下的相对分泌能力变化。"


def _delta_text(value: Any) -> str:
    try:
        if value is None or str(value) == "":
            return "无可行目标值"
        return f"{float(value):.6g}"
    except (TypeError, ValueError):
        return str(value)


def _text(value: Any, fallback: str = "") -> str:
    if value is None:
        return fallback
    text = str(value)
    if text.lower() in {"", "nan", "none", "nat"}:
        return fallback
    return text


def summarize_report_outputs(result: ReportOutputResult) -> dict[str, Any]:
    return {
        "output_dir": str(result.output_dir),
        "summary_path": str(result.summary_path),
        "report_path": str(result.report_path),
        "candidate_table_path": str(result.candidate_table_path),
        "tradeoff_path": str(result.tradeoff_path),
        "result_status": result.result_status,
        "matlab_alignment_status": result.matlab_alignment_status,
        "written_files": [str(path) for path in result.written_files],
    }


def _safe_prefix(value: str) -> str:
    cleaned = "".join(ch if ch.isalnum() or ch in ("-", "_") else "_" for ch in value.strip())
    return cleaned or "pcsec_pichia"


def _csv_value(value: Any) -> Any:
    if isinstance(value, (list, tuple, dict)):
        return json.dumps(value, ensure_ascii=False)
    return value


__all__ = [
    "CANDIDATE_COLUMNS",
    "CORE_CANDIDATE_EXPLANATION_COLUMNS",
    "ReportOutputResult",
    "TRADEOFF_COLUMNS",
    "build_markdown_report",
    "build_candidate_table",
    "build_candidate_interpretation",
    "build_summary_payload",
    "audit_candidate_report_consistency",
    "candidate_table_row",
    "classify_candidate_effect",
    "classify_secretory_process",
    "format_candidate_rows",
    "format_tradeoff_rows",
    "normalize_candidate_explanation_row",
    "summarize_report_outputs",
    "write_candidate_table",
    "write_outputs",
    "write_simulation_outputs",
    "write_tradeoff_table",
]
