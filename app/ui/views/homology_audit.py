from __future__ import annotations

from typing import Any

import pandas as pd
import streamlit as st

from app.services.pichia_external_reference_service import (
    export_external_reference_rows,
)
from app.services.pichia_homology_audit_service import (
    export_homology_audit_rows,
    load_homology_audit_browser_data,
)
from app.ui.common import HOMOLOGY_AUDIT_PAGE, request_navigation
from app.ui.views.simulation_gene_text import merge_candidate_text


ALL_OPTION = "全部"
SIMULATION_PAGE = "仿真验证"
READY_RULE_TRANSFER_STATUS = "rule_transfer_ready"

REVIEW_STATUS_OPTIONS = [
    "model_ready_rbh_high_confidence",
    "rbh_not_in_model",
    "low_identity_review_required",
    "coverage_review_required",
    "paralog_risk_review_required",
    "no_reciprocal_hit",
    "unresolved_query_symbol",
    "manual_review_required",
]
NAME_STATUS_OPTIONS = [
    "name_confirmed_by_rbh",
    "alias_confirmed_by_rbh",
    "pichia_locus_confirmed_by_rbh",
    "sequence_name_conflict",
    "external_name_missing",
    "internal_name_missing",
    "paralog_risk_review_required",
]
RULE_TRANSFER_STATUS_OPTIONS = [
    "rule_transfer_ready",
    "rule_transfer_supported_not_model_operable",
    "rule_transfer_low_confidence",
    "rule_transfer_paralog_risk",
    "rule_transfer_unresolved",
    "rule_transfer_not_supported",
]
TRISTATE_OPTIONS = {ALL_OPTION: None, "是": True, "否": False}

# 同源审计页领域枚举 → 中文（未知值回退原文）；与候选面板 _CANDIDATE_ENUM_LABELS 的重叠项用词一致。
_HOMOLOGY_ENUM_LABELS = {
    # review_status / homology_review_status
    "model_ready_rbh_high_confidence": "RBH 高置信·模型就绪",
    "rbh_not_in_model": "RBH 命中但不在模型内",
    "low_identity_review_required": "低相似度·需复核",
    "coverage_review_required": "覆盖度·需复核",
    "paralog_risk_review_required": "旁系同源风险·需复核",
    "no_reciprocal_hit": "无双向最佳命中",
    "unresolved_query_symbol": "查询符号未解析",
    "manual_review_required": "需人工复核",
    # name_consistency_status
    "name_confirmed_by_rbh": "RBH 确认命名",
    "alias_confirmed_by_rbh": "RBH 确认别名",
    "pichia_locus_confirmed_by_rbh": "RBH 确认 Pichia 基因座",
    "sequence_name_conflict": "序列-命名冲突",
    "external_name_missing": "缺外部命名",
    "internal_name_missing": "缺内部命名",
    # rule_transfer_status
    "rule_transfer_ready": "规则迁移·就绪",
    "rule_transfer_supported_not_model_operable": "规则迁移支持·但模型不可操作",
    "rule_transfer_low_confidence": "规则迁移·低置信",
    "rule_transfer_paralog_risk": "规则迁移·旁系同源风险",
    "rule_transfer_unresolved": "规则迁移·未解析",
    "rule_transfer_not_supported": "规则迁移·不支持",
    # external_crosscheck_status / 通用
    "not_available": "不可用",
    "available": "可用",
}

NAME_AUDIT_COLUMNS = {
    "internal_common_name": "内部常用名",
    "internal_gene_id": "内部基因ID",
    "internal_sequence_id": "SCE序列ID",
    "external_gene_name": "Pichia基因名",
    "external_locus_tag": "Pichia locus tag",
    "external_accession": "外部accession",
    "external_crosscheck_status": "外部库校验状态",
    "external_crosscheck_sources": "外部库来源",
    "external_crosscheck_warnings": "外部库告警",
    "identity_pct": "序列一致性 %",
    "query_coverage": "查询覆盖 %",
    "subject_coverage": "目标覆盖 %",
    "evalue": "e-value",
    "is_rbh": "RBH",
    "in_model_gene_index": "在模型gene_index",
    "name_consistency_status": "命名状态",
    "review_status": "同源审计状态",
    "warnings": "警告",
}
RULE_TRANSFER_COLUMNS = {
    "internal_common_name": "内部常用名",
    "query_symbol": "SCE符号",
    "sce_orf": "SCE ORF",
    "pichia_gene_id": "Pichia同源候选",
    "pichia_model_gene_id": "Pichia模型基因ID",
    "identity_pct": "序列一致性 %",
    "query_coverage": "查询覆盖 %",
    "subject_coverage": "目标覆盖 %",
    "evalue": "e-value",
    "is_rbh": "RBH",
    "in_model_gene_index": "在模型gene_index",
    "homology_review_status": "同源审计状态",
    "rule_transfer_status": "规则迁移状态",
    "warnings": "警告",
}
EXTERNAL_REFERENCE_COLUMNS = {
    "evidence_kind": "证据类型",
    "source_database": "来源库",
    "source_version": "来源版本",
    "retrieved_at": "抓取时间",
    "gene_id": "基因ID",
    "gene_name": "基因名",
    "pichia_gene_id": "Pichia 基因",
    "query_gene_id": "查询基因",
    "protein_name": "蛋白名",
    "function_description": "功能注释",
    "evidence_confidence": "证据置信度",
    "external_model_sources": "外部模型来源",
    "gpr_source_priority": "GPR 来源优先级",
    "external_gpr_candidate_count": "GPR 候选数",
    "best_external_gpr_source": "最佳外部 GPR 来源",
    "source_reaction_id": "来源反应",
    "source_gene_rule": "来源基因规则",
    "mapped_model_reaction_id": "映射到的模型反应",
    "external_gpr_mapping_status": "外部 GPR 映射状态",
    "external_gpr_conflict_warnings": "GPR 冲突告警",
    "gpr_transfer_status": "GPR 迁移状态",
    "manual_review_reasons": "人工复核原因",
}


def apply_homology_audit_prefill(gene_id: str) -> None:
    """Pre-fill the text-search filter with gene_id and jump here.

    Called from other pages (e.g. the genome-wide screen result table) so a
    reviewer can go straight from "this candidate looked interesting" to its
    homology/naming evidence without retyping the gene id.
    """
    st.session_state["homology_query"] = gene_id
    request_navigation(HOMOLOGY_AUDIT_PAGE)
    st.rerun()


def render_homology_audit() -> None:
    st.header(HOMOLOGY_AUDIT_PAGE)
    st.caption("只读展示离线 BLAST/RBH cache 产物，用于研发复核基因命名和同源规则迁移边界。")
    st.markdown(
        """
        <div class="concept-box">
        BLAST/RBH 是序列同源证据，不等同于表型证据；本页不会自动修改 SECRETION_GENE_CATALOG，
        也不会自动证明 KO/OE 在当前模型或湿实验中有效。模型可执行性仍需结合 GPR 规则和实际求解行为复核。
        </div>
        """,
        unsafe_allow_html=True,
    )

    filters = _render_filters()
    payload = load_homology_audit_browser_data(**filters)
    cache_status = payload.get("cache_status", {})
    if not cache_status.get("cache_available"):
        _render_missing_cache(cache_status)
        return

    summary = payload.get("summary", {})
    name_rows = list(payload.get("name_audit_rows", []))
    rule_rows = list(payload.get("rule_transfer_audit_rows", []))
    external_rows = list(payload.get("external_reference_rows", []))

    _render_summary_metrics(summary, name_rows, rule_rows)
    name_tab, rule_tab, cache_tab = st.tabs(["命名标准化", "同源规则迁移评估", "缓存状态与导出"])
    with name_tab:
        _render_rows_table("命名标准化", name_rows, NAME_AUDIT_COLUMNS)
    with rule_tab:
        _render_rows_table("同源规则迁移评估", rule_rows, RULE_TRANSFER_COLUMNS)
        _render_add_to_simulation_controls(rule_rows)
    with cache_tab:
        _render_cache_and_export(cache_status, name_rows, rule_rows, external_rows)


def _render_filters() -> dict[str, Any]:
    st.subheader("筛选")
    query = st.text_input("文本搜索", placeholder="KAR2、YJL034W、PAS_chr2-1_0140", key="homology_query")
    col_review, col_name, col_rule = st.columns(3)
    review_status = _status_select(col_review, "同源审计状态", REVIEW_STATUS_OPTIONS, "homology_review_status")
    name_status = _status_select(col_name, "命名状态", NAME_STATUS_OPTIONS, "homology_name_status")
    rule_status = _status_select(col_rule, "规则迁移状态", RULE_TRANSFER_STATUS_OPTIONS, "homology_rule_status")
    col_rbh, col_model, col_identity = st.columns(3)
    is_rbh = _tristate_select(col_rbh, "RBH", "homology_is_rbh")
    in_model = _tristate_select(col_model, "在模型gene_index", "homology_in_model")
    min_identity = col_identity.slider("最小 identity %", min_value=0.0, max_value=100.0, value=0.0, step=1.0)
    return {
        "query": query,
        "review_status": review_status,
        "name_consistency_status": name_status,
        "rule_transfer_status": rule_status,
        "is_rbh": is_rbh,
        "in_model_gene_index": in_model,
        "min_identity": min_identity if min_identity > 0 else None,
    }


def _status_select(column: Any, label: str, options: list[str], key: str) -> str | None:
    selected = column.selectbox(label, [ALL_OPTION, *options], index=0, key=key)
    return None if selected == ALL_OPTION else selected


def _tristate_select(column: Any, label: str, key: str) -> bool | None:
    selected = column.selectbox(label, list(TRISTATE_OPTIONS), index=0, key=key)
    return TRISTATE_OPTIONS[str(selected)]


def _render_missing_cache(cache_status: dict[str, Any]) -> None:
    st.warning("尚未找到可用的同源审计 cache。请先离线生成 cache 后再回到本页查看。")
    st.caption("页面打开时不会运行 BLAST、不会联网，也不会写入 catalog。")
    command = str(cache_status.get("recommended_build_command") or "python scripts\\build_pichia_homology_cache.py --catalog-only")
    st.code(command, language="powershell")
    missing_files = cache_status.get("missing_files") or []
    if missing_files:
        st.info("缺失文件：" + ", ".join(str(item) for item in missing_files))


def _render_summary_metrics(
    summary: dict[str, Any],
    name_rows: list[dict[str, Any]],
    rule_rows: list[dict[str, Any]],
) -> None:
    values = _summary_metric_values(summary, name_rows, rule_rows)
    columns = st.columns(6)
    columns[0].metric("总行数", values["total"])
    columns[1].metric("规则迁移就绪", values["ready"])
    columns[2].metric("支持但模型不可操作", values["supported_not_model_operable"])
    columns[3].metric("低置信", values["low_confidence"])
    columns[4].metric("未解析", values["unresolved"])
    columns[5].metric("旁系风险/不支持", values["paralog_or_not_supported"])


def _summary_metric_values(
    summary: dict[str, Any],
    name_rows: list[dict[str, Any]],
    rule_rows: list[dict[str, Any]],
) -> dict[str, int]:
    raw_counts = summary.get("rule_transfer_status_counts")
    counts = raw_counts if isinstance(raw_counts, dict) else _count_rows(rule_rows, "rule_transfer_status")
    total = int(
        summary.get("rule_transfer_row_count")
        or summary.get("name_audit_row_count")
        or len(rule_rows)
        or len(name_rows)
    )
    return {
        "total": total,
        "ready": int(counts.get("rule_transfer_ready", 0)),
        "supported_not_model_operable": int(counts.get("rule_transfer_supported_not_model_operable", 0)),
        "low_confidence": int(counts.get("rule_transfer_low_confidence", 0)),
        "unresolved": int(counts.get("rule_transfer_unresolved", 0)),
        "paralog_or_not_supported": int(counts.get("rule_transfer_paralog_risk", 0))
        + int(counts.get("rule_transfer_not_supported", 0)),
    }


def _count_rows(rows: list[dict[str, Any]], field: str) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in rows:
        value = str(row.get(field) or "")
        if value:
            counts[value] = counts.get(value, 0) + 1
    return counts


def _render_rows_table(title: str, rows: list[dict[str, Any]], columns: dict[str, str]) -> None:
    st.subheader(f"{title}（{len(rows)} 行）")
    if not rows:
        st.info("当前筛选条件下没有可显示的行。")
        return
    st.dataframe(_rows_to_frame(rows, columns), width='stretch', hide_index=True)


def _render_add_to_simulation_controls(rule_rows: list[dict[str, Any]]) -> None:
    st.markdown("**加入仿真验证输入**")
    st.caption(
        "这里只把已通过规则迁移评估、且存在 Pichia 模型 gene_id 的同源候选加入仿真输入；"
        "这不等于推荐、不等于表型证据，也不会自动运行仿真。"
    )
    ready_rows = _ready_rule_transfer_rows(rule_rows)
    if not ready_rows:
        st.info("当前筛选结果中没有可加入仿真验证的模型可操作同源候选。")
        return

    options = [_rule_transfer_option_label(row) for row in ready_rows]
    selected_labels = st.multiselect(
        "选择要加入仿真验证的模型基因",
        options,
        key="homology_ready_rule_transfer_selection",
    )
    selected_rows = [row for row, label in zip(ready_rows, options) if label in selected_labels]
    add_ko_col, add_oe_col = st.columns(2)
    with add_ko_col:
        if st.button("添加到 KO 输入并跳转仿真验证", key="homology_add_ready_ko", disabled=not selected_rows):
            _apply_rule_transfer_selection(selected_rows, action="ko")
    with add_oe_col:
        if st.button("添加到 OE 输入并跳转仿真验证", key="homology_add_ready_oe", disabled=not selected_rows):
            _apply_rule_transfer_selection(selected_rows, action="oe")


def _ready_rule_transfer_rows(rule_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    ready: list[dict[str, Any]] = []
    seen_gene_ids: set[str] = set()
    for row in rule_rows:
        model_gene_id = str(row.get("pichia_model_gene_id") or "").strip()
        if (
            row.get("rule_transfer_status") == READY_RULE_TRANSFER_STATUS
            and bool(row.get("in_model_gene_index"))
            and model_gene_id
            and model_gene_id not in seen_gene_ids
        ):
            seen_gene_ids.add(model_gene_id)
            ready.append(row)
    return ready


def _apply_rule_transfer_selection(rows: list[dict[str, Any]], *, action: str) -> None:
    model_gene_ids = _selected_model_gene_ids(rows)
    if not model_gene_ids:
        st.info("所选候选没有可加入仿真验证的 Pichia 模型 gene_id。")
        return
    key = f"pichia_draft_{action}_genes"
    st.session_state[key] = merge_candidate_text(str(st.session_state.get(key, "")), model_gene_ids)
    action_label = "KO" if action == "ko" else "OE"
    st.toast(f"已加入 {action_label} 输入：{', '.join(model_gene_ids)}；正在跳转到仿真验证。")
    request_navigation(SIMULATION_PAGE)
    st.rerun()


def _selected_model_gene_ids(rows: list[dict[str, Any]]) -> list[str]:
    return [str(row.get("pichia_model_gene_id") or "").strip() for row in _ready_rule_transfer_rows(rows)]


def _rule_transfer_option_label(row: dict[str, Any]) -> str:
    model_gene_id = str(row.get("pichia_model_gene_id") or "").strip()
    name = str(row.get("internal_common_name") or row.get("query_symbol") or row.get("sce_orf") or "").strip()
    return f"{model_gene_id} — {name}" if name else model_gene_id


def _rows_to_frame(rows: list[dict[str, Any]], columns: dict[str, str]) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {label: _display_value(row.get(key)) for key, label in columns.items()}
            for row in rows
        ]
    )


def _display_value(value: Any) -> Any:
    if isinstance(value, bool):
        return "是" if value else "否"
    if isinstance(value, list):
        return "; ".join(str(item) for item in value)
    if isinstance(value, dict):
        return "; ".join(f"{key}:{item}" for key, item in sorted(value.items()))
    if value is None:
        return ""
    if isinstance(value, str) and value in _HOMOLOGY_ENUM_LABELS:
        return _HOMOLOGY_ENUM_LABELS[value]
    return value


def _render_cache_and_export(
    cache_status: dict[str, Any],
    name_rows: list[dict[str, Any]],
    rule_rows: list[dict[str, Any]],
    external_rows: list[dict[str, Any]] | None = None,
) -> None:
    st.subheader("缓存状态")
    st.markdown(
        f"""
        - 缓存根目录: `{cache_status.get("cache_root", "")}`
        - 生成时间: `{cache_status.get("generated_at", "")}`
        - 行数: `{cache_status.get("row_count", 0)}`
        """
    )
    resolved_external_rows = external_rows or []
    _render_external_cache_status(cache_status, name_rows, resolved_external_rows)
    st.caption("导出的是当前筛选结果；导出动作不会重新运行 BLAST，也不会写回任何模型或 catalog。")
    export_kind = st.selectbox("导出表", ["命名标准化", "同源规则迁移评估", "外部数据库证据"], key="homology_export_kind")
    export_format = st.selectbox("导出格式", ["TSV", "CSV"], key="homology_export_format")
    rows = name_rows if export_kind == "命名标准化" else rule_rows
    exporter = export_homology_audit_rows
    if export_kind == "外部数据库证据":
        rows = resolved_external_rows
        exporter = export_external_reference_rows
    suffix = export_format.lower()
    st.download_button(
        "下载当前筛选结果",
        exporter(rows, file_format=suffix),
        file_name=f"pichia_homology_{_export_slug(export_kind)}.{suffix}",
        mime="text/tab-separated-values" if suffix == "tsv" else "text/csv",
    )


def _render_external_cache_status(
    cache_status: dict[str, Any],
    name_rows: list[dict[str, Any]],
    external_rows: list[dict[str, Any]],
) -> None:
    available = "可用" if cache_status.get("external_cache_available") else "不可用"
    source_counts = cache_status.get("external_source_counts")
    source_text = _format_source_counts(source_counts if isinstance(source_counts, dict) else {})
    st.markdown(
        f"""
        **外部命名参考缓存**

        - 状态: `{available}`
        - 缓存路径: `{cache_status.get("external_cache_path", "")}`
        - 生成时间: `{cache_status.get("external_cache_generated_at", "")}`
        - 参考条目数: `{cache_status.get("external_reference_count", 0)}`
        - 来源: `{source_text}`
        """
    )
    command = str(cache_status.get("recommended_external_build_command") or "")
    if command:
        st.code(command, language="powershell")
    warnings = cache_status.get("external_cache_warnings") or []
    if warnings:
        st.info("外部缓存告警: " + "; ".join(str(item) for item in warnings))

    counts = _count_rows(name_rows, "external_crosscheck_status")
    if counts:
        count_lines = "\n".join(
            f"- {_HOMOLOGY_ENUM_LABELS.get(status, status)}: {count}" for status, count in sorted(counts.items())
        )
    else:
        count_lines = "- 不可用: 0"
    st.markdown(
        "**外部库校验状态计数**\n\n"
        f"{count_lines}\n\n"
        "外部库校验只是离线命名参考：不改变 RBH 判定、KO/OE 推荐分级或表型证据，也不会在 Streamlit 里联网查询。"
    )
    _render_external_reference_cache_status(cache_status, external_rows)


def _render_external_reference_cache_status(
    cache_status: dict[str, Any],
    external_rows: list[dict[str, Any]],
) -> None:
    status = cache_status.get("external_reference_cache")
    if not isinstance(status, dict):
        status = {}
    source_counts = status.get("source_counts") if isinstance(status.get("source_counts"), dict) else {}
    type_counts = status.get("record_type_counts") if isinstance(status.get("record_type_counts"), dict) else {}
    retrieved = status.get("retrieved_at_range") if isinstance(status.get("retrieved_at_range"), dict) else {}
    gpr_priority = status.get("gpr_source_priority") if isinstance(status.get("gpr_source_priority"), dict) else {}
    st.markdown(
        f"""
        **外部参考缓存**

        - 状态: `{"可用" if status.get("cache_available") else "不可用"}`
        - 记录路径: `{status.get("records_path", "")}`
        - 记录数: `{status.get("record_count", 0)}`
        - 来源: `{_format_source_counts(source_counts)}`
        - 记录类型: `{_format_source_counts(type_counts)}`
        - 抓取时间: `{retrieved.get("first", "")}` 至 `{retrieved.get("last", "")}`
        - 外部模型来源: `{_display_value(status.get("external_model_sources") or [])}`
        - 外部 GPR 候选数: `{status.get("external_gpr_candidate_count", 0)}`
        - 最佳外部 GPR 来源: `{status.get("best_external_gpr_source", "")}`
        - 最佳 GPR 优先级档: `{gpr_priority.get("best_priority_tier", "")}`
        """
    )
    command = str(status.get("recommended_refresh_command") or "")
    if command:
        st.code(command, language="powershell")
    warnings = status.get("warnings") or []
    if warnings:
        st.info("外部参考告警: " + "; ".join(str(item) for item in warnings))
    if external_rows:
        st.dataframe(_rows_to_frame(external_rows, EXTERNAL_REFERENCE_COLUMNS), width='stretch', hide_index=True)
    else:
        st.info("当前外部参考缓存中没有可显示的基因功能或 GPR 候选行。")


def _export_slug(export_kind: str) -> str:
    if export_kind == "命名标准化":
        return "name_audit"
    if export_kind == "同源规则迁移评估":
        return "rule_transfer_audit"
    return "external_reference_evidence"


def _format_source_counts(source_counts: dict[str, Any]) -> str:
    if not source_counts:
        return "none"
    return ", ".join(f"{source}:{count}" for source, count in sorted(source_counts.items()))


__all__ = [
    "HOMOLOGY_AUDIT_PAGE",
    "NAME_AUDIT_COLUMNS",
    "RULE_TRANSFER_COLUMNS",
    "EXTERNAL_REFERENCE_COLUMNS",
    "render_homology_audit",
]
