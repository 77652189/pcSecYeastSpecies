from __future__ import annotations

from typing import Any

import pandas as pd
import streamlit as st

from app.services.pichia_gene_catalog_service import (
    hlf_opn_candidate_gene_summary,
    hlf_opn_executable_candidate_inputs,
    load_hlf_opn_candidate_genes,
    load_hlf_opn_gpr_overlay_review,
)
from app.ui.views.simulation_gene_text import merge_candidate_text


ALL_FILTER = "全部"


# hLF/OPN 候选 + overlay 复核面板的领域枚举 → 中文（本地 domain 字典，未知值回退原文）。
_CANDIDATE_ENUM_LABELS = {
    # operability_status
    "model_ko_executable": "模型可执行 KO",
    "model_oe_proxy_executable": "模型可执行 OE（反应代理）",
    "model_explain_only_complex_subunit": "仅解释：复合体亚基",
    "not_in_model": "不在模型内",
    "unresolved_name": "命名未解析",
    "manual_review_required": "需人工复核",
    # recommended_intervention
    "KO": "敲除(KO)",
    "OE": "过表达(OE)",
    "review_only": "仅复核",
    # evidence_type
    "homology_auxiliary": "同源辅助证据",
    "model_reaction_proxy": "模型反应代理",
    "phenotype_or_literature_evidence": "表型/文献证据",
    # evidence_confidence
    "curated_evidence_model_executable": "策展证据·模型可执行",
    "homology_supported_review_required": "同源支持·需复核",
    "model_executable_annotation_supported": "模型可执行·注释支持",
    "reaction_proxy_review_required": "反应代理·需复核",
    "high_exact_multi_source": "高·多源精确匹配",
    "unresolved": "未解析",
    # homology_review_status
    "low_identity_review_required": "低相似度·需复核",
    "model_ready_rbh_high_confidence": "RBH 高置信·模型就绪",
    "no_reciprocal_hit": "无双向最佳命中",
    "rbh_not_in_model": "RBH 命中但不在模型内",
    # rule_transfer_status
    "rule_transfer_low_confidence": "规则迁移·低置信",
    "rule_transfer_not_supported": "规则迁移·不支持",
    "rule_transfer_ready": "规则迁移·就绪",
    "rule_transfer_supported_not_model_operable": "规则迁移支持·但模型不可操作",
    # review_status / recommended_action（overlay）
    "shared": "共享",
    "candidate_gpr_overlay_review": "候选 GPR overlay 复核",
    "eligible_for_overlay_if_all_complex_subunits_are_confirmed": "全部复合体亚基确认后可 overlay",
    "manual_locus_review_required": "需人工核对基因座",
}


def _label(value: object) -> str:
    text = str(value or "")
    return _CANDIDATE_ENUM_LABELS.get(text, text)


def render_hlf_opn_candidate_panel(target_id: str) -> None:
    context = target_context_for_hlf_opn_candidates(target_id)
    with st.expander("hLF / OPN 候选基因", expanded=context is not None):
        if context is None:
            st.caption("当前目标不属于 hLF / OPN 专用候选上下文。")
            return

        candidates = load_hlf_opn_candidate_genes(target_context=context, include_shared=True)
        executable = hlf_opn_executable_candidate_inputs(target_context=context, include_shared=True)
        overlay_rows = load_hlf_opn_gpr_overlay_review(target_context=context, include_shared=True)
        summary = hlf_opn_candidate_gene_summary()

        _render_candidate_metrics(context, candidates, executable, overlay_rows, summary)
        filtered_candidates = _filter_candidates_for_display(candidates, context)
        candidate_frame = _candidate_frame(filtered_candidates)
        st.dataframe(candidate_frame, use_container_width=True, hide_index=True)
        st.download_button(
            f"导出 {context} 候选表",
            candidate_frame.to_csv(index=False).encode("utf-8-sig"),
            file_name=f"{context.lower()}_candidate_genes.csv",
            mime="text/csv",
            key=f"hlf_opn_export_candidates_{context}",
        )

        disabled = not (executable.get("ko_gene_ids") or executable.get("oe_gene_ids"))
        if st.button(
            f"加入 {context} 可执行 KO/OE 输入",
            key=f"hlf_opn_add_executable_{context}",
            disabled=disabled,
        ):
            added = _apply_executable_candidate_inputs(executable)
            _toast_added_inputs(context, added, int(executable.get("excluded_count") or 0))

        for warning in executable.get("warnings") or []:
            st.caption(str(warning))

        if overlay_rows:
            with st.expander("模型外候选 overlay 复核", expanded=False):
                overlay_frame = _overlay_frame(overlay_rows)
                st.dataframe(overlay_frame, use_container_width=True, hide_index=True)
                st.download_button(
                    f"导出 {context} overlay 复核表",
                    overlay_frame.to_csv(index=False).encode("utf-8-sig"),
                    file_name=f"{context.lower()}_overlay_review.csv",
                    mime="text/csv",
                    key=f"hlf_opn_export_overlay_{context}",
                )


def target_context_for_hlf_opn_candidates(target_id: str) -> str | None:
    text = str(target_id or "").strip()
    upper = text.upper()
    if upper == "HLF" or "HLF" in upper:
        return "hLF"
    if upper == "OPN" or "OPN" in upper:
        return "OPN"
    return None


def _render_candidate_metrics(
    context: str,
    candidates: list[dict[str, Any]],
    executable: dict[str, Any],
    overlay_rows: list[dict[str, Any]],
    summary: dict[str, Any],
) -> None:
    metrics = st.columns(5)
    metrics[0].metric(f"{context} 候选", len(candidates))
    metrics[1].metric("KO 可执行", len(executable.get("ko_gene_ids") or []))
    metrics[2].metric("OE 代理可执行", len(executable.get("oe_gene_ids") or []))
    metrics[3].metric("仅复核/模型外", int(executable.get("excluded_count") or 0))
    metrics[4].metric("模型外候选", len(overlay_rows))
    target_counts = summary.get("target_candidate_counts") if isinstance(summary.get("target_candidate_counts"), dict) else {}
    if target_counts:
        st.caption(
            "候选统计："
            f"hLF {target_counts.get('hLF', 0)}；"
            f"OPN {target_counts.get('OPN', 0)}。"
        )


def _filter_candidates_for_display(
    candidates: list[dict[str, Any]],
    context: str,
) -> list[dict[str, Any]]:
    statuses = [ALL_FILTER, *sorted({str(row.get("operability_status") or "") for row in candidates if row})]
    selected = st.selectbox(
        "模型可操作性",
        statuses,
        key=f"hlf_opn_operability_filter_{context}",
    )
    if selected == ALL_FILTER:
        return candidates
    return [row for row in candidates if str(row.get("operability_status") or "") == selected]


def _apply_executable_candidate_inputs(executable: dict[str, Any]) -> dict[str, list[str]]:
    ko_gene_ids = [str(item).strip() for item in executable.get("ko_gene_ids") or [] if str(item).strip()]
    oe_gene_ids = [str(item).strip() for item in executable.get("oe_gene_ids") or [] if str(item).strip()]
    if ko_gene_ids:
        st.session_state["pichia_draft_ko_genes"] = merge_candidate_text(
            str(st.session_state.get("pichia_draft_ko_genes", "")),
            ko_gene_ids,
        )
    if oe_gene_ids:
        st.session_state["pichia_draft_oe_genes"] = merge_candidate_text(
            str(st.session_state.get("pichia_draft_oe_genes", "")),
            oe_gene_ids,
        )
    return {"ko": ko_gene_ids, "oe": oe_gene_ids}


def _toast_added_inputs(context: str, added: dict[str, list[str]], excluded_count: int) -> None:
    parts: list[str] = []
    if added["ko"]:
        parts.append(f"KO {len(added['ko'])}")
    if added["oe"]:
        parts.append(f"OE {len(added['oe'])}")
    detail = "，".join(parts) if parts else "0"
    st.toast(f"已加入 {context} 可执行候选：{detail}；{excluded_count} 个候选保留为复核/导出。")


def _candidate_frame(rows: list[dict[str, Any]]) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "目标": _label(row.get("target_context", "")),
                "基因": row.get("gene_id", ""),
                "标准命名": row.get("display_name", ""),
                "标准符号": row.get("standard_symbol", ""),
                "蛋白注释": row.get("protein_name", ""),
                "外部ID": _format_external_ids(row.get("external_ids")),
                "模型可操作性": _label(row.get("operability_status", "")),
                "扰动": _label(row.get("recommended_intervention", "")),
                "证据类型": _label(row.get("evidence_type", "")),
                "证据置信度": _label(row.get("evidence_confidence", "")),
                "同源状态": _label(row.get("homology_review_status", "")),
                "规则迁移": _label(row.get("rule_transfer_status", "")),
                "限制": _join_values(row.get("warnings")),
            }
            for row in rows
        ]
    )


def _overlay_frame(rows: list[dict[str, Any]]) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "目标": _label(row.get("target_context", "")),
                "基因": row.get("gene_id", ""),
                "来源候选": row.get("source_common_name", ""),
                "复核状态": _label(row.get("review_status", "")),
                "证据置信度": _label(row.get("evidence_confidence", "")),
                "已有模型反应": _join_values(row.get("existing_model_reaction_ids")),
                "缺失模型反应": _join_values(row.get("missing_model_reaction_ids")),
                "建议动作": _label(row.get("recommended_action", "")),
                "风险": row.get("risk", ""),  # 自由文本说明（含模型反应名），不做字典映射
                "限制": _join_values(row.get("warnings")),
            }
            for row in rows
        ]
    )


def _format_external_ids(value: object) -> str:
    if not isinstance(value, dict):
        return ""
    return "; ".join(f"{key}:{item}" for key, item in value.items() if item)


def _join_values(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, (list, tuple)):
        return "; ".join(str(item) for item in value if str(item))
    return str(value)


__all__ = [
    "render_hlf_opn_candidate_panel",
    "target_context_for_hlf_opn_candidates",
]
