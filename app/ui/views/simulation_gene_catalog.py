from __future__ import annotations

import pandas as pd
import streamlit as st

from app.services.pichia_gene_catalog_service import (
    build_pichia_gene_evidence_cache,
    list_pichia_gene_rule_evidence,
    list_pichia_secretion_gene_evidence,
    list_verified_secretion_gene_library,
    load_pichia_full_model_gene_catalog,
)
from app.ui.views.simulation_gene_text import merge_candidate_text

KO_RUNNABLE_STATUS = "ko_runnable_gpr_gene_deletion"
OE_RUNNABLE_STATUS = "oe_runnable_reaction_proxy"
# Full-model catalog has ~1025 rows; cap how many query matches render at once so a broad
# keyword doesn't dump a huge table (narrow the search instead of scrolling/paginating).
FULL_MODEL_SEARCH_LIMIT = 30


def render_gene_lookup_panel() -> None:
    st.markdown("**基因 / 反应 ID 查找**")
    st.caption(
        "全基因组KO/OE筛查已经系统覆盖了全部1025个模型基因和策展库的反应级候选——"
        "找值得测试的候选，优先去「全基因组KO/OE筛查」结果里看，点候选行的"
        "「在仿真验证中核实」会自动填好下面的输入框。这里的搜索适合已经知道大致名字、"
        "想直接查到精确ID的场景。"
    )
    query = st.text_input(
        "搜索基因常用名 / locus / 模型基因ID / 反应ID",
        value=st.session_state.get("pichia_gene_lookup_query", ""),
        placeholder="例如：PDI1、ERO1、Kar2、PAS_chr2-2_0107、sec_Kar2p_complex_formation",
        key="pichia_gene_lookup_query",
    )
    if query.strip():
        _render_search_results(query)
    else:
        st.caption("输入关键词后显示匹配结果。")

    with st.expander("高级：外部证据 GPR overlay / 候选库维护", expanded=False):
        st.caption("外部证据 overlay 是实验性证据层，默认不进入仿真；当前无可执行补充规则时只作人工复核参考。")
        maintenance_col, overlay_col = st.columns([1.0, 1.0])
        with maintenance_col:
            if st.button(
                "刷新常用基因证据缓存",
                key="pichia_gene_refresh_lightweight_cache",
                help="刷新策展候选及反应代理缓存；不会重建全模型湿实验注释。",
            ):
                refreshed = list_pichia_secretion_gene_evidence("", force_refresh=True)
                st.success(f"常用基因证据缓存已刷新：{len(refreshed)} 条。")
            if st.button(
                "在线重建全模型湿实验注释缓存",
                key="pichia_gene_refresh_evidence",
                help="从 UniProt/KEGG 重新构建全部模型基因的湿实验注释；比较慢。",
            ):
                with st.spinner("正在从 UniProt/KEGG 构建湿实验注释缓存..."):
                    summary = build_pichia_gene_evidence_cache()
                st.success(
                    "湿实验注释缓存已更新："
                    f"{summary.get('database_supported_count', 0)} / {summary.get('total_genes', 0)} 个基因有数据库支持。"
                )
        with overlay_col:
            show_overlay = st.checkbox("显示外部证据 GPR overlay", value=False, key="pichia_gene_show_overlay")
        if show_overlay:
            _render_gene_rule_overlay_lookup(query)
        else:
            st.info("打开上方开关后查看 PDI1/ERO1/KAR2/OCH1/PEP4/PRB1 等 overlay 证据。")


def _render_search_results(query: str) -> None:
    rows, truncated_count = _collect_search_rows(query)
    if not rows:
        st.info("未找到匹配；可以尝试模型基因ID（如 PAS_chr...）或反应ID（如 sec_...）。")
        return
    st.dataframe(
        pd.DataFrame(
            [
                {
                    "来源": row["来源"],
                    "名称": row["名称"],
                    "类型": row["类型"],
                    "ID": row["ID"],
                    "可用于": row["可用于"],
                    "说明": row["说明"],
                }
                for row in rows
            ]
        ),
        use_container_width=True,
        hide_index=True,
    )
    if truncated_count:
        st.caption(f"全模型基因库还有 {truncated_count} 个匹配未显示，请输入更精确的关键词缩小范围。")

    options = [f"{row['名称']} — {row['ID']}" for row in rows]
    selected_labels = st.multiselect("选择候选（可多选，用于组合测试）", options, key="pichia_gene_search_sel")
    selected_rows = [row for row, label in zip(rows, options) if label in selected_labels]
    add_ko_col, add_oe_col = st.columns(2)
    with add_ko_col:
        if st.button("添加到敲除输入", key="pichia_gene_search_add_ko") and selected_rows:
            _apply_search_selection(selected_rows, action="ko")
    with add_oe_col:
        if st.button("添加到过表达输入", key="pichia_gene_search_add_oe") and selected_rows:
            _apply_search_selection(selected_rows, action="oe")
    message = st.session_state.pop("pichia_gene_catalog_message", "")
    if message:
        st.info(message)


def _collect_search_rows(query: str) -> tuple[list[dict[str, object]], int]:
    """Merge curated-library and full-model-gene matches into one result list.

    Both sources ultimately search the same underlying catalog data (curated rows already
    carry gene_id/ko_reaction_id/oe_reaction_id), so this is one search instead of the three
    separate browsers the page used to have.
    """
    rows: list[dict[str, object]] = []
    seen: set[tuple[str, str]] = set()

    def _add(source: str, name: str, kind: str, item_id: object, ko_ok: bool, oe_ok: bool, detail: str) -> None:
        item_id_str = str(item_id or "").strip()
        if not item_id_str or (kind, item_id_str) in seen:
            return
        seen.add((kind, item_id_str))
        usable = "/".join(label for label, ok in (("KO", ko_ok), ("OE", oe_ok)) if ok)
        rows.append(
            {
                "来源": source,
                "名称": name or item_id_str,
                "类型": "基因" if kind == "gene" else "反应",
                "ID": item_id_str,
                "可用于": usable or "仅解释",
                "说明": detail,
                "kind": kind,
            }
        )

    for row in list_verified_secretion_gene_library(query):
        name = str(row.get("display_name") or "")
        detail = str(row.get("function_annotation") or "")
        _add("策展库", name, "gene", row.get("model_gene_id"), True, True, detail)
        _add("策展库", name, "reaction", row.get("ko_reaction_id"), True, False, detail)
        _add("策展库", name, "reaction", row.get("oe_reaction_id"), False, True, detail)

    full_model_matches = _filter_full_model_gene_rows(load_pichia_full_model_gene_catalog(), query=query)
    truncated_count = max(0, len(full_model_matches) - FULL_MODEL_SEARCH_LIMIT)
    for gene in full_model_matches[:FULL_MODEL_SEARCH_LIMIT]:
        _add(
            "全模型",
            _full_model_gene_display_name(gene),
            "gene",
            gene.get("gene_id"),
            gene.get("ko_support_status") == KO_RUNNABLE_STATUS,
            gene.get("oe_support_status") == OE_RUNNABLE_STATUS,
            _full_model_gene_function_summary(gene),
        )
    return rows, truncated_count


def _partition_selection_by_kind(rows: list[dict[str, object]]) -> tuple[list[str], list[str]]:
    """Split selected search rows into (gene_ids, reaction_ids) for routing to the right input box."""
    genes = [str(row["ID"]) for row in rows if row.get("kind") == "gene" and row.get("ID")]
    reactions = [str(row["ID"]) for row in rows if row.get("kind") == "reaction" and row.get("ID")]
    return genes, reactions


def _apply_search_selection(rows: list[dict[str, object]], *, action: str) -> None:
    genes, reactions = _partition_selection_by_kind(rows)
    added: list[str] = []
    if genes:
        key = f"pichia_draft_{action}_genes"
        st.session_state[key] = merge_candidate_text(str(st.session_state.get(key, "")), genes)
        added.append(f"基因：{', '.join(genes)}")
    if reactions:
        key = f"pichia_draft_{action}_reactions"
        st.session_state[key] = merge_candidate_text(str(st.session_state.get(key, "")), reactions)
        added.append(f"反应：{', '.join(reactions)}")
    action_label = "敲除" if action == "ko" else "过表达"
    st.session_state["pichia_gene_catalog_message"] = (
        f"已加入{action_label}输入：" + "；".join(added) if added else "所选候选没有可用的 ID。"
    )
    st.rerun()


def _filter_full_model_gene_rows(
    rows: list[dict[str, object]],
    *,
    query: str = "",
    only_ko: bool = False,
    only_oe: bool = False,
    wet_lab_filter: str = "全部",
) -> list[dict[str, object]]:
    query_text = query.strip().lower()
    filtered_rows = [
        gene
        for gene in rows
        if (not only_ko or gene.get("ko_support_status") == KO_RUNNABLE_STATUS)
        and (not only_oe or gene.get("oe_support_status") == OE_RUNNABLE_STATUS)
        and _matches_wet_lab_filter(gene, wet_lab_filter)
    ]
    if not query_text:
        return filtered_rows
    return [gene for gene in filtered_rows if _full_model_gene_matches_query(gene, query_text)]


def _full_model_gene_matches_query(gene: dict[str, object], query_text: str) -> bool:
    return (
        query_text in str(gene.get("gene_id") or "").lower()
        or query_text in _full_model_gene_display_name(gene).lower()
        or query_text in _full_model_gene_function_summary(gene).lower()
        or query_text in str(gene.get("primary_category") or "").lower()
        or query_text in str(gene.get("processes") or "").lower()
        or query_text in " ".join(str(item) for item in gene.get("aliases") or []).lower()
        or query_text in str(gene.get("protein_name") or "").lower()
        or query_text in str(gene.get("function_annotation") or "").lower()
        or query_text in str(gene.get("ko_support_status") or "").lower()
        or query_text in str(gene.get("oe_support_status") or "").lower()
    )


def _full_model_gene_display_name(gene: dict[str, object]) -> str:
    display_name = str(gene.get("display_name") or "").strip()
    if display_name and display_name != str(gene.get("gene_id") or "").strip():
        return display_name
    standard_symbol = str(gene.get("standard_gene_symbol") or "").strip()
    if standard_symbol:
        return standard_symbol
    protein_name = str(gene.get("protein_name") or "").strip()
    if protein_name:
        return protein_name
    aliases = [str(item).strip() for item in gene.get("aliases") or [] if str(item).strip()]
    if aliases:
        return aliases[0]
    reactions = _reaction_tokens(gene)
    if reactions:
        return f"{'/'.join(reactions[:3])} 相关酶（未注释）"
    return "未注释模型基因"


def _full_model_gene_function_summary(gene: dict[str, object]) -> str:
    annotation = str(gene.get("function_annotation") or "").strip()
    if annotation:
        return annotation
    reactions = _reaction_tokens(gene)
    if reactions:
        joined = ", ".join(reactions[:5])
        suffix = " 等" if len(reactions) > 5 else ""
        return f"按模型 GPR 关联到反应：{joined}{suffix}。尚无外部基因名/功能注释。"
    return "尚无外部基因名/功能注释；仅可按模型 locus ID 和 GPR 关系使用。"


def _matches_wet_lab_filter(gene: dict[str, object], wet_lab_filter: str) -> bool:
    readiness = str(gene.get("wet_lab_readiness") or "model_only_not_experiment_ready")
    if wet_lab_filter == "可直接推进湿实验":
        return readiness == "database_supported_experiment_candidate"
    if wet_lab_filter == "需人工确认":
        return readiness == "manual_review_required"
    if wet_lab_filter == "仅模型级候选":
        return readiness == "model_only_not_experiment_ready"
    return True


def _reaction_tokens(gene: dict[str, object]) -> list[str]:
    raw_reactions = gene.get("sample_reactions") or gene.get("affected_reactions") or []
    tokens: list[str] = []
    seen: set[str] = set()
    for item in raw_reactions:
        token = str(item).strip()
        if not token:
            continue
        token = token.split("_no_", 1)[0]
        token = token.removesuffix("_fwd").removesuffix("_rvs")
        if token and token not in seen:
            seen.add(token)
            tokens.append(token)
    return tokens


def _render_gene_rule_overlay_lookup(query: str) -> None:
    rows = list_pichia_gene_rule_evidence(query)
    if not rows:
        st.info("外部证据 GPR overlay 中未找到匹配。")
        return
    executable_rows = [
        row
        for row in rows
        if str(row.get("rule_status") or "").startswith("overlay_executable")
    ]
    if not executable_rows:
        st.warning("暂无可执行补充规则；这些证据只能用于人工复核或反应级代理解释。")
    st.dataframe(
        pd.DataFrame(
            [
                {
                    "常用名": row.get("common_name") or "",
                    "候选 locus tag": row.get("candidate_locus_tag") or "未解析",
                    "蛋白 / 功能": row.get("protein_name") or "",
                    "证据置信度": row.get("confidence") or "",
                    "GPR 补充状态": _rule_overlay_status_label(row.get("rule_status")),
                    "推荐动作": _rule_overlay_action_label(row.get("recommended_action")),
                    "目标反应": ", ".join(str(item) for item in row.get("target_reaction_ids") or []),
                    "证据来源": _rule_evidence_source_label(row),
                }
                for row in rows
            ]
        ),
        use_container_width=True,
        hide_index=True,
    )


def _rule_evidence_source_label(row: dict[str, object]) -> str:
    if not row:
        return ""
    sources = ", ".join(str(item) for item in row.get("evidence_sources") or [])
    confidence = str(row.get("confidence") or "")
    external_ids = row.get("external_ids") if isinstance(row.get("external_ids"), dict) else {}
    ids = "; ".join(f"{key}: {value}" for key, value in external_ids.items() if value)
    parts = [part for part in (confidence, sources, ids) if part]
    return " | ".join(parts)


def _rule_overlay_status_label(status: object) -> str:
    labels = {
        "overlay_executable": "可执行 overlay（实验性）",
        "overlay_executable_complex_rule": "可执行复合体 overlay（实验性）",
        "display_only_requires_multi_source_confirmation": "仅展示：需要多源确认",
        "display_only_name_context_not_exact_kar2_locus": "仅展示：名称上下文不是精确 locus",
        "display_only_conflicts_with_existing_model_gene_annotation": "仅展示：与当前模型注释冲突",
        "display_only_multiple_or_indirect_candidates": "仅展示：多候选或间接证据",
        "not_executable": "不可执行：证据不足",
    }
    return labels.get(str(status or ""), str(status or ""))


def _rule_overlay_action_label(action: object) -> str:
    labels = {
        "enable_only_for_explicit_analysis": "仅在显式实验模式中使用",
        "keep_reaction_level_proxy_until_locus_is_confirmed": "确认 locus 前保留反应级代理",
        "manual_locus_review_required": "需要人工确认 K. phaffii locus",
        "do_not_replace_existing_model_gene_without_review": "不要在未复核前替换模型基因",
        "manual_review_required": "需要人工复核",
    }
    return labels.get(str(action or ""), str(action or ""))


__all__ = ["render_gene_lookup_panel"]
