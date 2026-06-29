from __future__ import annotations

import pandas as pd
import streamlit as st

from app.services.pichia_gene_catalog_service import (
    get_pichia_ko_genes_for_selection,
    get_pichia_ko_reactions_for_selection,
    get_pichia_oe_reactions_for_selection,
    list_curated_pichia_gene_catalog,
    list_curated_pichia_gene_catalog_by_category,
    load_pichia_full_model_gene_catalog,
)
from app.ui.views.simulation_gene_text import merge_candidate_text


def render_gene_lookup_panel() -> None:
    st.markdown("**毕赤酵母分泌基因库**")
    st.caption(
        "默认显示策略热点；也可以勾选查看全部模型基因。"
        "KO 直接按基因运行，过表达基因会按反应级代理模型运行。"
    )
    show_full = st.checkbox("显示全部模型基因（约 1025 个）", value=False, key="pichia_gene_show_full")
    query = st.text_input(
        "搜索",
        value=st.session_state.get("pichia_gene_lookup_query", ""),
        placeholder="例如：Kar2、ERAD、PAS_chr",
        key="pichia_gene_lookup_query",
    )
    if show_full:
        _render_full_model_gene_lookup(query)
    else:
        _render_curated_gene_lookup(query)


def _render_full_model_gene_lookup(query: str) -> None:
    full_rows = load_pichia_full_model_gene_catalog()
    if query.strip():
        query_text = query.lower()
        full_rows = [
            gene
            for gene in full_rows
            if query_text in gene["gene_id"].lower() or query_text in gene["primary_category"]
        ]
    if not full_rows:
        st.info("未找到匹配。")
        return
    st.dataframe(
        pd.DataFrame(
            [
                {
                    "基因ID": gene["gene_id"],
                    "分类": gene["primary_category"],
                    "通路": gene["processes"][:50],
                    "反应数": gene["n_reactions"],
                }
                for gene in full_rows
            ]
        ),
        use_container_width=True,
        hide_index=True,
    )
    selected = st.multiselect("选择模型基因 ID", [gene["gene_id"] for gene in full_rows], key="pichia_full_sel")
    add_ko_col, add_oe_col = st.columns(2)
    with add_ko_col:
        if st.button("添加到敲除基因") and selected:
            current = str(st.session_state.get("pichia_draft_ko_genes", ""))
            st.session_state["pichia_draft_ko_genes"] = merge_candidate_text(current, selected)
            st.rerun()
    with add_oe_col:
        if st.button("添加到过表达基因代理") and selected:
            current = str(st.session_state.get("pichia_draft_oe_genes", ""))
            st.session_state["pichia_draft_oe_genes"] = merge_candidate_text(current, selected)
            st.rerun()


def _render_curated_gene_lookup(query: str) -> None:
    curated = (
        list_curated_pichia_gene_catalog(query)
        if query.strip()
        else sum((list(value) for value in list_curated_pichia_gene_catalog_by_category().values()), [])
    )
    if not curated:
        st.info("未找到匹配。")
        return
    rows, current_category = [], None
    for entry in curated:
        if entry.category != current_category:
            current_category = entry.category
            rows.append({"基因": f"▸ {current_category}", "描述": "", "ID": "", "方式": ""})
        rows.append(
            {
                "基因": entry.common_name,
                "描述": entry.description,
                "ID": entry.gene_id or entry.oe_reaction_id[:30],
                "方式": entry.intervention.upper(),
            }
        )
    st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)
    selected = st.multiselect("选择基因", [entry.common_name for entry in curated if entry.common_name], key="pichia_gene_sel")
    if st.button("添加到敲除输入") and selected:
        genes = get_pichia_ko_genes_for_selection(selected)
        if genes:
            current = str(st.session_state.get("pichia_draft_ko_genes", ""))
            st.session_state["pichia_draft_ko_genes"] = merge_candidate_text(current, genes)
        else:
            reactions = get_pichia_ko_reactions_for_selection(selected)
            if reactions:
                current = str(st.session_state.get("pichia_draft_ko_reactions", ""))
                st.session_state["pichia_draft_ko_reactions"] = merge_candidate_text(current, reactions)
            else:
                st.toast("所选条目无对应的敲除模型基因或反应 ID")
        st.rerun()
    if st.button("添加到过表达反应代理") and selected:
        reactions = get_pichia_oe_reactions_for_selection(selected)
        if reactions:
            current = str(st.session_state.get("pichia_draft_oe_reactions", ""))
            st.session_state["pichia_draft_oe_reactions"] = merge_candidate_text(current, reactions)
        else:
            st.toast("所选条目无对应的过表达反应 ID（可能只支持敲除）")
        st.rerun()


__all__ = ["render_gene_lookup_panel"]
