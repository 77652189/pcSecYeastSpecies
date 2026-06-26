from __future__ import annotations

from dataclasses import dataclass

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
from app.services.pichia_screen_preview_service import preview_screen_inputs
from app.services.pichia_secretion_schema import SecretionRunRequest
from app.ui.common import PATHS


@dataclass(frozen=True)
class GenePerturbationFormState:
    ko_gene_text: str
    ko_reaction_text: str
    oe_gene_text: str
    oe_reaction_text: str
    candidate_limit: int

    @property
    def ko_gene_ids(self) -> tuple[str, ...]:
        return parse_candidate_text(self.ko_gene_text)

    @property
    def ko_reaction_ids(self) -> tuple[str, ...]:
        return parse_candidate_text(self.ko_reaction_text)

    @property
    def oe_gene_ids(self) -> tuple[str, ...]:
        return parse_candidate_text(self.oe_gene_text)

    @property
    def oe_reaction_ids(self) -> tuple[str, ...]:
        return parse_candidate_text(self.oe_reaction_text)


def parse_candidate_text(text: str) -> tuple[str, ...]:
    values: list[str] = []
    for raw in str(text or "").replace("，", ",").replace(";", ",").splitlines():
        for item in raw.split(","):
            cleaned = item.strip()
            if cleaned and cleaned not in values:
                values.append(cleaned)
    return tuple(values)


def render_gene_lookup_panel() -> None:
    st.markdown("**毕赤酵母分泌基因库**")
    st.caption("默认显示 37 条策展靶点；也可以勾选查看全部模型基因。KO 直接按基因运行，过表达基因会按反应级代理模型运行。")
    show_full = st.checkbox("显示全部模型基因（~1025 个）", value=False, key="pichia_gene_show_full")
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


def _merge_candidate_text(existing: str, additions: list[str]) -> str:
    values = list(parse_candidate_text(existing))
    for item in additions:
        if item and item not in values:
            values.append(item)
    return "\n".join(values)


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
            st.session_state["pichia_draft_ko_genes"] = _merge_candidate_text(current, selected)
            st.rerun()
    with add_oe_col:
        if st.button("添加到过表达基因代理") and selected:
            current = str(st.session_state.get("pichia_draft_oe_genes", ""))
            st.session_state["pichia_draft_oe_genes"] = _merge_candidate_text(current, selected)
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
            rows.append({"基因": f"⬇ {current_category}", "描述": "", "ID": "", "方式": ""})
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
            st.session_state["pichia_draft_ko_genes"] = _merge_candidate_text(current, genes)
        else:
            reactions = get_pichia_ko_reactions_for_selection(selected)
            if reactions:
                current = str(st.session_state.get("pichia_draft_ko_reactions", ""))
                st.session_state["pichia_draft_ko_reactions"] = _merge_candidate_text(current, reactions)
            else:
                st.toast("所选条目无对应的敲除模型基因或反应 ID", icon="⚠️")
        st.rerun()
    if st.button("添加到过表达反应代理") and selected:
        reactions = get_pichia_oe_reactions_for_selection(selected)
        if reactions:
            current = str(st.session_state.get("pichia_draft_oe_reactions", ""))
            st.session_state["pichia_draft_oe_reactions"] = _merge_candidate_text(current, reactions)
        else:
            st.toast("所选条目无对应的过表达反应 ID（可能只支持敲除）", icon="⚠️")
        st.rerun()


def render_gene_perturbation_form(target_id: str) -> GenePerturbationFormState:
    with st.expander("基因扰动", expanded=True):
        st.caption("可以从基因库选择，也可以手动输入。多个条目用逗号或换行分隔；单次最多 20 个候选。")
        render_gene_lookup_panel()
        with st.expander("输入说明", expanded=False):
            st.markdown(
                "- 敲除基因：填写模型基因 ID，例如 `PAS_chr2-2_0107`。\n"
                "- 敲除反应：用于没有明确基因 ID 的复合体级扰动，例如 `sec_Och1p_complex_formation`。\n"
                "- 过表达基因：会先解析到该基因参与的反应，再按 reaction-level OE proxy 运行。\n"
                "- 过表达反应：直接填写模型反应 ID，适合从上方策展库自动填入。"
            )
        ko_text = st.text_area(
            "敲除基因（KO gene）",
            height=60,
            key="pichia_draft_ko_genes",
            placeholder="例如：PAS_chr2-2_0107",
        )
        ko_rxn = st.text_area(
            "敲除反应（复合体级 KO）",
            height=60,
            key="pichia_draft_ko_reactions",
            placeholder="例如：sec_Och1p_complex_formation",
        )
        oe_text = st.text_area(
            "过表达基因（OE gene proxy）",
            height=60,
            key="pichia_draft_oe_genes",
            placeholder="例如：PAS_chr1-4_0586；会解析为反应级过表达代理",
        )
        oe_rxn = st.text_area(
            "过表达反应（高级 / OE reaction）",
            height=60,
            key="pichia_draft_oe_reactions",
            placeholder="例如：sec_Kar2p_complex_formation",
        )
        limit = int(st.number_input("候选数上限", 1, 20, 20, 1, key="pichia_limit"))
        state = GenePerturbationFormState(
            ko_gene_text=ko_text,
            ko_reaction_text=ko_rxn,
            oe_gene_text=oe_text,
            oe_reaction_text=oe_rxn,
            candidate_limit=limit,
        )
        if st.button("预检 KO/OE 输入", key="pichia_preview_screen_inputs"):
            _render_screen_input_preview(target_id, state)
        return state


def _render_screen_input_preview(target_id: str, state: GenePerturbationFormState) -> None:
    preview_request = SecretionRunRequest(
        target_source="builtin",
        target_id=target_id,
        ko_gene_ids=state.ko_gene_ids,
        ko_reaction_ids=state.ko_reaction_ids,
        oe_gene_ids=state.oe_gene_ids,
        oe_reaction_ids=state.oe_reaction_ids,
        screen_candidate_limit=state.candidate_limit,
    )
    with st.spinner("正在解析模型基因和反应 ID……"):
        preview = preview_screen_inputs(preview_request, PATHS)
    if preview.get("warnings"):
        for warning in preview["warnings"]:
            st.warning(warning)
    preview_rows = []
    for group_label, key in (
        ("敲除基因", "ko_genes"),
        ("敲除反应", "ko_reactions"),
        ("过表达基因代理", "oe_genes"),
        ("过表达反应", "oe_reactions"),
    ):
        for row in preview.get(key, []):
            preview_rows.append(
                {
                    "类别": group_label,
                    "输入": row.get("input_id"),
                    "状态": "已解析" if row.get("resolved") else "未解析",
                    "解析到的反应数": row.get("resolved_reaction_count"),
                    "反应预览": ", ".join(row.get("resolved_reactions_preview") or []),
                }
            )
    if preview_rows:
        st.dataframe(pd.DataFrame(preview_rows), use_container_width=True, hide_index=True)
    else:
        st.info("当前没有手动 KO/OE 输入；正式运行会使用小规模默认 smoke 候选。")


__all__ = [
    "GenePerturbationFormState",
    "parse_candidate_text",
    "render_gene_lookup_panel",
    "render_gene_perturbation_form",
]
