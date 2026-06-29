from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import pandas as pd
import streamlit as st

from app.services.pichia_screen_preview_service import preview_screen_inputs
from app.services.pichia_secretion_schema import SecretionRunRequest
from app.ui.common import PATHS
from app.ui.views.simulation_gene_catalog import render_gene_lookup_panel
from app.ui.views.simulation_gene_text import parse_candidate_text


MAPPING_LEVEL_LABELS = {
    "direct_gpr": "GPR 直接关联",
    "complex_subunit": "复合体亚基",
    "reaction_proxy": "反应代理",
    "metabolic_or_other": "代谢/其它",
    "unresolved": "未解析",
}

MAPPING_CONFIDENCE_LABELS = {
    "high": "高",
    "medium": "中",
    "low": "低",
    "unresolved": "未解析",
}


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


def render_gene_perturbation_form(target_id: str) -> GenePerturbationFormState:
    with st.expander("基因扰动", expanded=True):
        st.caption("可以从基因库选择，也可以手动输入。多个条目用逗号、分号或换行分隔；单次最多 20 个候选。")
        render_gene_lookup_panel()
        with st.expander("输入说明", expanded=False):
            st.markdown(
                "- 敲除基因：填写模型基因 ID，例如 `PAS_chr2-2_0107`。\n"
                "- 敲除反应：用于没有明确基因 ID 的复合体级扰动，例如 `sec_Och1p_complex_formation`。\n"
                "- 过表达基因：会先解析到该基因参与的反应，再按 reaction-level OE proxy 运行。\n"
                "- 过表达反应：直接填写模型反应 ID，适合从上方策略库自动填入。"
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


def gene_mapping_rows_for_display(rows: list[dict[str, Any]]) -> pd.DataFrame:
    display_rows: list[dict[str, Any]] = []
    for row in rows:
        display_rows.append(
            {
                "基因": row.get("gene_id", ""),
                "反应": row.get("reaction_id") or "未解析",
                "分泌环节": row.get("secretory_process", ""),
                "映射层级": MAPPING_LEVEL_LABELS.get(
                    str(row.get("mapping_level", "")),
                    row.get("mapping_level", ""),
                ),
                "置信度": MAPPING_CONFIDENCE_LABELS.get(
                    str(row.get("mapping_confidence", "")),
                    row.get("mapping_confidence", ""),
                ),
                "复合体": row.get("complex_id") or "",
                "复合体亚基": _join_values(row.get("complex_subunit_ids")),
                "亚基计量": _join_values(row.get("complex_subunit_stoichiometry")),
                "解释": row.get("interpretation", ""),
            }
        )
    return pd.DataFrame(display_rows)


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
    with st.spinner("正在解析模型基因和反应 ID..."):
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

    mapping_rows = list(preview.get("gene_mapping_rows") or [])
    if mapping_rows:
        st.markdown("**基因影响路径预览**")
        st.caption("该表只解释 gene -> reaction -> 分泌环节映射，不运行求解器，也不改变 KO/OE 数值结果。")
        st.dataframe(gene_mapping_rows_for_display(mapping_rows), use_container_width=True, hide_index=True)


def _join_values(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, (list, tuple)):
        return ", ".join(str(item) for item in value)
    return str(value)


__all__ = [
    "GenePerturbationFormState",
    "gene_mapping_rows_for_display",
    "parse_candidate_text",
    "render_gene_perturbation_form",
]
