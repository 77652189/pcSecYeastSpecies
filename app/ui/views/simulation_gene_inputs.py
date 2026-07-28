from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import pandas as pd
import streamlit as st

from app.services.pichia_screen_preview_service import preview_screen_inputs
from app.services.pichia_secretion_schema import SecretionRunRequest
from app.ui.common import PATHS
from app.ui.views.hlf_opn_candidate_panel import render_hlf_opn_candidate_panel
from app.ui.views.simulation_display import (
    GPR_ROLE_LABELS,
    KO_SUPPORT_STATUS_LABELS,
    OE_SUPPORT_STATUS_LABELS,
)
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
    enable_gene_rule_overlay: bool = False

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


def _has_text(session_key: str) -> bool:
    return bool(str(st.session_state.get(session_key) or "").strip())


def render_gene_perturbation_form(target_id: str) -> GenePerturbationFormState:
    """基因扰动输入。

    分层原则（2026-07-28 重排）：主输入（KO/OE **基因**）紧跟在页面顶部，反应级扰动降为"高级"、
    找基因的帮手收进一行折叠入口。此前是三个展开的辅助面板 + 长说明挡在输入框前面、四个框平权
    并列，研究员看不出该填哪个。

    **帮手面板必须留在输入框之前**（只是折叠起来）：它的"加入 KO/OE 输入"按钮会写
    `st.session_state["pichia_draft_ko_genes"]`，而 Streamlit 不允许在同名 key 的控件实例化之后
    再改它的 session_state。把它挪到输入框下面会让那个按钮报错——`test_hlf_opn_candidate_panel_
    is_embedded_before_gene_textareas` 锁的就是这条，别"顺手"调换顺序。

    session_state 的 key 一律不变——`apply_simulation_prefill` 靠它们从筛查页填入候选。
    """
    with st.expander("基因扰动", expanded=True):
        # 反应级候选（策展库 / 复合体假设）会被预填进下面的高级区；若真被填了就必须自动展开，
        # 否则用户从筛查页点"在仿真验证中核实"跳来只看到两个空的基因框，会以为跳转没生效。
        reaction_prefilled = _has_text("pichia_draft_ko_reactions") or _has_text("pichia_draft_oe_reactions")

        # 折叠成一行入口，展开才出内容——此前这两个面板默认展开、把输入框挤到屏幕外。
        # 位置必须在下面的输入框之前，原因见本函数 docstring（"加入 KO/OE 输入"要写 session_state）。
        with st.expander("不知道基因 ID？在这里查找 / 从候选库挑", expanded=False):
            render_hlf_opn_candidate_panel(target_id)
            render_gene_lookup_panel()

        st.caption(
            "最常用：填要敲除或过表达的**模型基因 ID**。从「全基因组KO/OE筛查」点“在仿真验证中核实”"
            "会自动填好这里，不用手抄。多个条目用逗号、分号或换行分隔；单次最多 20 个候选。"
        )
        ko_text = st.text_area(
            "敲除基因（KO gene）",
            height=60,
            key="pichia_draft_ko_genes",
            placeholder="例如：PAS_chr2-2_0107",
            help="正式基因级 KO：系统按 GPR 规则关闭会失活的反应。",
        )
        oe_text = st.text_area(
            "过表达基因（OE gene proxy）",
            height=60,
            key="pichia_draft_oe_genes",
            placeholder="例如：PAS_chr1-4_0586；会解析为反应级过表达代理",
            help=(
                "先做 GPR-aware 规划，只有单基因 / 同工酶等可解释场景才运行 reaction-level OE proxy；"
                "复合体亚基默认不做单基因 OE（没有亚基化学计量证据时不虚构容量提升）。"
            ),
        )

        with st.expander(
            "复合体 / 反应级扰动（高级）" + ("　←　已从筛查页填入，展开可见" if reaction_prefilled else ""),
            expanded=reaction_prefilled,
        ):
            st.caption(
                "分泌机器（chaperone、糖基化、COPII 等复合体）在模型里**没有 GPR**，只能在反应层面扰动，"
                "所以这类候选填的是反应 ID 而不是基因 ID。筛查页的复合体级候选会自动填到这里。"
            )
            ko_rxn = st.text_area(
                "敲除反应（复合体级 KO）",
                height=60,
                key="pichia_draft_ko_reactions",
                placeholder="例如：sec_Och1p_complex_formation",
                help="用于没有明确基因 ID 的复合体级扰动：直接把该复合体形成反应的流量压到 0。",
            )
            oe_rxn = st.text_area(
                "过表达反应（高级 / OE reaction）",
                height=60,
                key="pichia_draft_oe_reactions",
                placeholder="例如：sec_Kar2p_complex_formation",
                help="高级诊断入口：直接填模型反应 ID，把其 kcat 乘以倍数；不代表完整基因表达调控模型。",
            )
            enable_overlay = st.checkbox(
                "使用外部证据补充 GPR（实验性，默认关闭）",
                value=False,
                key="pichia_enable_gene_rule_overlay",
                help=(
                    "只在 KO/OE 预检和显式运行中使用 Python 侧证据 overlay；"
                    "不会写回原始模型，也不是 MATLAB 原始 GPR。"
                ),
            )
            limit = int(st.number_input("候选数上限", 1, 20, 20, 1, key="pichia_limit"))

        state = GenePerturbationFormState(
            ko_gene_text=ko_text,
            ko_reaction_text=ko_rxn,
            oe_gene_text=oe_text,
            oe_reaction_text=oe_rxn,
            candidate_limit=limit,
            enable_gene_rule_overlay=enable_overlay,
        )
        # 预检按钮紧贴它自己的输出（结果渲染在下方），故放在各输入区之后。
        if st.button("预检 KO/OE 输入", key="pichia_preview_screen_inputs"):
            _render_screen_input_preview(target_id, state)
        return state


def gene_mapping_rows_for_display(rows: list[dict[str, Any]]) -> pd.DataFrame:
    display_rows: list[dict[str, Any]] = []
    for row in rows:
        display_rows.append(
            {
                "基因": row.get("input_gene_id") or row.get("gene_id", ""),
                "模型基因": row.get("canonical_gene_id") or row.get("gene_id", ""),
                "显示名称": row.get("gene_display_name") or "",
                "标准符号": row.get("standard_symbol") or "",
                "蛋白名称": row.get("protein_name") or "",
                "命名置信度": row.get("annotation_confidence") or "",
                "命名状态": row.get("standard_name_status") or "",
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
        enable_gene_rule_overlay=state.enable_gene_rule_overlay,
    )
    with st.spinner("正在解析模型基因和反应 ID..."):
        preview = preview_screen_inputs(preview_request, PATHS)
    if preview.get("warnings"):
        for warning in preview["warnings"]:
            st.warning(warning)
    overlay = preview.get("gene_rule_overlay") if isinstance(preview.get("gene_rule_overlay"), dict) else {}
    if overlay:
        st.caption(
            "外部证据 GPR overlay：实验性 Python 分析副本；"
            f"可执行补充规则 {overlay.get('entry_count', 0)} 条，不写回原始模型。"
        )

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
                    "显示名称": row.get("gene_display_name") or "",
                    "标准符号": row.get("standard_symbol") or "",
                    "蛋白名称": row.get("protein_name") or "",
                    "命名状态": row.get("standard_name_status") or "",
                    "状态": "已解析" if row.get("resolved") else "未解析",
                    "解析到的反应数": row.get("resolved_reaction_count"),
                    "反应预览": ", ".join(row.get("resolved_reactions_preview") or []),
                    "KO 支持": KO_SUPPORT_STATUS_LABELS.get(
                        row.get("ko_support_status") or "", row.get("ko_support_status") or ""
                    ),
                    "OE 支持": OE_SUPPORT_STATUS_LABELS.get(
                        row.get("oe_support_status") or "", row.get("oe_support_status") or ""
                    ),
                    "GPR 角色": GPR_ROLE_LABELS.get(row.get("gpr_role") or "", row.get("gpr_role") or ""),
                    "置信度": MAPPING_CONFIDENCE_LABELS.get(
                        row.get("confidence") or row.get("mapping_confidence") or "",
                        row.get("confidence") or row.get("mapping_confidence") or "",
                    ),
                    "缺失信息": ", ".join(str(item) for item in row.get("missing_information") or []),
                }
            )
    if preview_rows:
        st.dataframe(pd.DataFrame(preview_rows), width='stretch', hide_index=True)
    else:
        st.info("当前没有手动 KO/OE 输入；正式运行会使用小规模默认 smoke 候选。")

    mapping_rows = list(preview.get("gene_mapping_rows") or [])
    if mapping_rows:
        st.markdown("**基因影响路径预览**")
        st.caption("该表只解释 gene -> reaction -> 分泌环节映射，不运行求解器，也不改变 KO/OE 数值结果。")
        st.dataframe(gene_mapping_rows_for_display(mapping_rows), width='stretch', hide_index=True)


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
