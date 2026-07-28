"""E2 统一候选选择器（[ADR-007](../../../docs/adr/007-secretory-machinery-gene-complex-reachability.md)）。

**为什么需要它**：模型内部有两条互不相交的候选路径——基因级（走 GPR、只覆盖代谢）和复合体/反应级
（分泌机器 2793 个复合体形成反应、零基因关联）。此前界面把这条分裂直接摊给用户：四个输入框
（KO基因/KO反应/OE基因/OE反应）平权并列，用户得先知道"我要动的东西在模型里算基因还是算反应"
才能填对框。但研究员的心智是"我要动 PDI1"——是基因还是反应属**模型实现细节**，不该由用户承担。

本模块把两者统一成"一行一个改造对象"：勾选 → 系统按对象类型路由到正确的输入框。

分层：纯数据整形（`build_unified_candidate_rows` / `route_selected_candidates`）与渲染分开，前者可测；
判定哪些候选算"可执行 vs 待复核"仍归引擎（`executable_inputs_for_hlf_opn_candidates`），本模块只搬运。
"""

from __future__ import annotations

from typing import Any

import pandas as pd
import streamlit as st

from app.services.pichia_gene_catalog_service import load_hlf_opn_candidate_genes
from app.ui.views.simulation_gene_text import merge_candidate_text


KIND_GENE = "基因"
KIND_COMPLEX = "复合体"
ALL_PROCESSES = "全部环节"

# 勾选后各类对象要写进哪个输入框。与 simulation.py 的 _prefill_field_values 同一套路由规则：
# 基因走基因框（经 GPR 解析），复合体走反应框（直接是模型反应 id、无基因可解析）。
_FIELD_BY_ROUTE = {
    ("KO", KIND_GENE): "pichia_draft_ko_genes",
    ("OE", KIND_GENE): "pichia_draft_oe_genes",
    ("KO", KIND_COMPLEX): "pichia_draft_ko_reactions",
    ("OE", KIND_COMPLEX): "pichia_draft_oe_reactions",
}

_GENE_EXECUTABLE_STATUSES = frozenset({"model_ko_executable", "model_oe_proxy_executable"})


def build_unified_candidate_rows(candidates: list[dict[str, Any]]) -> pd.DataFrame:
    """把策展候选整形成"一行一个可勾选的改造对象"，基因与复合体并列。

    每行的 `作用对象` 决定路由：能按基因跑的走基因，否则用它的复合体反应。两者都没有的行不出现在
    可选择列表里（选了也跑不了）；`review_only` 同理剔除——它没有建议的扰动方向。
    """
    rows: list[dict[str, Any]] = []
    for candidate in candidates:
        intervention = str(candidate.get("recommended_intervention") or "").upper()
        if intervention not in {"KO", "OE"}:
            continue
        gene_id = str(candidate.get("gene_id") or "").strip()
        status = str(candidate.get("operability_status") or "")
        reactions = [
            str(item).strip()
            for item in (
                candidate.get("executable_ko_reactions")
                or candidate.get("executable_oe_proxy_reactions")
                or candidate.get("review_reactions")
                or []
            )
            if str(item).strip()
        ]
        if status in _GENE_EXECUTABLE_STATUSES and gene_id:
            kind, model_object = KIND_GENE, gene_id
            confidence = "基因可直接跑"
        elif reactions:
            # 复合体反应在模型里可跑；不确定的是"实验室该动哪几个基因"（同源证据待复核）。
            kind, model_object = KIND_COMPLEX, reactions[0]
            confidence = "反应可跑·基因归属待复核"
        else:
            continue
        common_name = str(candidate.get("source_common_name") or "").strip()
        rows.append(
            {
                "选择": False,
                "候选": common_name or model_object,
                "分泌环节": str(candidate.get("source_category") or "未分类"),
                "改造方式": intervention,
                "作用对象": kind,
                "模型对象": model_object,
                "把握": confidence,
            }
        )
    frame = pd.DataFrame(rows)
    if frame.empty:
        return frame
    # 同一个复合体反应可能被多条策展条目指向（如 PDI1/ERO1/ERV2 共用一个复合体）——去重成一行，
    # 否则用户会看到几行长得一样、勾哪个都一样的候选。
    frame = frame.drop_duplicates(subset=["改造方式", "模型对象"], keep="first").reset_index(drop=True)
    return frame.sort_values(["分泌环节", "作用对象", "候选"]).reset_index(drop=True)


def route_selected_candidates(selected: pd.DataFrame) -> dict[str, list[str]]:
    """把勾选的行按 (改造方式 × 作用对象) 分到四个输入框，保持勾选顺序、去重。"""
    routed: dict[str, list[str]] = {field: [] for field in _FIELD_BY_ROUTE.values()}
    if selected is None or selected.empty:
        return routed
    for _, row in selected.iterrows():
        field = _FIELD_BY_ROUTE.get((str(row.get("改造方式") or "").upper(), str(row.get("作用对象") or "")))
        model_object = str(row.get("模型对象") or "").strip()
        if field is None or not model_object:
            continue
        if model_object not in routed[field]:
            routed[field].append(model_object)
    return routed


def apply_routed_candidates(routed: dict[str, list[str]]) -> int:
    """写进 session_state（merge、不覆盖用户已填内容），返回实际加入的条目数。

    必须在四个 text_area 控件实例化**之前**调用——Streamlit 不允许在同名 key 的控件创建后再改它的
    session_state（见 simulation_gene_inputs.render_gene_perturbation_form 的 docstring）。
    """
    added = 0
    for field, items in routed.items():
        if not items:
            continue
        st.session_state[field] = merge_candidate_text(str(st.session_state.get(field, "")), items)
        added += len(items)
    return added


def _attach_lab_gene_hints(frame: pd.DataFrame) -> pd.DataFrame:
    """给复合体行补「实验时对应基因」（E3 反向映射，ADR-007）。

    无策展映射数据时整列为 "—"——这是**预期状态**（策展待拍板），不是错误；此时研究员仍能跑复合体
    反应，只是"该动哪几个基因"要自己查。策展数据放进 Data/pcSecPichia/gene_complex_mapping.json 即生效。
    """
    if frame.empty:
        return frame
    try:
        from app.services.gene_complex_mapping_service import lab_gene_hint_for_complex
    except Exception:  # noqa: BLE001 - 映射层缺席不该拖垮选择器
        return frame
    hints: list[str] = []
    for _, row in frame.iterrows():
        if str(row.get("作用对象")) != KIND_COMPLEX:
            hints.append("—")  # 基因行本身就是基因，无需翻译
            continue
        try:
            hint = lab_gene_hint_for_complex(str(row.get("模型对象") or ""))
        except Exception:  # noqa: BLE001
            hint = ""
        hints.append(hint or "—")
    frame = frame.copy()
    frame["实验时对应基因"] = hints
    return frame


def render_candidate_selector(target_id: str, *, target_context: str) -> None:
    """勾选式候选选择器。放在输入框之前（见 apply_routed_candidates 的 session_state 约束）。"""
    candidates = load_hlf_opn_candidate_genes(target_context=target_context, include_shared=True)
    frame = build_unified_candidate_rows(candidates)
    if frame.empty:
        st.caption("当前目标没有可直接选择的策展候选。")
        return

    st.caption(
        "勾选想改造的对象 → 点下面的按钮，系统会自动填进对应输入框（**基因和复合体都在这里，"
        "不用自己判断该填哪个框**）。「把握」列说明模型能给到什么程度。"
    )
    frame = _attach_lab_gene_hints(frame)
    processes = [ALL_PROCESSES, *sorted(frame["分泌环节"].unique())]
    chosen_process = st.selectbox(
        "按分泌环节筛选",
        processes,
        key=f"candidate_selector_process_{target_context}",
        help="折叠/糖基化/转运等环节。想解决折叠瓶颈就看「二硫键」「ER 折叠与分子伴侣」。",
    )
    view_frame = frame if chosen_process == ALL_PROCESSES else frame[frame["分泌环节"] == chosen_process]

    edited = st.data_editor(
        view_frame,
        width="stretch",
        hide_index=True,
        disabled=[column for column in view_frame.columns if column != "选择"],
        column_config={
            "选择": st.column_config.CheckboxColumn("选择", help="勾选后点下面的按钮加入 KO/OE 输入"),
            "模型对象": st.column_config.TextColumn("模型对象（实际算的东西）"),
            "实验时对应基因": st.column_config.TextColumn(
                "实验时对应基因",
                help=(
                    "复合体在模型里没有基因关联，这一列来自人工策展的「基因↔复合体」映射（ADR-007）。"
                    "显示 — 表示还没有策展映射，需自行确认该复合体对应哪几个基因。"
                ),
            ),
        },
        key=f"candidate_selector_editor_{target_context}_{chosen_process}",
    )
    selected = edited[edited["选择"]] if "选择" in edited.columns else pd.DataFrame()

    if st.button(
        f"加入选中的 {len(selected)} 个候选到 KO/OE 输入",
        key=f"candidate_selector_apply_{target_context}",
        disabled=selected.empty,
        type="primary",
    ):
        routed = route_selected_candidates(selected)
        added = apply_routed_candidates(routed)
        genes = len(routed["pichia_draft_ko_genes"]) + len(routed["pichia_draft_oe_genes"])
        reactions = len(routed["pichia_draft_ko_reactions"]) + len(routed["pichia_draft_oe_reactions"])
        st.toast(f"已加入 {added} 个候选：基因 {genes} 个、复合体反应 {reactions} 个。")
        st.rerun()  # 让下面的输入框重新按新 session_state 渲染


__all__ = [
    "ALL_PROCESSES",
    "KIND_COMPLEX",
    "KIND_GENE",
    "apply_routed_candidates",
    "build_unified_candidate_rows",
    "render_candidate_selector",
    "route_selected_candidates",
]
