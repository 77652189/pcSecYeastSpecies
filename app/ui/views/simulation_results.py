from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import streamlit as st

from app.services.pichia_background_tasks import (
    load_latest_completed_background_result,
    load_last_result,
    poll_background_simulation,
    save_last_result,
)
from app.ui.common import PATHS
from app.ui.views.candidate_path_graph import render_secretion_path_graph
from app.ui.views.simulation_display import (
    CANDIDATE_DISPLAY_COLUMNS,
    candidate_effect_counts,
    candidate_row_label,
    normalise_candidate_frame_for_display,
)


def render_pichia_results() -> None:
    tsp = st.session_state.get("pichia_draft_task_status_path")
    if tsp:
        status, msg, result = poll_background_simulation(tsp)
        if status == "done" and result:
            st.session_state["last_pichia_secretion_draft_response"] = result
            save_last_result(result, PATHS)
            st.session_state.pop("pichia_draft_task_status_path", None)
            st.session_state.pop("pichia_draft_task_id", None)
            st.session_state["pichia_switch_to_results"] = True
            st.rerun()
        elif status == "error":
            st.error(f"仿真失败：{msg}")
            st.session_state.pop("pichia_draft_task_status_path", None)
            st.session_state.pop("pichia_draft_task_id", None)
        elif status in ("pending", "running"):
            st.info(f"⏳ {msg}")
            st.caption("仿真在后台继续执行。点击下方按钮刷新状态；完成后会显示结果摘要和输出文件。")
            if st.button("刷新仿真状态", key="pichia_refresh_task_status_button"):
                st.rerun()
            return
        elif status in ("lost", "stale"):
            latest = load_latest_completed_background_result(PATHS)
            if latest:
                st.warning(f"{msg} 已找到最近完成的仿真结果。")
                st.session_state["last_pichia_secretion_draft_response"] = latest
                save_last_result(latest, PATHS)
                st.session_state.pop("pichia_draft_task_status_path", None)
                st.session_state.pop("pichia_draft_task_id", None)
                st.session_state["pichia_switch_to_results"] = True
                st.rerun()
            st.error(f"任务状态异常：{msg}")
            st.caption("可以回到仿真构建页重新提交；旧任务不会继续阻塞页面。")
            st.session_state.pop("pichia_draft_task_status_path", None)
            st.session_state.pop("pichia_draft_task_id", None)

    if "last_pichia_secretion_draft_response" not in st.session_state:
        cached = load_last_result(PATHS)
        if not cached:
            cached = load_latest_completed_background_result(PATHS)
            if cached:
                save_last_result(cached, PATHS)
        if cached:
            st.session_state["last_pichia_secretion_draft_response"] = cached

    data = st.session_state.get("last_pichia_secretion_draft_response")
    if not data:
        st.info("👈 在「仿真构建」页面选择目标并运行。")
        return

    st.success("✅ 仿真完成") if data.get("success") else st.error("❌ 仿真失败")
    c1, c2, c3 = st.columns([2, 1, 1])
    with c1:
        st.metric("📦 分泌通量", data.get("objective_value", "—"), help="相对比较值，非实际产量")
    with c2:
        st.metric("目标蛋白", data.get("target_id", "—"))
    with c3:
        st.metric("MATLAB 对齐", data.get("matlab_alignment_status", "—"))
    st.caption("不同构建之间可横向对比。不代表实际发酵产量。")

    with st.expander("参数", expanded=False):
        medium_condition = data.get("medium_condition") if isinstance(data.get("medium_condition"), dict) else {}
        st.dataframe(
            pd.DataFrame(
                [
                    {"参数": key, "值": value}
                    for key, value in {
                        "目标": data.get("target_id"),
                        "状态": data.get("result_status"),
                        "MATLAB": data.get("matlab_alignment_status"),
                        "目标值": data.get("objective_value"),
                        "培养基条件": medium_condition.get("condition_id") if isinstance(medium_condition, dict) else None,
                        "碳源": medium_condition.get("carbon_source_id") if isinstance(medium_condition, dict) else None,
                        "培养基状态": medium_condition.get("status") if isinstance(medium_condition, dict) else None,
                        "科学解释状态": medium_condition.get("scientific_status") if isinstance(medium_condition, dict) else None,
                    }.items()
                    if value
                ]
            ),
            use_container_width=True,
            hide_index=True,
        )

    files = {
        "摘要": data.get("summary_path"),
        "报告": data.get("report_path"),
        "候选表": data.get("candidate_table_path"),
        "权衡": data.get("tradeoff_path"),
    }
    st.write("**输出文件**")
    st.dataframe(
        pd.DataFrame([{"文件": key, "路径": value} for key, value in files.items() if value]),
        use_container_width=True,
        hide_index=True,
    )

    warns = data.get("warnings") or []
    if warns:
        with st.expander("注意事项", expanded=False):
            for warning in warns:
                st.warning(warning)
    errors = data.get("errors") or []
    if errors:
        with st.expander("错误", expanded=True):
            for error in errors:
                st.error(error)

    protein_cost = _protein_cost_payload(data)
    if protein_cost:
        _render_protein_cost_analysis(protein_cost)
    target_growth = _target_growth_payload(data)
    if target_growth:
        _render_target_growth_analysis(target_growth)
    yield_recommendations = _yield_recommendation_payload(data)
    if yield_recommendations:
        _render_yield_improvement_recommendations(yield_recommendations)

    candidate_path = data.get("candidate_table_path")
    if candidate_path and Path(candidate_path).exists():
        with st.expander("候选表与分泌路径", expanded=True):
            _render_candidate_outputs(str(candidate_path), data.get("summary_path"))
    tradeoff_path = data.get("tradeoff_path")
    if tradeoff_path and Path(tradeoff_path).exists():
        with st.expander("生长权衡", expanded=False):
            st.dataframe(pd.read_csv(tradeoff_path), use_container_width=True)
    with st.expander("原始响应", expanded=False):
        st.caption("调试用")
        st.json(data)


def _render_candidate_outputs(candidate_path: str, summary_path: str | None) -> None:
    try:
        frame = pd.read_csv(candidate_path)
    except Exception as exc:
        st.warning(f"候选表读取失败：{exc}")
        return
    if frame.empty:
        st.info("候选表为空。")
        return
    frame = normalise_candidate_frame_for_display(frame)
    if "delta_objective" in frame.columns:
        frame = frame.sort_values("delta_objective", ascending=False, na_position="last")
    effect_counts = candidate_effect_counts(frame)
    if effect_counts:
        st.markdown("**候选分类汇总**")
        cols = st.columns(len(effect_counts))
        for col, (label, value) in zip(cols, effect_counts.items()):
            col.metric(label, value)
        if effect_counts.get("约束不可行"):
            st.info("「约束不可行」表示当前固定生长速率和约束组合下没有可行解，不等同于真实发酵条件必然不可行。")

    display = frame.rename(columns={key: value for key, value in CANDIDATE_DISPLAY_COLUMNS.items() if key in frame.columns})
    st.dataframe(display, use_container_width=True, hide_index=True)

    if len(frame) > 0:
        selectable = frame.reset_index(drop=True)
        row_labels = [candidate_row_label(idx, row) for idx, row in selectable.iterrows()]
        selected_index = st.selectbox("选择一行查看分泌路径图", range(len(row_labels)), format_func=lambda i: row_labels[i])
        row = selectable.iloc[int(selected_index)]
        summary = None
        if summary_path and Path(summary_path).exists():
            try:
                summary = json.loads(Path(summary_path).read_text(encoding="utf-8"))
            except Exception:
                summary = {}
        render_secretion_path_graph(row.to_dict(), summary or {})


def _protein_cost_payload(data: dict[str, object]) -> dict[str, object]:
    payload = data.get("protein_cost_analysis")
    if isinstance(payload, dict) and payload:
        return payload
    summary_path = data.get("summary_path")
    if not summary_path or not Path(str(summary_path)).exists():
        return {}
    try:
        summary = json.loads(Path(str(summary_path)).read_text(encoding="utf-8"))
    except Exception:
        return {}
    payload = summary.get("protein_cost_analysis") if isinstance(summary, dict) else None
    return payload if isinstance(payload, dict) else {}


def _target_growth_payload(data: dict[str, object]) -> dict[str, object]:
    payload = data.get("target_growth_analysis")
    if isinstance(payload, dict) and payload:
        return payload
    summary_path = data.get("summary_path")
    if not summary_path or not Path(str(summary_path)).exists():
        return {}
    try:
        summary = json.loads(Path(str(summary_path)).read_text(encoding="utf-8"))
    except Exception:
        return {}
    payload = summary.get("target_growth_analysis") if isinstance(summary, dict) else None
    return payload if isinstance(payload, dict) else {}


def _yield_recommendation_payload(data: dict[str, object]) -> dict[str, object]:
    payload = data.get("yield_improvement_recommendations")
    if isinstance(payload, dict) and payload:
        return payload
    summary_path = data.get("summary_path")
    if not summary_path or not Path(str(summary_path)).exists():
        return {}
    try:
        summary = json.loads(Path(str(summary_path)).read_text(encoding="utf-8"))
    except Exception:
        return {}
    payload = summary.get("yield_improvement_recommendations") if isinstance(summary, dict) else None
    return payload if isinstance(payload, dict) else {}


def _render_protein_cost_analysis(protein_cost: dict[str, object]) -> None:
    with st.expander("目标蛋白成本分析", expanded=True):
        st.caption("解释型相对评分，不代表真实发酵产量、培养成本或湿实验结果。")
        c1, c2, c3 = st.columns([1, 2, 2])
        with c1:
            st.metric("总相对成本分", protein_cost.get("total_relative_score", "—"))
        with c2:
            st.write("**主要成本类别**")
            st.write(", ".join(str(item) for item in protein_cost.get("dominant_cost_categories") or []) or "—")
        with c3:
            st.write("**状态**")
            st.write(protein_cost.get("result_status", "draft_explanatory"))

        items = protein_cost.get("cost_items") or []
        if items:
            frame = pd.DataFrame(items)
            display_columns = {
                "category": "类别",
                "label": "成本项",
                "relative_score": "相对分",
                "basis": "依据",
                "interpretation": "解释",
                "raw_value": "原始值",
            }
            columns = [key for key in display_columns if key in frame.columns]
            st.dataframe(
                frame[columns].rename(columns=display_columns),
                use_container_width=True,
                hide_index=True,
            )
        lp_attribution = protein_cost.get("lp_attribution")
        if isinstance(lp_attribution, dict) and lp_attribution:
            _render_lp_attribution(lp_attribution)
        cost_slope = protein_cost.get("cost_slope_compatibility")
        if isinstance(cost_slope, dict) and cost_slope:
            _render_cost_slope_compatibility(cost_slope)
        for warning in protein_cost.get("warnings") or []:
            st.warning(str(warning))


def _render_lp_attribution(lp_attribution: dict[str, object]) -> None:
    st.markdown("**LP 级归因证据**")
    st.caption("Python draft LP sensitivity，基于 SciPy HiGHS marginals；不是 MATLAB/SoPlex fully aligned shadow price。")
    objective = lp_attribution.get("objective_evidence") if isinstance(lp_attribution.get("objective_evidence"), dict) else {}
    c1, c2, c3 = st.columns([1, 1, 1])
    with c1:
        st.metric("LP 归因状态", lp_attribution.get("result_status", "—"))
    with c2:
        st.metric("目标反应", objective.get("objective_reaction", "—") if isinstance(objective, dict) else "—")
    with c3:
        st.metric("分泌通量", objective.get("secretion_flux", "—") if isinstance(objective, dict) else "—")

    sections = (
        ("主导约束块", "dominant_constraint_blocks"),
        ("Top constraint marginals", "top_constraint_marginals"),
        ("Top bound marginals", "top_bound_marginals"),
        ("目标相关 flux", "target_related_fluxes"),
    )
    for title, key in sections:
        rows = lp_attribution.get(key) or []
        if rows:
            st.write(f"**{title}**")
            st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)
    counts = lp_attribution.get("active_bound_counts")
    if isinstance(counts, dict) and counts:
        st.write("**Active bound marginal counts**")
        st.dataframe(
            pd.DataFrame([{"项目": key, "数量": value} for key, value in counts.items()]),
            use_container_width=True,
            hide_index=True,
        )
    for warning in lp_attribution.get("warnings") or []:
        st.warning(str(warning))


def _render_cost_slope_compatibility(cost_slope: dict[str, object]) -> None:
    st.markdown("**MATLAB-compatible 蛋白成本 slope（可选）**")
    st.caption(
        "当前默认路线是固定生长率、corrected medium、最大化目标蛋白分泌通量；"
        "历史 MATLAB 成本路线是固定目标分泌比例和生长率、优化 Ex_glc_D，用于 Protein_cost_TP 定义对比。"
    )
    c1, c2, c3 = st.columns([1, 1, 1])
    with c1:
        st.metric("开启状态", str(cost_slope.get("enabled", False)))
    with c2:
        st.metric("结果状态", cost_slope.get("result_status", "—"))
    with c3:
        st.metric("成功", str(cost_slope.get("success", "—")))
    _render_cost_slope_ratio_policy(cost_slope)
    st.caption(
        f"培养基兼容模式: {cost_slope.get('medium_compatibility_mode', 'corrected')}；"
        "该设置只影响可选 cost slope 对比，不改变默认分泌仿真。"
    )
    overrides = cost_slope.get("medium_bound_overrides") or []
    if overrides:
        st.write("**MATLAB legacy medium bound overrides**")
        st.dataframe(pd.DataFrame(overrides), use_container_width=True, hide_index=True)

    sections = (
        ("Glucose cost slopes", "glucose_cost_slopes"),
        ("Ribosome cost slopes", "ribosome_cost_slopes"),
        ("Cost slope rows", "rows"),
    )
    for title, key in sections:
        rows = cost_slope.get(key) or []
        if rows:
            st.write(f"**{title}**")
            st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)
    comparison_scope = cost_slope.get("comparison_scope")
    if isinstance(comparison_scope, dict) and comparison_scope:
        with st.expander("对比定义", expanded=False):
            st.json(comparison_scope)
    for warning in cost_slope.get("warnings") or []:
        st.warning(str(warning))


def _render_cost_slope_ratio_policy(cost_slope: dict[str, object]) -> None:
    policy = str(cost_slope.get("secretion_ratio_policy") or "explicit_absolute_ratios")
    capacity = cost_slope.get("capacity_reference")
    fractions = tuple(cost_slope.get("capacity_fractions") or ())
    if policy == "capacity_fraction_ratios":
        fraction_text = ", ".join(f"{float(value):.0%}" for value in fractions)
        st.info(
            "目标分泌比例来源: 未提供实验或用户指定比例，因此按当前 corrected 分泌 capacity "
            f"{capacity} 的 {fraction_text} 自动生成成本斜率网格。"
        )
    elif policy == "explicit_absolute_ratios":
        st.info("目标分泌比例来源: 使用请求中显式提供的绝对分泌比例，作为历史 MATLAB-style 固定需求。")
    else:
        st.warning("目标分泌比例来源: 当前 capacity 不可用，已退回历史绝对比例；该结果只适合作诊断参考。")


def _render_target_growth_analysis(target_growth: dict[str, object]) -> None:
    with st.expander("目标蛋白生长分析", expanded=True):
        st.caption("解释型 small-grid tradeoff，不代表真实发酵生长预测。")
        best_flux = target_growth.get("best_secretion_point") or {}
        best_per_biomass = target_growth.get("best_secretion_per_biomass_point") or {}
        c1, c2, c3 = st.columns([1, 1, 1])
        with c1:
            st.metric("趋势标签", target_growth.get("growth_sensitivity_label", "—"))
            reason = target_growth.get("growth_sensitivity_reason")
            if reason:
                st.caption(f"原因: {reason}")
        with c2:
            st.metric("最高分泌通量 mu", best_flux.get("mu", "—") if isinstance(best_flux, dict) else "—")
        with c3:
            st.metric("最高单位生物量 mu", best_per_biomass.get("mu", "—") if isinstance(best_per_biomass, dict) else "—")

        points = target_growth.get("tradeoff_points") or []
        if points:
            frame = pd.DataFrame(points)
            display_columns = {
                "mu": "生长速率 mu",
                "success": "成功",
                "secretion_flux": "分泌通量",
                "secretion_per_biomass": "单位生物量分泌",
                "status": "求解状态",
                "interpretation": "解释",
            }
            columns = [key for key in display_columns if key in frame.columns]
            st.dataframe(
                frame[columns].rename(columns=display_columns),
                use_container_width=True,
                hide_index=True,
            )
        for warning in target_growth.get("warnings") or []:
            st.warning(str(warning))


def _render_yield_improvement_recommendations(payload: dict[str, object]) -> None:
    with st.expander("目标蛋白产量提升推荐", expanded=True):
        st.caption("Python corrected draft 模型内推荐，不代表真实发酵产量或实验成功率。")
        counts = payload.get("summary_counts") if isinstance(payload.get("summary_counts"), dict) else {}
        c1, c2, c3 = st.columns([1, 1, 1])
        with c1:
            st.metric("推荐候选", counts.get("recommended", 0) if isinstance(counts, dict) else 0)
        with c2:
            st.metric("需人工/暂不推荐", counts.get("not_recommended", 0) if isinstance(counts, dict) else 0)
        with c3:
            st.metric("未解析", counts.get("unresolved", 0) if isinstance(counts, dict) else 0)

        rows = payload.get("recommended_candidates") or []
        if rows:
            frame = pd.DataFrame(rows)
            display_columns = {
                "recommendation_tier": "证据分级",
                "recommendation_label": "推荐等级",
                "display_name": "候选",
                "intervention_type": "扰动",
                "execution_mode": "执行模式",
                "delta_objective": "Δobjective",
                "secretory_process": "分泌环节",
                "database_annotation_sources": "数据库注释来源",
                "model_gpr_executable": "模型 GPR 可执行",
                "oe_reaction_proxy": "OE 反应代理",
                "evidence_tier": "证据等级",
                "recommendation_score": "推荐分",
                "rationale": "推荐理由",
            }
            columns = [key for key in display_columns if key in frame.columns]
            st.dataframe(frame[columns].rename(columns=display_columns), use_container_width=True, hide_index=True)
        else:
            st.info("当前候选没有进入模型内提升推荐。")
        st.caption("gene-level KO 与 reaction-level OE proxy 是不同证据层级；OE proxy 不能直接等同于湿实验基因过表达。")
        for warning in payload.get("warnings") or []:
            st.warning(str(warning))


__all__ = ["render_pichia_results"]
