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
    # Only rendered when 启用蛋白成本斜率对比 was on (pipeline.py sets protein_cost_analysis
    # to None otherwise) - everything shown here is LP-solve-derived, not a heuristic score.
    with st.expander("目标蛋白成本分析", expanded=True):
        st.caption(f"状态: {protein_cost.get('result_status', 'draft_cost_slope_analysis')}")
        lp_attribution = protein_cost.get("lp_attribution")
        if isinstance(lp_attribution, dict) and lp_attribution:
            _render_lp_attribution(lp_attribution)
        solver_robustness = protein_cost.get("solver_robustness")
        if isinstance(solver_robustness, dict) and solver_robustness:
            _render_solver_robustness(solver_robustness)
        oe_dose_response = protein_cost.get("oe_dose_response")
        if isinstance(oe_dose_response, dict) and oe_dose_response:
            _render_oe_dose_response(oe_dose_response)
        cost_slope = protein_cost.get("cost_slope_compatibility")
        if isinstance(cost_slope, dict) and cost_slope:
            _render_cost_slope_compatibility(cost_slope)
        for warning in protein_cost.get("warnings") or []:
            st.warning(str(warning))


def _render_solver_robustness(solver_robustness: dict[str, object]) -> None:
    st.markdown("**求解器稳健性（瓶颈归因是否跨求解器稳定）**")
    st.caption(
        "对偶解在退化最优解处不唯一，换求解算法可能把瓶颈归到不同资源。"
        "ranking-insensitive-to-solver 表示归因稳定；ranking-sensitive-to-solver 表示归因是数值假象，不是生物学结论。"
    )
    classification = str(solver_robustness.get("classification", "—"))
    c1, c2, c3 = st.columns([2, 1, 1])
    with c1:
        st.metric("稳健性分类", classification)
    with c2:
        st.metric("瓶颈跨求解器一致", str(solver_robustness.get("top_bottleneck_stable", "—")))
    with c3:
        st.metric("主导块跨求解器一致", str(solver_robustness.get("top_block_stable", "—")))
    if classification == "ranking-sensitive-to-solver":
        st.warning(str(solver_robustness.get("detail", "")))
    else:
        st.caption(str(solver_robustness.get("detail", "")))
    per_method = solver_robustness.get("per_method") or []
    if per_method:
        st.dataframe(pd.DataFrame(per_method), use_container_width=True, hide_index=True)
    for warning in solver_robustness.get("warnings") or []:
        st.warning(str(warning))


_OE_DOSE_RESPONSE_SHAPE_LABELS = {
    "saturating": "饱和型（适度过表达就够，再加收益递减）",
    "linear": "线性型（还在涨，值得进一步加大表达）",
    "threshold": "阈值型（要超过某个最小倍数才起效）",
    "flat_no_response": "无响应（任何倍数都几乎没提升，别过表达）",
    "non_monotonic_numerical_artifact": "非单调（数值假象，不可作结论）",
    "insufficient_points": "数据点不足",
}


def _render_oe_dose_response(oe_dose_response: dict[str, object]) -> None:
    st.markdown("**OE 剂量响应形状（过表达越多，分泌是持续上升还是很快到顶）**")
    st.caption(
        "对候选反应扫描一组过表达倍数，把分泌响应的形状分类——替代只在固定 2× 一个点看提升。"
        "这是相对形状信号，只说明趋势，不产出绝对产量或最优倍数。"
    )
    baseline = oe_dose_response.get("baseline_objective")
    baseline_txt = f"{float(baseline):.4g}" if isinstance(baseline, (int, float)) else "—"
    factors = oe_dose_response.get("tested_factors") or []
    factors_txt = "、".join(f"{float(f):g}×" for f in factors) or "—"
    c1, c2 = st.columns([1, 2])
    with c1:
        st.metric("无 OE 基线分泌目标", baseline_txt)
    with c2:
        st.metric("扫描的过表达倍数", factors_txt)

    def _pct(value: object) -> str:
        return f"{value * 100:.2f}%" if isinstance(value, (int, float)) else "—"

    shapes = [s for s in (oe_dose_response.get("reaction_shapes") or []) if isinstance(s, dict)]
    rows = [
        {
            "反应": shape.get("reaction_id", "—"),
            "形状": _OE_DOSE_RESPONSE_SHAPE_LABELS.get(str(shape.get("shape")), str(shape.get("shape"))),
            "最大相对增益": _pct(shape.get("max_relative_gain")),
            "半增益倍数(拐点)": (
                f"{float(shape['half_gain_factor']):g}×" if shape.get("half_gain_factor") is not None else "—"
            ),
            "最强剂量处相对增益": _pct(shape.get("relative_gain_at_max_factor")),
            "单调不减": shape.get("monotonic_non_decreasing", "—"),
        }
        for shape in shapes
    ]
    if rows:
        st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)
    else:
        st.caption("（没有可用的剂量响应结果）")
    # A decrease that OE cannot truly cause is a degenerate-optimum artifact: flag it explicitly.
    for shape in shapes:
        if shape.get("shape") == "non_monotonic_numerical_artifact":
            st.warning(f"{shape.get('reaction_id')}：{shape.get('detail', '')}")
    for warning in oe_dose_response.get("warnings") or []:
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

    oe_actionable = lp_attribution.get("oe_actionable_bottlenecks") or []
    floor_only = lp_attribution.get("floor_constraints_not_oe_addressable") or []
    st.write("**OE 可缓解瓶颈（binding 上限，按复合体）**")
    st.caption(
        "只列当前解处 binding 的上限产能约束——这是 OE（放宽产能上限）真能缓解的瓶颈线索，按复合体给出。"
        "注意这只是线索不是保证：耦合结构下放宽一个上限会让瓶颈转移，需与真实 reaction_oe_tradeoff 交叉验证。"
    )
    if oe_actionable:
        st.dataframe(pd.DataFrame(oe_actionable), use_container_width=True, hide_index=True)
    else:
        st.caption("（当前解没有 binding 的上限产能约束）")
    if floor_only:
        st.write("**下界/floor 约束（OE 动不了，仅供参考）**")
        st.caption(
            "这些是下界（最低要求类）约束，OE 放宽的是上限、对它们无效，绝不能当作 OE 靶点。"
            "此前 PDI1 单体、核糖体装配就是落在这里的假阳性。"
        )
        st.dataframe(pd.DataFrame(floor_only), use_container_width=True, hide_index=True)

    sections = (
        ("主导约束块", "dominant_constraint_blocks"),
        ("Top constraint marginals", "top_constraint_marginals"),
        ("Top bound marginals（未按 bound_type 分离的原始表）", "top_bound_marginals"),
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
                "standard_symbol": "标准符号",
                "gene_display_name": "标准显示名",
                "protein_name": "蛋白名称",
                "annotation_confidence": "命名置信度",
                "standard_name_status": "命名状态",
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
