from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import plotly.express as px
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
    value_of_information = _value_of_information_payload(data)
    if value_of_information:
        _render_value_of_information(value_of_information)

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


def _value_of_information_payload(data: dict[str, object]) -> dict[str, object]:
    payload = data.get("value_of_information")
    if isinstance(payload, dict) and payload:
        return payload
    summary_path = data.get("summary_path")
    if not summary_path or not Path(str(summary_path)).exists():
        return {}
    try:
        summary = json.loads(Path(str(summary_path)).read_text(encoding="utf-8"))
    except Exception:
        return {}
    payload = summary.get("value_of_information") if isinstance(summary, dict) else None
    return payload if isinstance(payload, dict) else {}


_VOI_AMBIGUITY_LABELS = {
    "near_tie": "近似并列（模型分不清谁更好）",
    "capacity_flip": "跨容量假设翻转",
    "solver_flip": "跨求解器翻转",
}


def _render_value_of_information(payload: dict[str, object]) -> None:
    items = [item for item in (payload.get("information_items") or []) if isinstance(item, dict)]
    ranked = [row for row in (payload.get("ranked_candidates") or []) if isinstance(row, dict)]
    if not items and not ranked:
        return
    with st.expander("排序可信度 & 该测什么（价值-of-information）", expanded=True):
        st.caption(
            "模型给的是相对排序、不是绝对产量。这里标出顶部名次里模型分不清的候选（排序不可信），"
            "并给出最能消解歧义的最小湿实验——只排测量优先级，不预测结果、不自动认定谁更好。"
        )
        if payload.get("has_actionable_ambiguity"):
            st.warning(
                f"顶部排序有 {len(items)} 处歧义（近似并列或跨假设翻转），当前名次不完全可信——见下方建议测量。"
            )
        else:
            st.success("顶部候选分数分离明显，当前相对排序较可信（未检测到近似并列或翻转）。")

        frame = pd.DataFrame(
            [
                {"候选": str(row.get("candidate_id", "?")), "推荐分数": float(row.get("score") or 0.0), "名次": int(row.get("rank") or 0)}
                for row in ranked
                if row.get("score") is not None
            ]
        )
        if not frame.empty:
            head = frame.sort_values("名次").head(10)
            figure = px.bar(
                head,
                x="推荐分数",
                y="候选",
                orientation="h",
                title="候选排序：分数越接近越难区分（越该做湿实验测）",
            )
            figure.update_layout(
                xaxis_title="推荐分数（相对，非绝对产量）",
                yaxis_title="",
                yaxis={"categoryorder": "total ascending"},
            )
            st.plotly_chart(figure, use_container_width=True)

        if items:
            rows = [
                {
                    "优先级": item.get("priority_rank"),
                    "歧义类型": _VOI_AMBIGUITY_LABELS.get(str(item.get("ambiguity_kind")), str(item.get("ambiguity_kind"))),
                    "涉及候选": "、".join(str(candidate) for candidate in (item.get("candidates") or [])),
                    "建议测量": item.get("suggested_measurement"),
                }
                for item in items
            ]
            st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)
        for warning in payload.get("warnings") or []:
            st.caption(str(warning))


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

# Chinese labels for the LP secretory-process / resource-layer codes. Kept local because the
# Streamlit facade may not import the engine (see test_streamlit_ui_does_not_import_engine_directly).
_SECRETORY_PROCESS_LABELS = {
    "ribosome": "翻译（核糖体）",
    "proteasome_degradation": "蛋白降解（蛋白酶体）",
    "disulfide_folding": "二硫键折叠 / DSB",
    "n_glycan_processing": "N-糖基化",
    "o_glycan_processing": "O-糖基化",
    "chaperone_folding": "分子伴侣折叠",
    "erad_misfolding": "错误折叠 / ERAD",
    "er_translocation": "ER 转运",
    "er_to_golgi_transport": "ER→Golgi 转运",
    "golgi_surface_transport": "Golgi→胞外运输",
    "secretory_capacity": "分泌容量",
    "metabolic_or_other": "代谢 / 其它",
    "unknown": "未解析",
}


def _oe_dose_response_curve_frame(oe_dose_response: dict[str, object]) -> pd.DataFrame:
    """Flatten the per-reaction dose-response points into a long frame for a line chart.

    One row per (reaction, factor): factor on x, relative gain (%) on y. The factor-1.0 baseline
    is included so each curve starts at 0%. Points with a missing factor/gain are dropped.
    """
    rows: list[dict[str, object]] = []
    for shape in oe_dose_response.get("reaction_shapes") or []:
        if not isinstance(shape, dict):
            continue
        reaction = str(shape.get("reaction_id", "?"))
        shape_label = _OE_DOSE_RESPONSE_SHAPE_LABELS.get(str(shape.get("shape")), str(shape.get("shape")))
        legend = f"{reaction}｜{shape_label}"
        for point in shape.get("point_deltas") or []:
            if not isinstance(point, dict):
                continue
            factor = point.get("factor")
            gain = point.get("relative_gain")
            if factor is None or gain is None:
                continue
            rows.append(
                {
                    "过表达倍数": float(factor),
                    "分泌相对提升(%)": float(gain) * 100.0,
                    "反应｜形状": legend,
                }
            )
    return pd.DataFrame(rows)


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

    curve = _oe_dose_response_curve_frame(oe_dose_response)
    if not curve.empty:
        figure = px.line(
            curve.sort_values("过表达倍数"),
            x="过表达倍数",
            y="分泌相对提升(%)",
            color="反应｜形状",
            markers=True,
            title="OE 剂量响应曲线：过表达倍数越高，分泌相对提升怎么走",
        )
        figure.update_layout(
            xaxis_title="过表达倍数（×，1 = 不过表达）",
            yaxis_title="分泌相对提升（%）",
            legend_title_text="反应｜形状",
        )
        st.plotly_chart(figure, use_container_width=True)

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


def _lp_oe_bottleneck_frame(lp_attribution: dict[str, object]) -> pd.DataFrame:
    """Frame of OE-actionable bottlenecks for a horizontal bar chart (which resource layer binds).

    Only the OE-actionable (binding upper-bound) entries are charted — the lower-bound floors are
    not OE targets and are deliberately excluded so the chart cannot suggest acting on them.
    """
    rows: list[dict[str, object]] = []
    for entry in lp_attribution.get("oe_actionable_bottlenecks") or []:
        if not isinstance(entry, dict):
            continue
        marginal = entry.get("abs_marginal")
        if marginal is None:
            marginal = abs(float(entry.get("marginal") or 0.0))
        rows.append(
            {
                "反应": str(entry.get("reaction_id", "?")),
                "影子价格(绝对值)": abs(float(marginal or 0.0)),
                "分泌资源层": _SECRETORY_PROCESS_LABELS.get(
                    str(entry.get("secretory_process", "unknown")), str(entry.get("secretory_process", "unknown"))
                ),
            }
        )
    return pd.DataFrame(rows)


def _lp_floor_bottleneck_frame(lp_attribution: dict[str, object]) -> pd.DataFrame:
    """Frame of the largest binding LOWER-bound floors — the 'why is it limited' signal.

    These are minimum-requirement constraints (folding, translation, ERAD) that OE cannot relax,
    but they carry the largest shadow prices and are the honest answer to which resource layer
    dominates a target's limitation, so they are charted separately from the OE-actionable ceilings.
    """
    rows: list[dict[str, object]] = []
    for entry in lp_attribution.get("floor_constraints_not_oe_addressable") or []:
        if not isinstance(entry, dict):
            continue
        marginal = entry.get("abs_marginal")
        if marginal is None:
            marginal = abs(float(entry.get("marginal") or 0.0))
        rows.append(
            {
                "反应": str(entry.get("reaction_id", "?")),
                "影子价格(绝对值)": abs(float(marginal or 0.0)),
                "分泌资源层": _SECRETORY_PROCESS_LABELS.get(
                    str(entry.get("secretory_process", "unknown")), str(entry.get("secretory_process", "unknown"))
                ),
            }
        )
    return pd.DataFrame(rows)


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
        bottleneck_frame = _lp_oe_bottleneck_frame(lp_attribution)
        if not bottleneck_frame.empty:
            figure = px.bar(
                bottleneck_frame.sort_values("影子价格(绝对值)"),
                x="影子价格(绝对值)",
                y="反应",
                color="分泌资源层",
                orientation="h",
                title="OE 可缓解瓶颈：哪一层最限制分泌（影子价格绝对值越大越紧）",
            )
            figure.update_layout(
                xaxis_title="影子价格绝对值（越大越限制分泌）",
                yaxis_title="",
                legend_title_text="分泌资源层",
            )
            st.plotly_chart(figure, use_container_width=True)
        st.dataframe(pd.DataFrame(oe_actionable), use_container_width=True, hide_index=True)
    else:
        st.caption("（当前解没有 binding 的上限产能约束）")
    if floor_only:
        st.write("**为什么受限：最强约束层（下界/最低要求，OE 动不了）**")
        st.caption(
            "这些下界（折叠 / 翻译 / ERAD 等最低要求）承载最大的影子价格，是本靶点「为什么受限、卡在哪一层」的答案；"
            "但 OE 放宽的是上限、对它们无效，绝不能当 OE 靶点（PDI1 单体、核糖体装配就是这里的经典假阳性）。"
        )
        floor_frame = _lp_floor_bottleneck_frame(lp_attribution)
        if not floor_frame.empty:
            figure = px.bar(
                floor_frame.sort_values("影子价格(绝对值)"),
                x="影子价格(绝对值)",
                y="反应",
                color="分泌资源层",
                orientation="h",
                title="为什么受限：最强约束层（影子价格越大越限制；OE 动不了）",
            )
            figure.update_layout(
                xaxis_title="影子价格绝对值（越大越限制分泌）",
                yaxis_title="",
                legend_title_text="分泌资源层",
            )
            st.plotly_chart(figure, use_container_width=True)
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
