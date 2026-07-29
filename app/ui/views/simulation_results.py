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
from app.core.i18n import sim_result_column_label, sim_result_value_label, sim_result_warning_label
from app.services.per_strain_oe_candidate_run import run_next_oe_candidate_analysis
from app.services.per_strain_shortlist_run import build_modified_strain_shortlist, recompute_stale_candidates
from app.ui.common import PATHS
from app.ui.views.candidate_path_graph import render_secretion_path_graph
from app.ui.views.simulation_display import (
    CANDIDATE_DISPLAY_COLUMNS,
    candidate_effect_counts,
    candidate_row_label,
    normalise_candidate_frame_for_display,
)


def _localized_frame(records: object, value_columns: tuple[str, ...] = ()) -> pd.DataFrame:
    """Raw payload records -> display DataFrame with Chinese column names + mapped enum values.

    Single funnel so the engine payload's English field names / enum codes never leak onto the
    results page. Values in ``value_columns`` are translated via the central i18n value dict; all
    column names via the central column dict. Unknown keys/values fall back to their raw text.
    """
    frame = pd.DataFrame(list(records or []))
    if frame.empty:
        return frame
    for column in value_columns:
        if column in frame.columns:
            frame[column] = frame[column].map(sim_result_value_label)
    return frame.rename(columns={column: sim_result_column_label(column) for column in frame.columns})


def _fmt_num(value: object) -> str:
    """数值 → 3 位有效数字；真的不是数才给 —（结果页最常看的顶部指标统一格式，避免长浮点串）。

    必须接受**数字字符串**：后台任务把 objective_value 原样透传，而它在响应里是字符串
    （如 '0.00429116390133443'，为保留完整精度）。此前只认 int/float，导致"分泌通量"顶部指标
    长期显示 —，而同一个数在下方分析里却正常出现（2026-07-29 用户报出）。
    """
    if isinstance(value, bool):
        return "—"
    if isinstance(value, (int, float)):
        return f"{value:.3g}"
    if isinstance(value, str):
        try:
            return f"{float(value.strip()):.3g}"
        except ValueError:
            return "—"
    return "—"


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
        st.metric("分泌通量", _fmt_num(data.get("objective_value")), help="相对比较值，非实际产量")
    with c2:
        st.metric("目标蛋白", data.get("target_id", "—"))
    with c3:
        # 自定义构建根本没有 MATLAB 参照物，显示"待处理"会让人以为还有一步没做。
        alignment_status = str(data.get("matlab_alignment_status") or "")
        is_custom_target = not bool((data.get("target_metadata") or {}).get("is_builtin", True))
        if alignment_status == "pending" and is_custom_target:
            st.metric("与历史实现对照", "不适用")
        else:
            st.metric("与历史实现对照", sim_result_value_label(alignment_status))
    st.caption(
        "不同构建之间可横向对比。不代表实际发酵产量。"
        "「与历史实现对照」只说明这次结果能否与旧 MATLAB 版逐点核对，与结果好坏无关；自定义构建没有对照物。"
    )

    # 默认只给主结论。参数/注意事项/生长权衡/原始数据这类诊断信息此前各占一个折叠区常驻页面——
    # 折叠区本身也是认知成本（得先猜里面有什么才决定点不点）。合并成一个总开关，默认关。
    show_diagnostics = st.checkbox(
        "显示诊断细节（参数 · 注意事项 · 生长权衡 · 原始数据）",
        value=False,
        key="pichia_results_show_diagnostics",
    )

    if show_diagnostics:
        with st.expander("参数", expanded=False):
            medium_condition = data.get("medium_condition") if isinstance(data.get("medium_condition"), dict) else {}
            st.dataframe(
                pd.DataFrame(
                    [
                        {"参数": key, "值": sim_result_value_label(value)}
                        for key, value in {
                            "目标": data.get("target_id"),
                            "状态": data.get("result_status"),
                            "MATLAB": data.get("matlab_alignment_status"),
                            "目标值": data.get("objective_value"),
                            "培养基条件": medium_condition.get("condition_id") if isinstance(medium_condition, dict) else None,
                            "碳源": medium_condition.get("carbon_source_id") if isinstance(medium_condition, dict) else None,
                            "培养基状态": medium_condition.get("status") if isinstance(medium_condition, dict) else None,
                            "科学解释状态": medium_condition.get("scientific_status") if isinstance(medium_condition, dict) else None,
                            "碳源标定档": (medium_condition.get("carbon_source_formulation") or {}).get("formulation_status") if isinstance(medium_condition, dict) else None,
                        }.items()
                        if value
                    ]
                ),
                width='stretch',
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
            width='stretch',
            hide_index=True,
        )

    warns = data.get("warnings") or []
    if warns and show_diagnostics:
        with st.expander("注意事项", expanded=False):
            for warning in warns:
                st.warning(sim_result_warning_label(warning))
    errors = data.get("errors") or []
    if errors:
        with st.expander("错误", expanded=True):
            for error in errors:
                st.error(error)

    protein_cost = _protein_cost_payload(data)
    if protein_cost:
        _render_protein_cost_analysis(protein_cost)
    _render_next_oe_candidates(data)
    _render_modified_strain_shortlist(data)
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
    if show_diagnostics and tradeoff_path and Path(tradeoff_path).exists():
        with st.expander("生长权衡", expanded=False):
            st.dataframe(
                pd.read_csv(tradeoff_path).rename(columns=sim_result_column_label),
                width='stretch',
                hide_index=True,
            )
    if show_diagnostics:
        with st.expander("完整结果数据（高级）", expanded=False):
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
    st.dataframe(display, width='stretch', hide_index=True)

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


def _render_value_of_information(payload: dict[str, object]) -> None:
    items = [item for item in (payload.get("information_items") or []) if isinstance(item, dict)]
    ranked = [row for row in (payload.get("ranked_candidates") or []) if isinstance(row, dict)]
    if not items and not ranked:
        return
    with st.expander("排序可信度 & 该测什么", expanded=True):
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
                {"候选": _short_reaction_label(str(row.get("candidate_id", "?"))), "推荐分数": float(row.get("score") or 0.0), "名次": int(row.get("rank") or 0)}
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
                text="推荐分数",
                color_discrete_sequence=["#0F766E"],
                title="候选排序：分数越接近越难区分（越该做湿实验测）",
            )
            figure.update_traces(texttemplate="%{text:.3g}", textposition="outside", cliponaxis=False)
            figure.update_layout(
                xaxis_title="推荐分数（相对，非绝对产量）",
                yaxis_title="",
                yaxis={"categoryorder": "total ascending"},
            )
            st.plotly_chart(figure, width='stretch')

        if items:
            rows = []
            for item in items:
                candidate_text = "、".join(
                    _short_reaction_label(str(candidate)) for candidate in (item.get("candidates") or [])
                )
                # Rebuild the measurement suggestion in Chinese from structured fields — the engine
                # emits it in English (prioritize_value_of_information); the UI owns the localized wording.
                measurement = (
                    f"对候选 {candidate_text} 做靶点特异的分泌定量湿实验，消解它们的相对次序"
                    + ("（影响 top 名次，优先做）" if item.get("resolves_top_of_ranking") else "")
                )
                rows.append(
                    {
                        "优先级": item.get("priority_rank"),
                        "歧义类型": sim_result_value_label(item.get("ambiguity_kind")),
                        "涉及候选": candidate_text,
                        "建议测量": measurement,
                    }
                )
            st.dataframe(pd.DataFrame(rows), width='stretch', hide_index=True)
        for warning in payload.get("warnings") or []:
            st.caption(sim_result_warning_label(warning))


def _render_next_oe_candidates(data: dict[str, object]) -> None:
    """针对当前改造菌株的"下一步 OE 候选"：opt-in 重解 + 有界剂量响应（ADR-004 #1 迭代候选）。

    读上次运行时暂存的 KO/OE 改造参数（`pichia_last_run_modifications`），用户点按钮后才触发
    额外求解（改造后重解 → 瓶颈复合体 → 有界 OE 剂量响应 → 排序），结果缓存在 session。
    """
    mods = st.session_state.get("pichia_last_run_modifications")
    with st.expander("下一步 OE 候选（针对当前改造菌株）", expanded=False):
        st.caption(
            "把本次构建里设定的 KO/OE 应用到模型后重新求解，找改造后这一株当前卡在上限的产能瓶颈复合体，"
            "再对它们扫有界过表达剂量响应，按真实相对效应排出下一步该 OE 谁。瓶颈会随改造转移——每改一轮都应重跑。"
            "相对信号、复合体级、非绝对产量。"
        )
        if not isinstance(mods, dict) or not mods.get("target_id"):
            st.info("先在「仿真构建」页运行一次仿真，这里就能基于同一株算下一步 OE 候选。")
            return
        if not mods.get("target_is_builtin"):
            # 说清楚"为什么不能"和"怎么才能"，别只说不支持——三段式成为默认后这条会经常出现。
            st.info(
                "这一步需要按目标 ID 重新装配模型，而三段式 / 自定义 JSON 构建的目标没有可复原的内置定义，"
                "所以暂时算不了。想用它：在「① 目标蛋白」里改用**快速选择（内置模板）**跑一次即可。"
            )
            return

        oe_rx = [str(r) for r in (mods.get("oe_reaction_ids") or [])]
        ko_rx = [str(r) for r in (mods.get("ko_reaction_ids") or [])]
        gene_mods = [str(g) for g in (mods.get("oe_gene_ids") or [])] + [str(g) for g in (mods.get("ko_gene_ids") or [])]
        st.markdown(
            f"**改造后菌株**（将被重解）：OE 反应 `{len(oe_rx)}` 个、KO 反应 `{len(ko_rx)}` 个；过表达按 2× 产能建模。"
        )
        if not oe_rx and not ko_rx:
            st.caption("本次没有反应级 KO/OE —— 将按无额外改造分析，得到的是野生型当前瓶颈。")
        if gene_mods:
            st.caption(f"注意：{len(gene_mods)} 个基因级改造暂不纳入重解（只应用反应/复合体级）。")

        if st.button("计算下一步 OE 候选（会额外求解，约数十秒）", key="pichia_next_oe_candidates_run_button"):
            with st.spinner("重解改造后菌株并扫描瓶颈复合体剂量响应…"):
                st.session_state["pichia_next_oe_candidates_result"] = run_next_oe_candidate_analysis(
                    target_id=str(mods.get("target_id")),
                    ko_reaction_ids=tuple(ko_rx),
                    oe_reaction_ids=tuple(oe_rx),
                    mu=float(mods.get("mu") or 0.10),
                    media_type=int(mods.get("media_type") or 4),
                    carbon_source_id=str(mods.get("carbon_source_id") or "glucose"),
                    enable_ribosome_translation_constraint=bool(mods.get("enable_ribosome")),
                    enable_misfolding_constraint=bool(mods.get("enable_misfolding")),
                )

        readout = st.session_state.get("pichia_next_oe_candidates_result")
        if not isinstance(readout, dict):
            return
        if readout.get("error"):
            st.error(f"计算失败：{readout.get('error')}")
            return
        if not readout.get("modified_solve_success"):
            st.warning("改造后菌株在当前约束下无可行解（或无 LP 灵敏度），无法给出瓶颈候选。")
        _render_next_oe_candidates_result(readout)


def _render_next_oe_candidates_result(readout: dict[str, object]) -> None:
    c1, c2, c3 = st.columns([1, 1, 1])
    with c1:
        st.metric("改造后分泌目标", _fmt_num(readout.get("modified_objective_value")))
    with c2:
        st.metric("排序依据", "真实剂量效应" if readout.get("dose_response_available") else "限制强度")
    with c3:
        st.metric("碳源", sim_result_value_label(readout.get("carbon_source_id")))

    candidates = [c for c in (readout.get("candidates") or []) if isinstance(c, dict)]
    if candidates:
        dose_available = bool(readout.get("dose_response_available"))
        rows: list[dict[str, object]] = []
        for rank, candidate in enumerate(candidates, 1):
            row: dict[str, object] = {
                "排名": rank,
                "OE 目标(复合体)": _short_reaction_label(str(candidate.get("reaction", "—"))),
                "分泌资源层": _resource_layer_label(candidate.get("layer")),
                "限制强度": _fmt_num(candidate.get("shadow_price")),
            }
            if dose_available:
                effect = candidate.get("effect")
                row["剂量响应形状"] = sim_result_value_label(candidate.get("shape")) if candidate.get("shape") else "—"
                row["最大相对增益"] = f"{effect * 100:.2f}%" if isinstance(effect, (int, float)) else "—"
                row["半增益倍数(拐点)"] = (
                    f"{float(candidate['half_gain_factor']):g}×" if candidate.get("half_gain_factor") is not None else "—"
                )
            rows.append(row)
        st.dataframe(pd.DataFrame(rows), width='stretch', hide_index=True)
        if not dose_available:
            st.caption("未取得剂量响应（瓶颈复合体无有界扫描结果）——暂按限制强度排序，只表示谁在卡着（生效约束），不表示松开后能涨多少。")
    else:
        st.caption("（改造后未发现过表达可缓解的、卡在上限的瓶颈复合体。）")

    for caveat in readout.get("caveats") or []:
        st.caption(f"• {caveat}")
    floors = readout.get("floor_constraints_not_oe_addressable") or []
    if floors:
        st.caption(f"另有 {len(floors)} 个下界约束（最低需求，如折叠/翻译）——过表达松不动，不列为候选。")
    for warning in readout.get("modification_warnings") or []:
        st.warning(str(warning))


_REUSE_MODULE_LABELS = {
    "folding": "折叠", "glycosylation": "糖基化", "transport": "转运", "translation": "翻译",
    "degradation": "降解", "secretory_capacity": "分泌容量", "metabolic": "代谢/其它", "unknown": "未解析",
}


def _reuse_module_label(module: object) -> str:
    return _REUSE_MODULE_LABELS.get(str(module or ""), str(module or "—"))


def _reuse_status_label(reuse_status: object, recompute_status: object) -> str:
    if recompute_status == "recomputed":
        return "🔄 已重算"
    if recompute_status == "recompute_failed":
        return "⚠ 重算失败(回退野生型)"
    if reuse_status == "reusable":
        return "✅ 可复用"
    if reuse_status == "stale":
        return "⚠ 已失效(待重算)"
    return "—"


def _render_modified_strain_shortlist(data: dict[str, object]) -> None:
    """改造后候选短名单（OE+KO · 分层复用 · ADR-004 #1 迭代2 D5）。

    复用同口径**野生型全基因组基线**的短名单，只对受改造影响的分泌层重算、其余复用（L1 打标 → L2 重算）。
    复用 `pichia_last_run_modifications` 暂存；未命中该口径的后台基线 → 诚实指引先跑后台构建。
    """
    mods = st.session_state.get("pichia_last_run_modifications")
    with st.expander("改造后候选短名单（OE + KO · 分层复用）", expanded=False):
        st.caption(
            "复用同口径**野生型全基因组基线**的 OE+KO 短名单，只对受改造影响的分泌层重算、其余直接复用——"
            "给改造后菌株一份带「可复用 / 已失效」标注的下一步候选。复用是近似：只对与瓶颈无关的**分泌专属层**"
            "（折叠/糖基化/转运等）干净有效，**代谢层保守重算**。相对信号、非绝对产量。"
        )
        if not isinstance(mods, dict) or not mods.get("target_id"):
            st.info("先在「仿真构建」页运行一次仿真，这里就能基于同一株算改造后候选短名单。")
            return
        if not mods.get("target_is_builtin"):
            st.info(
                "同上：改造后短名单要按目标 ID 复原野生型基线，三段式 / 自定义 JSON 目标没有内置定义可复原。"
                "改用**快速选择（内置模板）**跑一次即可使用。"
            )
            return

        oe_rx = [str(r) for r in (mods.get("oe_reaction_ids") or [])]
        ko_rx = [str(r) for r in (mods.get("ko_reaction_ids") or [])]
        if st.button("构建改造后候选短名单（复用地基 + 打标，约数十秒）", key="pichia_modified_shortlist_build_button"):
            with st.spinner("读野生型基线 + 解野生型/改造后瓶颈 + 受影响层打标…"):
                st.session_state["pichia_modified_shortlist_result"] = build_modified_strain_shortlist(
                    target_id=str(mods.get("target_id")),
                    ko_reaction_ids=tuple(ko_rx),
                    oe_reaction_ids=tuple(oe_rx),
                    mu=float(mods.get("mu") or 0.10),
                    media_type=int(mods.get("media_type") or 4),
                    carbon_source_id=str(mods.get("carbon_source_id") or "glucose"),
                    enable_ribosome_translation_constraint=bool(mods.get("enable_ribosome")),
                    enable_misfolding_constraint=bool(mods.get("enable_misfolding")),
                )
                st.session_state.pop("pichia_modified_shortlist_l2", None)  # 新 L1 → 清掉旧 L2

        readout = st.session_state.get("pichia_modified_shortlist_result")
        if not isinstance(readout, dict):
            return
        if not readout.get("available"):
            st.warning(
                "该口径下还没有**野生型全基因组后台基线**，无法复用。请先在此口径跑一次后台构建"
                "（`tools/run_genome_wide_ko_oe_screen_parallel.py`，hour-scale、跑完自动落口径指纹缓存），再回来构建短名单。"
            )
            return

        l2 = st.session_state.get("pichia_modified_shortlist_l2")
        shown = l2 if isinstance(l2, dict) else readout
        _render_modified_shortlist_result(shown)

        stale_total = int(readout.get("oe_stale_count", 0)) + int(readout.get("ko_stale_count", 0))
        if not isinstance(l2, dict) and stale_total:
            if st.button(
                f"重算已失效候选（{stale_total} 个，会额外求解）", key="pichia_modified_shortlist_l2_button"
            ):
                with st.spinner("在改造后菌株上重算已失效候选（OE 剂量响应 + KO 扰动）…"):
                    st.session_state["pichia_modified_shortlist_l2"] = recompute_stale_candidates(readout)
                st.rerun()


def _render_modified_shortlist_result(readout: dict[str, object]) -> None:
    layer = str(readout.get("layer") or "L1")
    c1, c2, c3 = st.columns(3)
    with c1:
        st.metric("层级", "L2（已重算失效项）" if layer == "L2" else "L1（即时复用）")
    with c2:
        st.metric("改造后瓶颈层", "、".join(_reuse_module_label(m) for m in (readout.get("modified_bottleneck_modules") or [])) or "—")
    with c3:
        st.metric("受影响层", "、".join(_reuse_module_label(m) for m in (readout.get("affected_modules") or [])) or "—")

    for title, key, reusable_key, stale_key in (
        ("OE 候选", "oe_candidates", "oe_reusable_count", "oe_stale_count"),
        ("KO 候选", "ko_candidates", "ko_reusable_count", "ko_stale_count"),
    ):
        candidates = [c for c in (readout.get(key) or []) if isinstance(c, dict)]
        st.markdown(f"**{title}**（可复用 {readout.get(reusable_key, 0)} · 已失效 {readout.get(stale_key, 0)}）")
        if not candidates:
            st.caption("（无有实质提升的候选。）")
            continue
        rows: list[dict[str, object]] = []
        for rank, candidate in enumerate(candidates, 1):
            effect = candidate.get("effective_effect", candidate.get("wildtype_effect"))
            rows.append(
                {
                    "排名": rank,
                    "候选": candidate.get("candidate") or candidate.get("gene_id") or "—",
                    "分泌层": _reuse_module_label(candidate.get("reuse_module")),
                    "复用状态": _reuse_status_label(candidate.get("reuse_status"), candidate.get("recompute_status")),
                    # 本模型的相对效应常是亚百分比（0.0x%）——用有效数字而非固定 2 位小数，否则全显示 0.00%、看不出排序。
                    "相对效应": f"{float(effect) * 100:.3g}%" if isinstance(effect, (int, float)) else "—",
                }
            )
        st.dataframe(pd.DataFrame(rows), width="stretch", hide_index=True)

    if readout.get("caveat"):
        st.caption(f"• {readout['caveat']}")
    for warning in readout.get("recompute_warnings") or []:
        st.warning(str(warning))


def _render_protein_cost_analysis(protein_cost: dict[str, object]) -> None:
    # Only rendered when 启用蛋白成本斜率对比 was on (pipeline.py sets protein_cost_analysis
    # to None otherwise) - everything shown here is LP-solve-derived, not a heuristic score.
    with st.expander("目标蛋白成本分析", expanded=True):
        st.caption(f"状态: {sim_result_value_label(protein_cost.get('result_status', 'draft_cost_slope_analysis'))}")
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
            st.warning(sim_result_warning_label(warning))


def _render_solver_robustness(solver_robustness: dict[str, object]) -> None:
    st.markdown("**求解器稳健性（瓶颈归因是否跨求解器稳定）**")
    st.caption(
        "对偶解在退化最优解处不唯一，换求解算法可能把瓶颈归到不同资源。"
        "ranking-insensitive-to-solver 表示归因稳定；ranking-sensitive-to-solver 表示归因是数值假象，不是生物学结论。"
    )
    classification = str(solver_robustness.get("classification", "—"))
    c1, c2, c3 = st.columns([2, 1, 1])
    with c1:
        st.metric("稳健性分类", sim_result_value_label(classification))
    with c2:
        st.metric("瓶颈跨求解器一致", "是" if solver_robustness.get("top_bottleneck_stable") else "否")
    with c3:
        st.metric("主导块跨求解器一致", "是" if solver_robustness.get("top_block_stable") else "否")
    if classification == "ranking-sensitive-to-solver":
        st.warning(str(solver_robustness.get("detail", "")))
    else:
        st.caption(str(solver_robustness.get("detail", "")))
    per_method = solver_robustness.get("per_method") or []
    if per_method:
        st.dataframe(
            _localized_frame(per_method, value_columns=("result_status", "top_dominant_block")),
            width='stretch',
            hide_index=True,
        )
    for warning in solver_robustness.get("warnings") or []:
        st.warning(sim_result_warning_label(warning))


def _short_reaction_label(reaction_id: str, max_len: int = 34) -> str:
    """图表轴/图例上的反应标签：优先显示**可读名称**，没有名称才退回压缩后的 id。

    名称解析统一走 `app.services.display_naming`（模型自带 rxnNames → 策展俗名 → 借基因名），
    别在这里再造一套——此前各页面各自截断 id，同一个反应在不同页面长得不一样、也永远看不到名字。
    """
    resolved = str(reaction_id)
    try:
        from app.services.display_naming import reaction_display_name

        name = reaction_display_name(resolved)
    except Exception:  # noqa: BLE001 - 命名是显示增强，失败就退回原 id
        name = ""
    if name:
        # 图表空间有限：有名字时只放名字（完整 id 在表格与悬浮里仍可见）。
        return name if len(name) <= max_len else f"{name[: max_len - 1]}…"
    label = resolved
    for suffix in ("_complex_formation", "_formation", "_complex"):
        if label.endswith(suffix):
            label = label[: -len(suffix)]
            break
    if len(label) > max_len:
        label = f"{label[: max_len - 12]}…{label[-11:]}"
    return label or resolved


def _resource_layer_label(process: object) -> str:
    """Chinese resource-layer label for an LP-attribution ``secretory_process`` code.

    The engine now classifies every reaction -- including the ``sec_*`` secretory-machine
    complexes (e.g. ``sec_Pdi1p_complex_formation`` -> ``disulfide_folding``) -- into the
    shared process vocabulary, and every one of those codes lives in the central i18n value
    dictionary, so this is a straight localized lookup. The earlier name-based inference
    fallback (which guessed the folding layer from the reaction id when the engine said
    ``unknown``) is no longer needed and was removed; unmapped codes still degrade to '未解析'.
    """
    code = str(process or "unknown")
    return sim_result_value_label(code)


def _oe_dose_response_curve_frame(oe_dose_response: dict[str, object]) -> pd.DataFrame:
    """Flatten the per-reaction dose-response points into a long frame for a line chart.

    One row per (reaction, factor): factor on x, relative gain (%) on y. The factor-1.0 baseline
    is included so each curve starts at 0%. Points with a missing factor/gain are dropped.
    """
    rows: list[dict[str, object]] = []
    for shape in oe_dose_response.get("reaction_shapes") or []:
        if not isinstance(shape, dict):
            continue
        reaction = _short_reaction_label(str(shape.get("reaction_id", "?")))
        shape_label = sim_result_value_label(shape.get("shape"))
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
            color_discrete_sequence=px.colors.qualitative.Set2,
            title="OE 剂量响应曲线：过表达倍数越高，分泌相对提升怎么走",
        )
        figure.update_layout(
            xaxis_title="过表达倍数（×，1 = 不过表达）",
            yaxis_title="分泌相对提升（%）",
            legend_title_text="反应｜形状",
        )
        st.plotly_chart(figure, width='stretch')

    def _pct(value: object) -> str:
        return f"{value * 100:.2f}%" if isinstance(value, (int, float)) else "—"

    shapes = [s for s in (oe_dose_response.get("reaction_shapes") or []) if isinstance(s, dict)]
    rows = [
        {
            "反应": _short_reaction_label(str(shape.get("reaction_id", "—"))),
            "形状": sim_result_value_label(shape.get("shape")),
            "最大相对增益": _pct(shape.get("max_relative_gain")),
            "半增益倍数(拐点)": (
                f"{float(shape['half_gain_factor']):g}×" if shape.get("half_gain_factor") is not None else "—"
            ),
            "最强剂量处相对增益": _pct(shape.get("relative_gain_at_max_factor")),
            "单调不减": sim_result_value_label(shape.get("monotonic_non_decreasing")),
        }
        for shape in shapes
    ]
    if rows:
        st.dataframe(pd.DataFrame(rows), width='stretch', hide_index=True)
    else:
        st.caption("（没有可用的剂量响应结果）")
    # A decrease that OE cannot truly cause is a degenerate-optimum artifact: flag it explicitly.
    for shape in shapes:
        if shape.get("shape") == "non_monotonic_numerical_artifact":
            st.warning(f"{shape.get('reaction_id')}：{shape.get('detail', '')}")
    for warning in oe_dose_response.get("warnings") or []:
        st.warning(sim_result_warning_label(warning))


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
        reaction_id = str(entry.get("reaction_id", "?"))
        rows.append(
            {
                "反应": _short_reaction_label(reaction_id),
                "限制强度": abs(float(marginal or 0.0)),
                "分泌资源层": _resource_layer_label(entry.get("secretory_process")),
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
        reaction_id = str(entry.get("reaction_id", "?"))
        rows.append(
            {
                "反应": _short_reaction_label(reaction_id),
                "限制强度": abs(float(marginal or 0.0)),
                "分泌资源层": _resource_layer_label(entry.get("secretory_process")),
            }
        )
    return pd.DataFrame(rows)


def _render_lp_attribution(lp_attribution: dict[str, object]) -> None:
    st.markdown("**约束归因证据（进阶：模型为什么这么排）**")
    st.caption("进阶证据：约束灵敏度草稿（“限制强度”＝该约束有多吃紧）；只看相对大小与名次，不要当精确数值用。")
    objective = lp_attribution.get("objective_evidence") if isinstance(lp_attribution.get("objective_evidence"), dict) else {}
    c1, c2, c3 = st.columns([1, 1, 1])
    with c1:
        st.metric("LP 归因状态", sim_result_value_label(lp_attribution.get("result_status")))
    with c2:
        st.metric("目标反应", objective.get("objective_reaction", "—") if isinstance(objective, dict) else "—")
    with c3:
        st.metric("分泌通量", objective.get("secretion_flux", "—") if isinstance(objective, dict) else "—")

    oe_actionable = lp_attribution.get("oe_actionable_bottlenecks") or []
    floor_only = lp_attribution.get("floor_constraints_not_oe_addressable") or []
    st.write("**过表达可缓解的瓶颈（卡在上限，按复合体）**")
    st.caption(
        "只列当前解里卡在上限的产能约束（生效约束）——这是过表达（放宽产能上限）真能缓解的瓶颈线索，按复合体给出。"
        "注意只是线索不是保证：耦合结构下放宽一个上限会让瓶颈转移，需与真实逐候选权衡交叉验证。"
    )
    if oe_actionable:
        bottleneck_frame = _lp_oe_bottleneck_frame(lp_attribution)
        if not bottleneck_frame.empty:
            figure = px.bar(
                bottleneck_frame.nlargest(8, "限制强度").sort_values("限制强度"),
                x="限制强度",
                y="反应",
                color="分泌资源层",
                orientation="h",
                text="限制强度",
                color_discrete_sequence=px.colors.qualitative.Set2,
                title="过表达可缓解的瓶颈：哪一层最卡分泌（限制强度越大越紧）",
            )
            figure.update_traces(texttemplate="%{text:.2g}", textposition="outside", cliponaxis=False)
            figure.update_layout(
                xaxis_title="限制强度（越大越卡分泌）",
                yaxis_title="",
                legend_title_text="分泌资源层",
                yaxis={"categoryorder": "total ascending"},
            )
            st.plotly_chart(figure, width='stretch')
        st.dataframe(
            _localized_frame(oe_actionable, value_columns=("bound_type", "secretory_process", "oe_actionable")),
            width='stretch',
            hide_index=True,
        )
    else:
        st.caption("（当前解没有卡在上限的产能约束）")
    if floor_only:
        st.write("**为什么受限：最强约束层（下界/最低需求，过表达松不动）**")
        st.caption(
            "这些下界（折叠 / 翻译 / ERAD 等最低需求）承载最大的限制强度，是本靶点「为什么受限、卡在哪一层」的答案；"
            "但过表达松的是上限、对它们无效，绝不能当过表达靶点（PDI1 单体、核糖体装配就是这里的经典假阳性）。"
        )
        floor_frame = _lp_floor_bottleneck_frame(lp_attribution)
        if not floor_frame.empty:
            figure = px.bar(
                floor_frame.nlargest(8, "限制强度").sort_values("限制强度"),
                x="限制强度",
                y="反应",
                color="分泌资源层",
                orientation="h",
                text="限制强度",
                color_discrete_sequence=px.colors.qualitative.Set2,
                title="为什么受限：最强约束层（限制强度越大越卡；过表达松不动）",
            )
            figure.update_traces(texttemplate="%{text:.0f}", textposition="outside", cliponaxis=False)
            figure.update_layout(
                xaxis_title="限制强度（越大越卡分泌）",
                yaxis_title="",
                legend_title_text="分泌资源层",
                yaxis={"categoryorder": "total ascending"},
            )
            st.plotly_chart(figure, width='stretch')
        st.dataframe(
            _localized_frame(floor_only, value_columns=("bound_type", "secretory_process", "oe_actionable")),
            width='stretch',
            hide_index=True,
        )

    sections = (
        ("主导约束块（按资源层汇总）", "dominant_constraint_blocks"),
        ("约束级限制强度 Top", "top_constraint_marginals"),
        ("边界级限制强度 Top（未按边界类型分离的原始表）", "top_bound_marginals"),
        ("目标相关通量", "target_related_fluxes"),
    )
    for title, key in sections:
        rows = lp_attribution.get(key) or []
        if rows:
            st.write(f"**{title}**")
            st.dataframe(
                _localized_frame(rows, value_columns=("constraint_type", "block", "bound_type", "secretory_process")),
                width='stretch',
                hide_index=True,
            )
    counts = lp_attribution.get("active_bound_counts")
    if isinstance(counts, dict) and counts:
        st.write("**生效边界限制强度计数**")
        st.dataframe(
            pd.DataFrame([{"项目": sim_result_column_label(key), "数量": value} for key, value in counts.items()]),
            width='stretch',
            hide_index=True,
        )
    for warning in lp_attribution.get("warnings") or []:
        st.warning(sim_result_warning_label(warning))


def _render_cost_slope_compatibility(cost_slope: dict[str, object]) -> None:
    st.markdown("**MATLAB 兼容的蛋白成本斜率对比（可选）**")
    st.caption(
        "当前默认路线是固定生长率、corrected 培养基、最大化目标蛋白分泌通量；"
        "历史 MATLAB 成本路线是固定目标分泌比例和生长率、优化葡萄糖摄取 Ex_glc_D，用于 Protein_cost_TP 定义对比。"
    )
    c1, c2, c3 = st.columns([1, 1, 1])
    with c1:
        st.metric("开启状态", sim_result_value_label(cost_slope.get("enabled", False)))
    with c2:
        st.metric("结果状态", sim_result_value_label(cost_slope.get("result_status")))
    with c3:
        st.metric("成功", sim_result_value_label(cost_slope.get("success")))
    _render_cost_slope_ratio_policy(cost_slope)
    st.caption(
        f"培养基兼容模式: {cost_slope.get('medium_compatibility_mode', 'corrected')}；"
        "该设置只影响可选 cost slope 对比，不改变默认分泌仿真。"
    )
    overrides = cost_slope.get("medium_bound_overrides") or []
    if overrides:
        st.write("**MATLAB 历史培养基边界覆盖**")
        st.dataframe(_localized_frame(overrides), width='stretch', hide_index=True)

    sections = (
        ("葡萄糖成本斜率", "glucose_cost_slopes"),
        ("核糖体成本斜率", "ribosome_cost_slopes"),
        ("成本斜率明细", "rows"),
    )
    for title, key in sections:
        rows = cost_slope.get(key) or []
        if rows:
            st.write(f"**{title}**")
            st.dataframe(
                _localized_frame(rows, value_columns=("status", "glucose_cost_status", "success", "cost_key")),
                width='stretch',
                hide_index=True,
            )
    comparison_scope = cost_slope.get("comparison_scope")
    if isinstance(comparison_scope, dict) and comparison_scope:
        with st.expander("对比定义", expanded=False):
            st.json(comparison_scope)
    for warning in cost_slope.get("warnings") or []:
        st.warning(sim_result_warning_label(warning))


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
        st.caption("在一小组生长速率上试算分泌量，看“长得快”和“分泌多”怎么权衡；不代表真实发酵生长预测。")
        best_flux = target_growth.get("best_secretion_point") or {}
        best_per_biomass = target_growth.get("best_secretion_per_biomass_point") or {}
        c1, c2, c3 = st.columns([1, 1, 1])
        with c1:
            # 这些枚举此前直出英文（decreasing / monotonic_decreasing_successful_grid），走中央字典汉化。
            st.metric("生长越快分泌怎么变", sim_result_value_label(target_growth.get("growth_sensitivity_label", "—")))
            reason = target_growth.get("growth_sensitivity_reason")
            if reason:
                st.caption(f"判断依据：{sim_result_value_label(reason)}")
        with c2:
            st.metric("分泌最高时的生长速率", best_flux.get("mu", "—") if isinstance(best_flux, dict) else "—")
        with c3:
            st.metric("单位菌体分泌最高时的生长速率", best_per_biomass.get("mu", "—") if isinstance(best_per_biomass, dict) else "—")

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
                width='stretch',
                hide_index=True,
            )
        for warning in target_growth.get("warnings") or []:
            st.warning(sim_result_warning_label(warning))


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
            # 单元格里的英文枚举码统一走集中字典翻译（表头由 display_columns 负责）
            for enum_column in (
                "recommendation_tier",
                "recommendation_label",
                "intervention_type",
                "execution_mode",
                "wet_lab_readiness",
                "standard_name_status",
                "model_gpr_executable",
                "oe_reaction_proxy",
            ):
                if enum_column in frame.columns:
                    frame[enum_column] = frame[enum_column].map(sim_result_value_label)
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
                "delta_objective": "Δ目标值",
                "secretory_process": "分泌环节",
                "database_annotation_sources": "数据库注释来源",
                "model_gpr_executable": "模型 GPR 可执行",
                "oe_reaction_proxy": "OE 反应代理",
                "evidence_tier": "证据等级",
                "recommendation_score": "推荐分",
                "rationale": "推荐理由",
            }
            columns = [key for key in display_columns if key in frame.columns]
            st.dataframe(frame[columns].rename(columns=display_columns), width='stretch', hide_index=True)
        else:
            st.info("当前候选没有进入模型内提升推荐。")
        st.caption("gene-level KO 与 reaction-level OE proxy 是不同证据层级；OE proxy 不能直接等同于湿实验基因过表达。")
        for warning in payload.get("warnings") or []:
            st.warning(sim_result_warning_label(warning))


__all__ = ["render_pichia_results"]
