"""全基因组 KO/OE 分泌-生长权衡筛查页面。

历史设计取舍见 pcSecYeastSpecies/docs/archive/pichia_ko_oe_genome_screen_design_2026-07-02.md。
任何人都能直接打开这个页面触发筛查；筛查本身是小时级的后台任务，跑在独立
子进程里（不是 Streamlit 内部线程），这样刷新/重启页面不会杀掉正在跑的任务。
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import plotly.express as px
import streamlit as st

from app.core.i18n import sim_result_column_label, sim_result_value_label, sim_result_warning_label
from app.services import genome_wide_screen_analysis as analysis
from app.services import genome_wide_screen_service as service
from app.services import genome_wide_screen_shortlist as shortlist_service
from app.services import screen_report_service
from app.services.genome_wide_screen_registry import (
    RunInfo,
    latest_runs_by_group,
    list_runs,
    older_runs_by_group,
    run_scope_family,
)
from app.services.pichia_secretion_service import discover_project_paths
from app.ui.views.homology_audit import apply_homology_audit_prefill
from app.ui.views.simulation import apply_simulation_prefill

DEFAULT_TARGETS = ["hLF", "OPN_ALPHA_FULL_PROJECT"]
QUEUE_STATE_KEY = "genome_wide_screen_queued_request"
CONFLICT_STATE_KEY = "genome_wide_screen_pending_conflict"


def _paths():
    return discover_project_paths(Path(__file__))


def render_genome_wide_screen() -> None:
    st.header("全基因组 KO/OE 分泌-生长权衡筛查")
    st.caption(
        "对模型里全部基因逐一评估敲除(KO)/过表达(OE)对目标蛋白分泌产量和细胞生长速率的影响。"
        "全量跑一次是小时级任务，启动后可以关闭页面，回来再看进度。"
    )

    paths = _paths()
    _maybe_launch_queued_request(paths)

    runs = list_runs(paths)
    active_runs = [run for run in runs if run.status in {"starting", "running"}]

    # Results first: this is an expensive, infrequent, hours-long computation -
    # once it's done, checking its output is the common case on every visit to
    # this page, while starting a new run and eyeballing old run history are
    # both rarer, so they sit lower and out of the way.
    _render_results_section(paths, runs)
    st.divider()
    _render_launch_controls(paths, active_runs)
    st.divider()
    _render_run_list(runs)


SCOPE_LABELS = {
    "gene": "全基因组（约1025个基因，小时级）",
    "gene_limited": "小规模基因试跑（smoke/pilot，分钟级）",
    "catalog": "策展复合体反应对照表（约30个反应，分钟级）",
    "complex_hypothesis": "复合体假设性整体过表达测试（源自已有KO筛查结果，分钟级）",
}


def _render_launch_controls(paths, active_runs: list[RunInfo]) -> None:
    st.subheader("启动新一轮筛查")
    scope = st.radio(
        "筛查范围",
        options=["gene", "catalog"],
        format_func=lambda value: SCOPE_LABELS[value],
        horizontal=True,
        key="genome_wide_screen_scope",
        help=(
            "全基因组：逐一评估模型里所有基因的KO/OE，覆盖面完整，但是小时级任务。"
            "策展复合体反应对照表：只测文献策展、有实验证据支持的约30个复合体级反应"
            "（KAR2/PDI1/PMT等），分钟级出结果，适合先看一轮已知靶点有没有生长代价。"
        ),
    )
    col_targets, col_mode, col_workers = st.columns(3)
    targets = col_targets.multiselect("靶点", options=DEFAULT_TARGETS, default=DEFAULT_TARGETS)
    mode = col_mode.selectbox("精度模式", options=["fast", "precise"], index=0, help="fast=3个生长速率采样点，precise=11个点，耗时更长")
    workers = col_workers.number_input("并行worker数", min_value=1, max_value=8, value=service.DEFAULT_WORKERS)
    gene_limit = None
    if scope == "gene":
        gene_limit_enabled = st.checkbox("只跑前N个基因（用于小规模验证，不勾选=全基因组）")
        gene_limit = st.number_input("N", min_value=1, value=50, disabled=not gene_limit_enabled) if gene_limit_enabled else None

    if st.button("启动筛查", type="primary", disabled=not targets):
        conflicts = service.check_for_conflicts(paths)
        request = {"targets": targets, "mode": mode, "workers": int(workers), "gene_limit": gene_limit, "scope": scope}
        if conflicts:
            st.session_state[CONFLICT_STATE_KEY] = request
        else:
            result = service.submit_screen(**request, paths=paths)
            st.success(f"已启动：{result.run_name}。可以离开页面，回来刷新查看进度。")
            st.rerun()

    pending = st.session_state.get(CONFLICT_STATE_KEY)
    if pending:
        _render_conflict_dialog(paths, pending, active_runs)


def _render_conflict_dialog(paths, pending: dict, active_runs: list[RunInfo]) -> None:
    st.warning(
        f"检测到 {len(active_runs)} 个正在运行的任务："
        + "、".join(f"{run.run_name}（{run.progress_label}）" for run in active_runs)
    )
    choice = st.radio(
        "新任务怎么处理？",
        ["取消，不启动", "排队，等当前任务跑完自动开始", "强制并发启动（会更慢，抢同一批CPU核）"],
        key="genome_wide_screen_conflict_choice",
    )
    if st.button("确认"):
        if choice.startswith("取消"):
            st.session_state.pop(CONFLICT_STATE_KEY, None)
            st.rerun()
        elif choice.startswith("排队"):
            st.session_state[QUEUE_STATE_KEY] = pending
            st.session_state.pop(CONFLICT_STATE_KEY, None)
            st.info("已加入队列，当前任务跑完、且没有其他运行中任务时会自动启动。")
            st.rerun()
        else:
            result = service.submit_screen(**pending, paths=paths)
            st.session_state.pop(CONFLICT_STATE_KEY, None)
            st.success(f"已并发启动：{result.run_name}")
            st.rerun()


def _maybe_launch_queued_request(paths) -> None:
    queued = st.session_state.get(QUEUE_STATE_KEY)
    if not queued:
        return
    if service.check_for_conflicts(paths):
        return  # still busy; try again next time this page loads
    result = service.submit_screen(**queued, paths=paths)
    st.session_state.pop(QUEUE_STATE_KEY, None)
    st.toast(f"排队的任务已自动启动：{result.run_name}")


def _run_table_frame(runs: list[RunInfo]) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "运行名": run.run_name,
                "状态": sim_result_value_label(run.status),
                "进度": run.progress_label,
                "靶点": ", ".join(run.targets),
                "模式": sim_result_value_label(run.mode),
                "更新时间": run.updated_at,
                "心跳超时未更新": sim_result_value_label(run.is_stale),
            }
            for run in runs
        ]
    )


def _render_run_list(runs: list[RunInfo]) -> None:
    st.subheader("运行记录")
    if not runs:
        st.caption("还没有跑过。")
        return

    latest = latest_runs_by_group(runs)
    older = older_runs_by_group(runs)

    # Grouped by analysis type (matching SCOPE_LABELS), then by target within each type -
    # researchers looking for "the" result of one analysis type shouldn't have to pick it
    # out of a flat list that also contains superseded re-runs of the same thing.
    for scope, label in SCOPE_LABELS.items():
        scope_runs = [run for run in latest if run_scope_family(run) == scope]
        if not scope_runs:
            continue
        st.markdown(f"**{label}**")
        scope_runs = sorted(scope_runs, key=lambda run: run.targets)
        st.dataframe(_run_table_frame(scope_runs), width='stretch', hide_index=True)

    unlabeled = [run for run in latest if run_scope_family(run) not in SCOPE_LABELS]
    if unlabeled:
        st.markdown("**其他**")
        st.dataframe(_run_table_frame(unlabeled), width='stretch', hide_index=True)

    if older:
        with st.expander(f"查看历史版本（{len(older)}个已被更新的运行取代，仍保留在磁盘上）"):
            st.dataframe(_run_table_frame(older), width='stretch', hide_index=True)


def _split_result_runs(runs: list[RunInfo]) -> tuple[list[RunInfo], list[RunInfo]]:
    latest_runs = latest_runs_by_group(runs)
    latest_run_names = {run.run_name for run in latest_runs}
    latest_done = [run for run in latest_runs if run.status == "done"]
    older_done = [run for run in runs if run.status == "done" and run.run_name not in latest_run_names]
    return latest_done, older_done


def _run_option_label(run: RunInfo) -> str:
    return f"{run.run_name} · {SCOPE_LABELS.get(run_scope_family(run), run.scope)}"


def _render_results_section(paths, runs: list[RunInfo]) -> None:
    st.subheader("结果查看")
    done_runs = [run for run in runs if run.status == "done"]
    if not done_runs:
        st.caption("还没有已完成的运行。")
        return

    latest_done, older_done = _split_result_runs(runs)

    # No dropdown, and no single flat list either: "latest" naturally has one
    # entry per scope family (gene/catalog/complex_hypothesis/...) x target,
    # e.g. overnight_hLF_full and overnight_OPN_full are both "the current
    # gene-scope result" at once. Picking "which analysis" and "which target
    # within it" as two small decisions stays readable; one flat radio with
    # every combination's full label does not.
    selected_run: RunInfo | None = None
    if latest_done:
        by_scope: dict[str, list[RunInfo]] = {}
        for run in latest_done:
            by_scope.setdefault(run_scope_family(run), []).append(run)
        scope_options = [scope for scope in SCOPE_LABELS if scope in by_scope]
        scope_options += [scope for scope in by_scope if scope not in SCOPE_LABELS]

        if len(scope_options) == 1:
            selected_scope = scope_options[0]
        else:
            selected_scope = st.radio(
                "查看哪类分析的结果",
                options=scope_options,
                format_func=lambda scope: SCOPE_LABELS.get(scope, scope),
                horizontal=True,
                key="genome_wide_screen_result_scope_select",
            )

        scope_runs = sorted(by_scope[selected_scope], key=lambda run: run.targets)
        if len(scope_runs) == 1:
            selected_run = scope_runs[0]
            st.caption(f"当前显示：{selected_run.run_name}（靶点：{', '.join(selected_run.targets)}）")
        else:
            selected_run = st.radio(
                "选择靶点",
                options=scope_runs,
                format_func=lambda run: f"{run.run_name}（靶点：{', '.join(run.targets)}）",
                horizontal=True,
                key="genome_wide_screen_result_run_select",
            )
    else:
        st.caption("当前没有最新的已完成运行；可在下方“历史版本”中查看已被新运行取代的完成结果。")

    if older_done:
        with st.expander(f"改看一个已被取代的历史版本（{len(older_done)} 个，仍保留在磁盘上）"):
            viewing_older = st.radio(
                "历史版本",
                options=[None, *older_done],
                format_func=lambda run: "（不查看，用上面选中的最新结果）" if run is None else _run_option_label(run),
                key="genome_wide_screen_older_run_select",
            )
            if viewing_older is not None:
                selected_run = viewing_older
                st.info(f"当前显示历史版本结果：{selected_run.run_name}（不是最新结果，可在此处切回“不查看”）")

    if selected_run is None:
        return

    run_name = selected_run.run_name
    # Older runs recorded before catalog scope existed always used gene_tradeoff_rows.csv;
    # csv_path in status.json is authoritative when present (it's scope-specific).
    csv_path = Path(selected_run.csv_path) if selected_run.csv_path else paths.local_runs_dir / run_name / "gene_tradeoff_rows.csv"
    if not csv_path.exists():
        st.error(f"找不到结果文件：{csv_path}")
        return

    frame = analysis.load_gene_tradeoff_csv(str(csv_path))
    target_ids = sorted(frame.target_id.dropna().unique().tolist())
    per_target_results = {target_id: analysis.analyze_single_target(frame, target_id) for target_id in target_ids}

    meta_cols = st.columns(4)
    meta_cols[0].metric("精度模式", selected_run.mode)
    meta_cols[1].metric("更新时间", selected_run.updated_at or "—")
    meta_cols[2].metric("靶点数", len(target_ids))
    meta_cols[3].metric("候选行数", len(frame))

    # R1 瓶颈读出是 per-target、单独计算并缓存的资产（不随本次筛查产生）；有就显示“为什么受限”，没有优雅降级。
    r1_readout_dir = paths.local_runs_dir / "r1_readout"
    # R2 剂量响应形状也是离线/后台单独扫描并缓存的（有界额外求解），面板只读缓存；缺失显示“未扫描”。
    dose_response_dir = paths.local_runs_dir / "candidate_shortlist_readout"
    tabs = st.tabs([f"靶点：{target_id}" for target_id in target_ids] + (["靶点间差异"] if len(target_ids) >= 2 else []))
    for tab, target_id in zip(tabs, target_ids):
        with tab:
            try:
                readout = shortlist_service.build_shortlist_readout(
                    frame, target_id, r1_readout_dir=r1_readout_dir, dose_response_dir=dose_response_dir
                )
                _render_shortlist_readout(readout, target_id=target_id)
                st.divider()
            except Exception as exc:  # noqa: BLE001 - 读出只是概览，任何失败都不该拖垮下方的明细表格
                st.caption(f"候选短名单读出暂不可用：{exc}")
            _render_dimension_tables(per_target_results[target_id])
    if len(target_ids) >= 2:
        with tabs[-1]:
            divergence = analysis.analyze_target_divergence(frame, target_ids)
            st.markdown("同一个基因的KO，在不同靶点上效应差异最大的候选：")
            st.dataframe(_localize_screen_frame(divergence), width='stretch', hide_index=True)

    st.divider()
    _render_llm_report_section(paths, selected_run, csv_path)


def _render_shortlist_readout(readout: dict, *, target_id: str) -> None:
    """候选短名单读出面板：为什么受限(R1) + OE 提升候选短名单(图表) + 该测什么(R4)。

    复用本次筛查的缓存解 + 单独缓存的 R1 读出，零新增求解；相对信号，非绝对产量。
    """
    shortlist = readout.get("oe_shortlist") or []
    floors = readout.get("why_limited_floors") or []
    voi = readout.get("value_of_information") or {}

    with st.expander("候选短名单读出：为什么受限 · OE 提升候选 · 该测什么", expanded=True):
        # 一句话结论
        if shortlist and readout.get("has_strong_oe_lever"):
            top = shortlist[0]
            headline = f"OE 提升候选里 **{top['candidate']}**（{top['layer']}）最强（+{float(top['effect']) * 100:.2f}%）"
        else:
            headline = "**没有强 OE 提升杠杆**（最高相对提升 < 1%）——这个靶点大概率不受限于可 OE 的分泌机器上限"
        top_floor = floors[0]["reaction_id"] if floors else None
        if top_floor:
            st.markdown(f"> 一句话：**{target_id}** 的分泌最强约束在 `{top_floor}`（下界/最低要求，OE 动不了）；{headline}。")
        else:
            st.markdown(f"> 一句话：**{target_id}** — {headline}。")
        st.caption("相对信号，非绝对产量 / mg·L⁻¹；复用本次筛查在固定倍数 OE、corrected 培养基下的模型解，零新增求解。")

        # 1. 为什么受限（R1）
        st.markdown("**为什么受限（R1 LP 影子价格 · 最强约束层）**")
        if floors:
            st.caption(
                "下界=最低要求类约束，承载最大影子价格，是“卡在哪一层”的答案；但 OE 放宽的是上限、对它们无效"
                "（floor≠可 OE 杠杆）。此 R1 读出对该靶点单独计算、单独缓存，不是本次筛查的产物。"
            )
            floor_frame = pd.DataFrame(
                [{"反应": _short_reaction(f["reaction_id"]), "影子价格(绝对值)": float(f["abs_marginal"])} for f in floors]
            )
            figure = px.bar(
                floor_frame.sort_values("影子价格(绝对值)"),
                x="影子价格(绝对值)",
                y="反应",
                orientation="h",
                text="影子价格(绝对值)",
                color_discrete_sequence=["#0F766E"],
                title="最强约束层：影子价格绝对值越大越限制分泌（下界，OE 动不了）",
            )
            figure.update_traces(texttemplate="%{text:.3g}", textposition="outside", cliponaxis=False)
            figure.update_layout(xaxis_title="影子价格绝对值（越大越限制分泌）", yaxis_title="", yaxis={"categoryorder": "total ascending"})
            st.plotly_chart(figure, width='stretch')
        else:
            st.caption(
                "未找到该靶点的 R1 瓶颈读出（local_runs/r1_readout/）。"
                "可先用 python_pichia/tools/run_target_bottleneck_lp_attribution_check.py 生成，再回来看“为什么受限”。"
            )

        # 2. OE 提升候选短名单
        st.markdown(f"**OE 提升候选短名单（按相对提升排序，top-{len(shortlist)}）**")
        if shortlist:
            dose_available = bool(readout.get("dose_response_available"))
            table_rows = []
            for row in shortlist:
                entry = {
                    "候选": str(row["candidate"]),
                    "资源层": str(row["layer"]),
                    "相对提升(%)": float(row["effect"]) * 100.0,
                    "生长保持": round(float(row["growth_retention"]), 3),
                    "证据置信度": sim_result_value_label(row["confidence"]) if row["confidence"] else "—",
                }
                if dose_available:
                    # 越加越好(线性) / 很快到顶(饱和) / 要过阈值 / 无响应——来自离线剂量响应扫描缓存
                    entry["剂量响应形状"] = sim_result_value_label(row.get("shape")) if row.get("shape") else "—"
                    max_gain = row.get("shape_max_gain")
                    entry["最大增益"] = f"+{float(max_gain) * 100:.2f}%" if isinstance(max_gain, (int, float)) else "—"
                    half = row.get("shape_half_gain_factor")
                    # 半增益倍数：达到最大增益一半所需的 OE 倍数，越小说明越早到顶（饱和越快）
                    entry["半增益倍数"] = f"{float(half):.2g}×" if isinstance(half, (int, float)) else "—"
                table_rows.append(entry)
            shortlist_frame = pd.DataFrame(table_rows)
            figure = px.bar(
                shortlist_frame.sort_values("相对提升(%)"),
                x="相对提升(%)",
                y="候选",
                color="资源层",
                orientation="h",
                text="相对提升(%)",
                color_discrete_sequence=px.colors.qualitative.Set2,
                title="OE 提升候选：相对野生型的分泌提升（越长越强，颜色=分泌资源层）",
            )
            figure.update_traces(texttemplate="%{text:.2f}%", textposition="outside", cliponaxis=False)
            figure.update_layout(
                xaxis_title="相对提升（%，相对野生型；非绝对产量）",
                yaxis_title="",
                legend_title_text="分泌资源层",
                yaxis={"categoryorder": "total ascending"},
            )
            st.plotly_chart(figure, width='stretch')
            st.dataframe(shortlist_frame, width="stretch", hide_index=True)
            risky = readout.get("growth_risky_candidates") or []
            if risky:
                st.caption("⚠️ 生长有代价（保持率 < 0.9）：" + "、".join(str(c) for c in risky))
            _render_dose_response_note(readout, dose_available)
            if dose_available:
                _render_dose_response_curve(readout, shortlist)
        else:
            st.caption("本次筛查在当前 OE 倍数下没有 ratio>1 的 OE 提升候选（不代表无解，可能是单基因、当前倍数强度不足以突破瓶颈）。")

        # 3. 该测什么（R4 价值-of-information）
        st.markdown("**该测什么（R4 价值-of-information）**")
        st.caption("模型给的是相对排序、不是绝对产量。这里标出顶部名次里模型分不清的候选，并给出最能消解歧义的最小湿实验——只排测量优先级，不预测结果。")
        _render_shortlist_voi(voi)


def _render_shortlist_voi(voi: dict) -> None:
    items = [item for item in (voi.get("information_items") or []) if isinstance(item, dict)]
    if voi.get("has_actionable_ambiguity") and items:
        st.warning(f"顶部排序有 {len(items)} 处近似并列（模型分不清谁更好），当前名次不完全可信——建议按下表做最小湿实验定序。")
        rows = []
        for item in items:
            candidate_text = "、".join(str(candidate) for candidate in (item.get("candidates") or []))
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
        st.dataframe(pd.DataFrame(rows), width="stretch", hide_index=True)
    elif voi.get("ranked_candidates"):
        st.success("顶部候选相对提升分离明显，当前排序较可信；优先验证榜首即可。")
    else:
        st.caption("无可排序的 OE 提升候选，暂无“该测什么”的明确建议。")
    for warning in voi.get("warnings") or []:
        st.caption(sim_result_warning_label(warning))


def _render_dose_response_note(readout: dict, dose_available: bool) -> None:
    if not dose_available:
        st.caption(
            "剂量响应形状：本靶点短名单尚未扫描——这一步要对每个候选重解多个过表达倍数（有界额外求解），"
            "离线/后台用 `run_shortlist_dose_response.py` 生成后，上表会多出“剂量响应形状”一列。当前先看相对提升排序即可。"
        )
        return
    dose = readout.get("dose_response") or {}
    factors = dose.get("tested_factors") or []
    factors_txt = "、".join(f"{float(f):g}×" for f in factors) or "—"
    st.caption(
        f"剂量响应形状（扫描倍数 {factors_txt}；相对形状信号、非绝对产量）：饱和型＝适度过表达就够、再加收益递减；"
        "线性型＝还在涨、值得加大表达；阈值型＝要超过某个最小倍数才起效；无响应＝任何倍数都几乎没提升、别在这过表达。"
    )
    for warning in dose.get("warnings") or []:
        st.caption(sim_result_warning_label(warning))


# 曲线只画能看出趋势的杠杆（最大增益≥此值），近乎平的小候选会糊在 x 轴附近、只添乱。
_DOSE_RESPONSE_CURVE_MIN_GAIN = 0.01  # 1%
_DOSE_RESPONSE_CURVE_MAX_LINES = 6


def _render_dose_response_curve(readout: dict, shortlist: list) -> None:
    """把短名单头部候选的剂量响应画成曲线：过表达倍数 vs 分泌相对提升(%)，曲线变平=很快到顶。"""
    dose = readout.get("dose_response") or {}
    shapes = [s for s in (dose.get("reaction_shapes") or []) if isinstance(s, dict)]
    if not shapes:
        return
    reaction_to_candidate = {str(row["reaction"]): str(row["candidate"]) for row in shortlist}
    ranked = sorted(shapes, key=lambda s: -(float(s.get("max_relative_gain") or 0.0)))
    rows: list[dict[str, object]] = []
    plotted = 0
    for shape in ranked:
        max_gain = shape.get("max_relative_gain")
        if not isinstance(max_gain, (int, float)) or max_gain < _DOSE_RESPONSE_CURVE_MIN_GAIN:
            continue
        if plotted >= _DOSE_RESPONSE_CURVE_MAX_LINES:
            break
        reaction = str(shape.get("reaction_id"))
        candidate = reaction_to_candidate.get(reaction, reaction)
        for point in shape.get("point_deltas") or []:
            if not isinstance(point, dict):
                continue
            factor = point.get("factor")
            gain = point.get("relative_gain")
            if factor is None or gain is None:
                continue
            rows.append(
                {"过表达倍数": float(factor), "分泌相对提升(%)": float(gain) * 100.0, "候选": candidate}
            )
        plotted += 1
    if not rows:
        return
    curve = pd.DataFrame(rows)
    figure = px.line(
        curve.sort_values("过表达倍数"),
        x="过表达倍数",
        y="分泌相对提升(%)",
        color="候选",
        markers=True,
        color_discrete_sequence=px.colors.qualitative.Set2,
        title="剂量响应曲线：过表达倍数越高分泌怎么走（曲线变平＝很快到顶、再加收益递减）",
    )
    figure.update_layout(
        xaxis_title="OE 过表达倍数（×，1＝不过表达）",
        yaxis_title="分泌相对提升（%，相对野生型；非绝对产量）",
        legend_title_text=f"候选（仅显示最大增益≥{_DOSE_RESPONSE_CURVE_MIN_GAIN:.0%} 的前 {_DOSE_RESPONSE_CURVE_MAX_LINES} 个）",
    )
    st.plotly_chart(figure, width='stretch')


_SCREEN_VALUE_COLUMNS = (
    "candidate_kind",
    "gpr_role",
    "intervention_type",
    "feasibility_interpretation",
    "annotation_confidence",
    "standard_name_status",
    "has_timeout",
    "secretory_process",
)


def _localize_screen_frame(df: pd.DataFrame) -> pd.DataFrame:
    """筛查表的展示用中文副本：枚举列映射成中文、列名走中央 i18n 字典。不改原 df——
    行选择逻辑仍按原始英文列名（`gene_id` 等）读，见 _render_verifiable_table。"""
    if len(df.columns) == 0:
        return df  # 真正无列的空 df 无需处理；有列的空 df 仍要汉化表头
    display = df.copy()
    for column in _SCREEN_VALUE_COLUMNS:
        if column in display.columns:
            display[column] = display[column].map(sim_result_value_label)
    return display.rename(columns={column: sim_result_column_label(column) for column in display.columns})


def _short_reaction(reaction_id: object) -> str:
    """去掉复合体反应名冗长的样板后缀，缩短图表 y 轴标签（模型自身实体名，非保密内容）。"""
    label = str(reaction_id)
    for suffix in ("_complex_formation", "_complex", "_formation"):
        if label.endswith(suffix):
            label = label[: -len(suffix)]
            break
    return label if len(label) <= 34 else label[:31] + "…"


def _render_dimension_tables(result: analysis.DimensionalResults) -> None:
    st.metric("必需基因数（KO后任何测试生长速率下都不可行）", len(result.essential_genes))
    with st.expander(f"必需基因清单（{len(result.essential_genes)}）"):
        st.dataframe(_localize_screen_frame(result.essential_genes), width='stretch', hide_index=True)

    if not result.solver_inconclusive_rows.empty:
        st.warning(
            f"有 {len(result.solver_inconclusive_rows)} 个候选行因为求解超时或求解器失败而无法证明可行性；"
            "这些不是必需基因结论。"
        )
    with st.expander(f"求解未定的候选行（{len(result.solver_inconclusive_rows)}）", expanded=False):
        st.dataframe(_localize_screen_frame(result.solver_inconclusive_rows), width='stretch', hide_index=True)

    with st.expander(f"求解器重试证据（{len(result.solver_retry_evidence)}）", expanded=False):
        st.dataframe(_localize_screen_frame(result.solver_retry_evidence), width='stretch', hide_index=True)

    st.markdown(f"**产量升高但生长受损的KO候选（{len(result.ko_yield_up_growth_cost)}）** — 需要额外生物学方法补救才具备实验可行性")
    _render_verifiable_table(
        result.ko_yield_up_growth_cost,
        target_id=result.target_id,
        intervention_type="KO",
        table_key=f"dimtable_ko_cost_{result.target_id}",
    )

    st.markdown(f"**产量升高且生长完全不受影响的KO候选（{len(result.ko_clean_wins)}）** — 零代价候选")
    _render_verifiable_table(
        result.ko_clean_wins,
        target_id=result.target_id,
        intervention_type="KO",
        table_key=f"dimtable_ko_clean_{result.target_id}",
    )

    with st.expander(f"KO降低产量的候选（{len(result.ko_yield_down)}）— 可行但比野生型更差", expanded=False):
        st.caption(
            "这些基因敲除后细胞还能活（不是必需基因），但分泌产量反而降低。"
            "很多是`gpr_role=complex_subunit`——敲掉多亚基复合体的一个亚基削弱了整体功能，"
            "间接拖累和分泌共享的资源预算。这类基因目前只有KO方向的数据："
            "复合体亚基默认不支持单基因OE（避免在没有验证亚基化学计量/限速步骤证据的情况下虚构容量提升），"
            "所以看不到对应的OE候选——不是没测，是模型判断单基因OE不可信。"
        )
        _render_verifiable_table(
            result.ko_yield_down,
            target_id=result.target_id,
            intervention_type="KO",
            table_key=f"dimtable_ko_down_{result.target_id}",
        )

    st.markdown(f"**产量升高的OE候选（{len(result.oe_yield_up)}）**")
    if result.oe_yield_up.empty:
        st.caption("当前OE倍数下没有找到能提升产量的基因——不代表代码有问题，可能是单基因、当前倍数的过表达强度不足以突破瓶颈。")
    _render_verifiable_table(
        result.oe_yield_up,
        target_id=result.target_id,
        intervention_type="OE",
        table_key=f"dimtable_oe_up_{result.target_id}",
    )

    if not result.complex_oe_hypothesis.empty:
        st.markdown(f"**假设性整体过表达测试（{len(result.complex_oe_hypothesis)}）** — 针对KO降低分泌、但因是复合体亚基而没有单基因OE数据的候选")
        st.warning(
            "**这是假设性结果，不是已验证的过表达方案：** " + str(result.complex_oe_hypothesis.iloc[0]["hypothesis_note"])
        )
        _render_verifiable_table(
            result.complex_oe_hypothesis.drop(columns=["hypothesis_note"]),
            target_id=result.target_id,
            intervention_type="OE",
            table_key=f"dimtable_complex_hypothesis_{result.target_id}",
        )


def _render_verifiable_table(df: pd.DataFrame, *, target_id: str, intervention_type: str, table_key: str) -> None:
    """Render a candidate table with row selection wired to the 仿真验证 cross-link.

    Selecting a row surfaces a button that pre-fills the simulation page's
    target + KO/OE gene inputs with this exact candidate and jumps there, so a
    reviewer can go from "this gene looked interesting" to a live simulation
    in one click instead of re-typing the target/gene IDs by hand.
    """
    if df.empty:
        st.dataframe(_localize_screen_frame(df), width='stretch', hide_index=True)
        return
    event = st.dataframe(
        _localize_screen_frame(df),  # 展示用中文副本
        width='stretch',
        hide_index=True,
        on_select="rerun",
        selection_mode="single-row",
        key=table_key,
    )
    selected_rows = event.selection.rows if event and event.selection else []
    if selected_rows:
        row = df.iloc[selected_rows[0]]  # 行选择按原始英文列 df、按位置对齐（展示副本行序一致）
        candidate_id = str(row["gene_id"])
        candidate_kind = str(row["candidate_kind"]) if "candidate_kind" in row else "gene"
        common_name = str(row["common_name"]) if "common_name" in row and row["common_name"] else ""
        standard_symbol = str(row.get("standard_symbol") or "").strip()
        gene_display_name = str(row.get("gene_display_name") or "").strip()
        friendly_name = standard_symbol or gene_display_name or common_name
        display_label = f"{friendly_name} ({candidate_id})" if friendly_name else candidate_id
        is_gene = candidate_kind == "gene"
        columns = st.columns(2 if is_gene else 1)
        if columns[0].button(
            f"在仿真验证中核实 {display_label}",
            key=f"{table_key}_verify_btn",
            type="primary",
        ):
            apply_simulation_prefill(target_id, candidate_id, intervention_type, candidate_kind)
        if is_gene:
            if columns[1].button(f"查看 {display_label} 的同源证据", key=f"{table_key}_homology_btn"):
                apply_homology_audit_prefill(candidate_id)


def _render_llm_report_section(paths, selected_run: RunInfo, csv_path: Path) -> None:
    st.subheader("研发建议报告（Writer LLM + Judge LLM）")
    st.caption(
        "页面只读取已有筛查结果生成 fact pack；点击按钮后才调用 LLM。"
        "LLM 只能读取 fact pack，最终报告必须通过程序校验和 Judge 审核。"
    )
    try:
        fact_pack = screen_report_service.build_fact_pack_for_runs(paths, csv_paths=(csv_path,))
        fact_summary = screen_report_service.summarize_fact_pack(fact_pack)
    except Exception as exc:  # noqa: BLE001 - user-facing diagnostic
        st.error(f"fact pack 生成失败：{exc}")
        return

    target_summary = fact_summary.get("targets", {})
    st.dataframe(
        _localize_screen_frame(
            pd.DataFrame(
                [
                    {"靶点": target, **counts}
                    for target, counts in target_summary.items()
                ]
            )
        ),
        width="stretch",
        hide_index=True,
    )
    st.download_button(
        "下载 fact_pack.json",
        data=json.dumps(fact_pack, ensure_ascii=False, indent=2),
        file_name=f"{selected_run.run_name}_fact_pack.json",
        mime="application/json",
    )

    latest_report_dirs = screen_report_service.latest_report_runs(paths, run_name=selected_run.run_name)
    latest_final = next((path / "final_report.md" for path in latest_report_dirs if (path / "final_report.md").exists()), None)
    if latest_final:
        with st.expander(f"查看最近一次已通过审核的报告：{latest_final.parent.name}", expanded=False):
            report_text = latest_final.read_text(encoding="utf-8")
            st.markdown(report_text)
            st.download_button(
                "下载 final_report.md",
                data=report_text,
                file_name="final_report.md",
                mime="text/markdown",
            )

    if st.button("生成研发建议报告", type="primary", key=f"screen_report_generate_{selected_run.run_name}"):
        with st.spinner("正在生成 fact pack、调用 Writer/Judge 并校验……"):
            try:
                result = screen_report_service.generate_judged_screen_report(paths, csv_paths=(csv_path,))
            except Exception as exc:  # noqa: BLE001 - missing LLM config should be visible, not fatal
                st.error(f"报告生成未完成：{exc}")
                st.info("如未配置 OPENAI_API_KEY，可先下载 fact_pack.json；测试环境使用 fake LLM，不依赖真实 API。")
                return
        if result.success and result.final_report_path:
            st.success(f"报告已生成：{result.final_report_path}")
            st.markdown(result.final_report_path.read_text(encoding="utf-8"))
        else:
            st.error("报告未通过程序校验或 Judge 审核，未生成 final_report.md。")
            st.json({"validator": result.validator_result, "judge": result.judge_result, "manifest": str(result.manifest_path)})


__all__ = ["render_genome_wide_screen"]
