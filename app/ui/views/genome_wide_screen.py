"""全基因组 KO/OE 分泌-生长权衡筛查页面。

设计取舍见 pcSecYeastSpecies/docs/pichia_ko_oe_genome_screen_design.md。
任何人都能直接打开这个页面触发筛查；筛查本身是小时级的后台任务，跑在独立
子进程里（不是 Streamlit 内部线程），这样刷新/重启页面不会杀掉正在跑的任务。
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import streamlit as st

from app.services import genome_wide_screen_analysis as analysis
from app.services import genome_wide_screen_service as service
from app.services.genome_wide_screen_registry import (
    RunInfo,
    latest_runs_by_group,
    list_runs,
    older_runs_by_group,
)
from app.services.llm_report_service import get_default_generator
from app.services.pichia_secretion_service import discover_project_paths
from app.ui.common import request_navigation

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

    _render_launch_controls(paths, active_runs)
    st.divider()
    _render_run_list(runs)
    st.divider()
    _render_results_section(paths, runs)


SCOPE_LABELS = {
    "gene": "全基因组（约1025个基因，小时级）",
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
                "run_name": run.run_name,
                "状态": run.status,
                "进度": run.progress_label,
                "靶点": ", ".join(run.targets),
                "模式": run.mode,
                "更新时间": run.updated_at,
                "心跳超时未更新": run.is_stale,
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
        scope_runs = [run for run in latest if run.scope == scope]
        if not scope_runs:
            continue
        st.markdown(f"**{label}**")
        scope_runs = sorted(scope_runs, key=lambda run: run.targets)
        st.dataframe(_run_table_frame(scope_runs), width='stretch', hide_index=True)

    unlabeled = [run for run in latest if run.scope not in SCOPE_LABELS]
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


def _render_results_section(paths, runs: list[RunInfo]) -> None:
    st.subheader("结果查看")
    done_runs = [run for run in runs if run.status == "done"]
    if not done_runs:
        st.caption("还没有已完成的运行。")
        return

    latest_done, older_done = _split_result_runs(runs)

    viewing_older: RunInfo | None = None
    if older_done:
        with st.expander(f"查看历史版本（{len(older_done)}个已被更新的运行取代，仍保留在磁盘上）"):
            viewing_older = st.selectbox(
                "选择一个历史版本查看（选中后会替换下面的结果）",
                options=[None, *older_done],
                format_func=lambda run: "（不查看，显示下方最新结果）" if run is None else f"{run.run_name}（{SCOPE_LABELS.get(run.scope, run.scope)}）",
                key="genome_wide_screen_older_run_select",
            )

    if viewing_older is not None:
        selected_run = viewing_older
        st.info(f"当前显示历史版本结果：{selected_run.run_name}（不是最新结果，上面的下拉框可以切回不查看）")
    else:
        if not latest_done:
            st.caption("当前没有最新的已完成运行；可在上方历史版本中查看已被新运行取代的完成结果。")
            return
        selected_run = st.selectbox(
            "选择要查看的运行", options=latest_done,
            format_func=lambda run: f"{run.run_name}（{SCOPE_LABELS.get(run.scope, run.scope)}）",
        )
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

    tabs = st.tabs([f"靶点：{target_id}" for target_id in target_ids] + (["靶点间差异"] if len(target_ids) >= 2 else []))
    for tab, target_id in zip(tabs, target_ids):
        with tab:
            _render_dimension_tables(per_target_results[target_id])
    if len(target_ids) >= 2:
        with tabs[-1]:
            divergence = analysis.analyze_target_divergence(frame, target_ids)
            st.markdown("同一个基因的KO，在不同靶点上效应差异最大的候选：")
            st.dataframe(divergence, width='stretch', hide_index=True)

    st.divider()
    _render_llm_report_section(run_name, per_target_results)


def _render_dimension_tables(result: analysis.DimensionalResults) -> None:
    st.metric("必需基因数（KO后任何测试生长速率下都不可行）", len(result.essential_genes))
    with st.expander(f"必需基因清单（{len(result.essential_genes)}）"):
        st.dataframe(result.essential_genes, width='stretch', hide_index=True)

    if not result.solver_inconclusive_rows.empty:
        st.warning(
            f"有 {len(result.solver_inconclusive_rows)} 个候选行因为求解超时或求解器失败而无法证明可行性；"
            "这些不是必需基因结论。"
        )
    with st.expander(f"求解未定的候选行（{len(result.solver_inconclusive_rows)}）", expanded=False):
        st.dataframe(result.solver_inconclusive_rows, width='stretch', hide_index=True)

    with st.expander(f"求解器重试证据（{len(result.solver_retry_evidence)}）", expanded=False):
        st.dataframe(result.solver_retry_evidence, width='stretch', hide_index=True)

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
        st.dataframe(df, width='stretch', hide_index=True)
        return
    event = st.dataframe(
        df,
        width='stretch',
        hide_index=True,
        on_select="rerun",
        selection_mode="single-row",
        key=table_key,
    )
    selected_rows = event.selection.rows if event and event.selection else []
    if selected_rows:
        row = df.iloc[selected_rows[0]]
        candidate_id = str(row["gene_id"])
        candidate_kind = str(row["candidate_kind"]) if "candidate_kind" in row else "gene"
        common_name = str(row["common_name"]) if "common_name" in row and row["common_name"] else ""
        display_label = f"{common_name} ({candidate_id})" if common_name else candidate_id
        if st.button(
            f"在仿真验证中核实 {display_label}（{intervention_type}，靶点 {target_id}）",
            key=f"{table_key}_verify_btn",
            type="primary",
        ):
            _apply_verify_prefill(target_id, candidate_id, intervention_type, candidate_kind)


def _verify_prefill_field_values(candidate_id: str, intervention_type: str, candidate_kind: str) -> dict[str, str]:
    """Pure routing decision behind _apply_verify_prefill - which draft field gets candidate_id.

    "gene" (full-genome screen rows) goes to the gene-ID inputs, which resolve through
    GPR; everything else ("catalog_reaction" from the curated catalog screen,
    "complex_oe_hypothesis" from the hypothetical whole-complex OE test, and any future
    non-gene candidate_kind) goes to the reaction-ID inputs, since those rows' gene_id
    field actually holds a direct reaction id with no gene to resolve. Allowlisting
    "gene" rather than denylisting "catalog_reaction" so a new non-gene candidate_kind
    fails safe (reaction routing) instead of silently landing in the wrong box.
    """
    is_gene = candidate_kind == "gene"
    return {
        "pichia_draft_ko_genes": candidate_id if (is_gene and intervention_type == "KO") else "",
        "pichia_draft_ko_reactions": candidate_id if (not is_gene and intervention_type == "KO") else "",
        "pichia_draft_oe_genes": candidate_id if (is_gene and intervention_type == "OE") else "",
        "pichia_draft_oe_reactions": candidate_id if (not is_gene and intervention_type == "OE") else "",
    }


def _apply_verify_prefill(target_id: str, candidate_id: str, intervention_type: str, candidate_kind: str = "gene") -> None:
    """Pre-fill 仿真验证 page's session_state and jump there.

    Replaces (not appends to) the KO/OE draft fields so the simulation isolates
    exactly this one candidate instead of mixing in whatever was left over from a
    previous exploration.
    """
    st.session_state["pichia_tab_selector"] = "仿真构建"
    st.session_state["pichia_draft_build_mode"] = "快速选择（内置模板）"
    st.session_state["pichia_template"] = target_id
    st.session_state.update(_verify_prefill_field_values(candidate_id, intervention_type, candidate_kind))
    request_navigation("仿真验证")
    st.rerun()


def _render_llm_report_section(run_name: str, per_target_results: dict[str, analysis.DimensionalResults]) -> None:
    st.subheader("生成总结报告")
    report_path = _report_cache_path(run_name)
    if report_path.exists():
        st.markdown(report_path.read_text(encoding="utf-8"))
        if st.button("重新生成报告"):
            report_path.unlink()
            st.rerun()
        return

    if st.button("用 LLM 生成总结和建议", type="primary"):
        with st.spinner("正在调用 LLM 总结……"):
            try:
                generator = get_default_generator()
                summaries = [result.to_summary_dict() for result in per_target_results.values()]
                report_text = generator.generate(
                    summaries, run_metadata={"run_name": run_name, "targets": list(per_target_results)}
                )
            except Exception as exc:  # noqa: BLE001 - surface the real error to the user, this is user-facing config
                st.error(f"报告生成失败：{exc}")
                return
        report_path.write_text(report_text, encoding="utf-8")
        st.rerun()


def _report_cache_path(run_name: str) -> Path:
    paths = _paths()
    directory = paths.local_runs_dir / run_name
    directory.mkdir(parents=True, exist_ok=True)
    return directory / "llm_report.md"


__all__ = ["render_genome_wide_screen"]
