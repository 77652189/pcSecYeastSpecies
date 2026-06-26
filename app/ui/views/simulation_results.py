from __future__ import annotations

import json
import time
from pathlib import Path

import pandas as pd
import streamlit as st

from app.services.pichia_background_tasks import (
    load_last_result,
    poll_background_simulation,
    save_last_result,
)
from app.ui.common import PATHS
from app.ui.views.candidate_path_graph import render_secretion_path_graph


def display_value(value: object, fallback: str = "") -> str:
    if value is None:
        return fallback
    try:
        if pd.isna(value):
            return fallback
    except TypeError:
        pass
    text = str(value)
    if text.lower() in {"nan", "none", "nat"}:
        return fallback
    return text


def target_semantics_label(value: object) -> str:
    labels = {
        "project_defined_hLF": "项目定义 hLF（用户提供序列）",
        "product_target": "产品草案目标",
        "opn_leader_candidate": "OPN leader 候选",
        "custom": "自定义目标",
        "native_signal_plus_mature_hLF": "人源天然信号肽 + mature hLF",
        "mature_secreted_with_leader_candidate": "成熟蛋白 + leader 候选",
        "mature_secreted": "成熟分泌蛋白",
        "custom_user_sequence": "用户自定义序列",
        "user_provided_as_provided": "用户提供序列，按原样使用",
        "as_provided": "按提供序列记录",
    }
    text = display_value(value)
    return labels.get(text, text or "未声明")


def render_pichia_results() -> None:
    tsp = st.session_state.get("pichia_draft_task_status_path")
    if tsp:
        status, msg, result = poll_background_simulation(tsp)
        if status == "done" and result:
            st.session_state["last_pichia_secretion_draft_response"] = result
            save_last_result(result, PATHS)
            st.session_state.pop("pichia_draft_task_status_path", None)
            st.session_state.pop("pichia_draft_task_id", None)
            st.rerun()
        elif status == "error":
            st.error(f"仿真失败：{msg}")
            st.session_state.pop("pichia_draft_task_status_path", None)
            st.session_state.pop("pichia_draft_task_id", None)
        elif status in ("pending", "running"):
            st.info(f"⏳ {msg}")
            time.sleep(2)
            st.rerun()
        elif status == "lost":
            st.error(f"状态丢失：{msg}")
            st.session_state.pop("pichia_draft_task_status_path", None)
            st.session_state.pop("pichia_draft_task_id", None)

    if "last_pichia_secretion_draft_response" not in st.session_state:
        cached = load_last_result(PATHS)
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
        st.dataframe(
            pd.DataFrame(
                [
                    {"参数": key, "值": value}
                    for key, value in {
                        "目标": data.get("target_id"),
                        "状态": data.get("result_status"),
                        "MATLAB": data.get("matlab_alignment_status"),
                        "目标值": data.get("objective_value"),
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


def _screen_status_label(raw_status: object, success: object) -> str:
    success_text = str(success).strip().lower()
    if success is True or success_text == "true":
        return "求解成功"
    status = display_value(raw_status)
    return {
        "2": "约束不可行",
        "3": "目标无界",
        "4": "求解器数值错误",
        "missing_reaction": "反应未找到",
        "unresolved_gene": "基因未解析",
        "unresolved_reaction": "反应未解析",
        "missing_objective": "目标反应未找到",
    }.get(status, "求解失败")


def _normalise_candidate_frame_for_display(frame: pd.DataFrame) -> pd.DataFrame:
    frame = frame.copy()
    if "solver_status_label" not in frame.columns:
        frame["solver_status_label"] = ""
    if "failure_reason" not in frame.columns:
        frame["failure_reason"] = ""
    for idx, row in frame.iterrows():
        label = _screen_status_label(row.get("status"), row.get("success"))
        if not display_value(row.get("solver_status_label")):
            frame.at[idx, "solver_status_label"] = label
        if label != "求解成功" and not display_value(row.get("failure_reason")):
            frame.at[idx, "failure_reason"] = label
        if display_value(row.get("status")) == "2" and display_value(row.get("effect_label")) in {"", "求解失败"}:
            frame.at[idx, "effect_label"] = "约束不可行"
    return frame


def _candidate_effect_counts(frame: pd.DataFrame) -> dict[str, int]:
    counts = {"提升分泌": 0, "降低分泌": 0, "约束不可行": 0, "求解失败": 0, "未解析": 0, "无明显变化": 0}
    for _, row in frame.iterrows():
        effect = display_value(row.get("effect_label"), "未解析")
        status = display_value(row.get("status"))
        solver_status = display_value(row.get("solver_status_label") or row.get("failure_reason"))
        if status == "2" or "不可行" in effect or "不可行" in solver_status:
            counts["约束不可行"] += 1
        elif "提升" in effect:
            counts["提升分泌"] += 1
        elif "降低" in effect:
            counts["降低分泌"] += 1
        elif "无明显" in effect:
            counts["无明显变化"] += 1
        elif "失败" in effect or "失败" in solver_status:
            counts["求解失败"] += 1
        else:
            counts["未解析"] += 1
    return {key: value for key, value in counts.items() if value}


def _render_candidate_outputs(candidate_path: str, summary_path: str | None) -> None:
    try:
        frame = pd.read_csv(candidate_path)
    except Exception as exc:
        st.warning(f"候选表读取失败：{exc}")
        return
    if frame.empty:
        st.info("候选表为空。")
        return
    frame = _normalise_candidate_frame_for_display(frame)
    if "delta_objective" in frame.columns:
        frame = frame.sort_values("delta_objective", ascending=False, na_position="last")
    effect_counts = _candidate_effect_counts(frame)
    if effect_counts:
        st.markdown("**候选分类汇总**")
        cols = st.columns(len(effect_counts))
        for col, (label, value) in zip(cols, effect_counts.items()):
            col.metric(label, value)
        if effect_counts.get("约束不可行"):
            st.info("「约束不可行」表示当前固定生长速率和约束组合下没有可行解，不等同于真实发酵条件必然不可行。")

    labels = {
        "input_gene_id": "输入基因",
        "resolved_reaction_id": "解析反应",
        "reaction_id": "反应",
        "intervention_type": "扰动类型",
        "effect_label": "效果",
        "objective_value": "目标值",
        "baseline_objective_value": "基线",
        "candidate_id": "候选ID",
        "delta_objective": "Δ目标值",
        "secretory_process": "分泌环节",
        "solver_status_label": "求解状态",
        "failure_reason": "失败原因",
    }
    display = frame.rename(columns={key: value for key, value in labels.items() if key in frame.columns})
    st.dataframe(display, use_container_width=True, hide_index=True)

    if len(frame) > 0:
        selectable = frame.reset_index(drop=True)
        row_labels: list[str] = []
        for idx, row in selectable.iterrows():
            gene = display_value(row.get("gene_id") or row.get("input_gene_id") or row.get("candidate_id"))
            effect = display_value(row.get("effect_label"), "未解析")
            delta = display_value(row.get("delta_objective"), "无可行目标值")
            row_labels.append(f"{idx + 1}. {gene} | {effect} | Δ={delta}")
        selected_index = st.selectbox("选择一行查看分泌路径图", range(len(row_labels)), format_func=lambda i: row_labels[i])
        row = selectable.iloc[int(selected_index)]
        summary = None
        if summary_path and Path(summary_path).exists():
            try:
                summary = json.loads(Path(summary_path).read_text(encoding="utf-8"))
            except Exception:
                summary = {}
        render_secretion_path_graph(row.to_dict(), summary or {})


__all__ = ["render_pichia_results", "target_semantics_label"]
