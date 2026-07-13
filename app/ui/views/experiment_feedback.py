from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Any

import streamlit as st

from app.services.pichia_experiment_feedback_service import (
    DEFAULT_EXPERIMENT_FEEDBACK_ROOT,
    export_experiment_feedback_issues,
    export_experiment_feedback_report,
    list_experiment_feedback_runs,
    list_prediction_fact_packs,
    load_experiment_feedback_run,
    submit_experiment_feedback_import,
)


SELECTED_RUN_KEY = "experiment_feedback_selected_run"
PENDING_SELECTED_RUN_KEY = "experiment_feedback_selected_run_pending"
LAST_IMPORT_KEY = "experiment_feedback_last_import"
RUN_NAME_KEY = "experiment_feedback_run_name"
EXPERIMENT_UPLOAD_KEY = "experiment_feedback_experiment_upload"
PREDICTION_UPLOAD_KEY = "experiment_feedback_prediction_upload"


def render_experiment_feedback() -> None:
    st.header("实验反馈闭环")
    st.caption("导入实验记录并不等于模型已校准；只有通过 validation、control matching 和唯一 prediction linkage 的记录才进入统计。")
    st.info("所有上传、冲突、linkage 与 calibration 产物只写入 ignored local_runs；原始实验记录不会进入 LLM。")

    runs = list_experiment_feedback_runs(DEFAULT_EXPERIMENT_FEEDBACK_ROOT)
    run_names = [str(row.get("run_name") or "") for row in runs if row.get("run_name")]
    if PENDING_SELECTED_RUN_KEY in st.session_state:
        st.session_state[SELECTED_RUN_KEY] = st.session_state.pop(PENDING_SELECTED_RUN_KEY)
    if st.session_state.get(SELECTED_RUN_KEY) not in run_names:
        st.session_state[SELECTED_RUN_KEY] = run_names[0] if run_names else ""
    selected_run = st.selectbox(
        "当前实验反馈 run",
        options=[""] + run_names,
        format_func=lambda value: value or "尚未选择 run",
        key=SELECTED_RUN_KEY,
    )
    selected_payload = (
        load_experiment_feedback_run(DEFAULT_EXPERIMENT_FEEDBACK_ROOT / selected_run)
        if selected_run
        else None
    )

    import_tab, validation_tab, linkage_tab, calibration_tab = st.tabs(
        ["导入 / 修正", "Validation / Conflicts", "Linkage", "Calibration"]
    )
    with import_tab:
        _render_import()
    with validation_tab:
        _render_validation(selected_payload)
    with linkage_tab:
        _render_linkage(selected_payload)
    with calibration_tab:
        _render_calibration(selected_payload)


def _render_import() -> None:
    st.subheader("导入脱敏实验 bundle")
    st.caption("支持 canonical JSONL，或 record_type + payload_json 的 CSV 导入适配格式。修正问题后重新上传为新的 run，原 run 不覆盖。")
    experiment_upload = st.file_uploader(
        "实验文件（CSV / JSONL）",
        type=["csv", "jsonl"],
        key=EXPERIMENT_UPLOAD_KEY,
    )
    prediction_paths = list_prediction_fact_packs("local_runs")
    prediction_path = st.selectbox(
        "已有 prediction fact pack（可选）",
        options=[""] + prediction_paths,
        format_func=lambda value: Path(value).parent.name if value else "不选择已有 fact pack",
        key="experiment_feedback_prediction_path",
    )
    prediction_upload = st.file_uploader(
        "或上传 prediction fact_pack.json（可选，优先于已有路径）",
        type=["json"],
        key=PREDICTION_UPLOAD_KEY,
    )
    if RUN_NAME_KEY not in st.session_state:
        st.session_state[RUN_NAME_KEY] = f"feedback-{datetime.now().strftime('%Y%m%d-%H%M%S')}"
    st.text_input("Run 名称", key=RUN_NAME_KEY)
    with st.expander("Calibration 配置", expanded=False):
        increase_threshold = st.number_input(
            "命中阈值（candidate/control ratio）",
            min_value=1.0001,
            value=1.05,
            step=0.01,
            key="experiment_feedback_increase_threshold",
        )
        decrease_threshold = st.number_input(
            "降低方向阈值",
            min_value=0.0001,
            max_value=0.9999,
            value=0.95,
            step=0.01,
            key="experiment_feedback_decrease_threshold",
        )
        baseline_hit_rate = st.number_input(
            "基线命中率",
            min_value=0.0001,
            max_value=1.0,
            value=0.10,
            step=0.01,
            key="experiment_feedback_baseline_hit_rate",
        )
        top_k_text = st.text_input(
            "Top-K（逗号分隔）",
            value="5,10",
            key="experiment_feedback_top_k",
        )
    if st.button(
        "导入并校验",
        type="primary",
        disabled=experiment_upload is None,
        key="experiment_feedback_submit_import",
    ):
        try:
            top_k = tuple(int(item.strip()) for item in top_k_text.split(",") if item.strip())
            result = submit_experiment_feedback_import(
                experiment_filename=experiment_upload.name,
                experiment_bytes=experiment_upload.getvalue(),
                prediction_filename=prediction_upload.name if prediction_upload else "",
                prediction_bytes=prediction_upload.getvalue() if prediction_upload else None,
                prediction_path=None if prediction_upload else prediction_path or None,
                run_name=str(st.session_state.get(RUN_NAME_KEY) or ""),
                output_root=DEFAULT_EXPERIMENT_FEEDBACK_ROOT,
                calibration_config={
                    "increase_threshold_ratio": float(increase_threshold),
                    "decrease_threshold_ratio": float(decrease_threshold),
                    "baseline_hit_rate": float(baseline_hit_rate),
                    "top_k": top_k,
                },
            )
        except Exception as exc:
            st.error(f"导入失败：{type(exc).__name__}: {exc}")
        else:
            st.session_state[LAST_IMPORT_KEY] = result
            st.session_state[PENDING_SELECTED_RUN_KEY] = result["run_name"]
            st.success(f"已创建 run：{result['run_name']}")
            st.rerun()


def _render_validation(payload: dict[str, Any] | None) -> None:
    st.subheader("Validation / Conflicts")
    if not payload:
        st.info("请选择或导入一个 run。")
        return
    validation = payload.get("validation") or {}
    errors = list(validation.get("errors") or [])
    warnings = list(validation.get("warnings") or [])
    col_valid, col_errors, col_warnings = st.columns(3)
    col_valid.metric("Schema 有效", "是" if validation.get("is_valid") else "否")
    col_errors.metric("Errors", len(errors))
    col_warnings.metric("Warnings", len(warnings))
    if errors:
        st.error("存在阻止正式校准的问题；修正后请作为新 run 重新上传。")
        st.dataframe(errors, use_container_width=True, hide_index=True)
    if warnings:
        st.warning("以下记录已保留，但需要复核。")
        st.dataframe(warnings, use_container_width=True, hide_index=True)
    run_dir = Path(str(payload.get("run_dir") or ""))
    conflict_bytes = export_experiment_feedback_issues(run_dir, issue_kind="conflicts")
    warning_bytes = export_experiment_feedback_issues(run_dir, issue_kind="warnings")
    export_conflicts, export_warnings = st.columns(2)
    export_conflicts.download_button(
        "导出 conflicts.jsonl",
        data=conflict_bytes,
        file_name=f"{payload.get('run_name', 'run')}_conflicts.jsonl",
        mime="application/x-ndjson",
        disabled=not conflict_bytes,
    )
    export_warnings.download_button(
        "导出 warnings.jsonl",
        data=warning_bytes,
        file_name=f"{payload.get('run_name', 'run')}_warnings.jsonl",
        mime="application/x-ndjson",
        disabled=not warning_bytes,
    )


def _render_linkage(payload: dict[str, Any] | None) -> None:
    st.subheader("Prediction-to-experiment Linkage")
    st.caption("状态包括 matched、ambiguous、missing_prediction、context_mismatch；common name 不会单独形成 matched。")
    if not payload:
        st.info("请选择或导入一个 run。")
        return
    linkage = payload.get("linkage") or {}
    columns = st.columns(4)
    columns[0].metric("Matched", int(linkage.get("matched_count", 0)))
    columns[1].metric("Ambiguous", int(linkage.get("ambiguous_count", 0)))
    columns[2].metric("Missing", int(linkage.get("missing_prediction_count", 0)))
    columns[3].metric("Context mismatch", int(linkage.get("context_mismatch_count", 0)))
    links = list(linkage.get("links") or [])
    if links:
        st.dataframe(links, use_container_width=True, hide_index=True)
    else:
        st.info("当前 run 没有非 control linkage 记录。")


def _render_calibration(payload: dict[str, Any] | None) -> None:
    st.subheader("hLF / OPN Calibration")
    st.caption("不可校准记录不会进入指标分母；阴性、失败、缺 control 和 ambiguous linkage 仍保留在 records 中。")
    if not payload:
        st.info("请选择或导入一个 run。")
        return
    calibration = payload.get("calibration") or {}
    if not calibration.get("available"):
        st.warning(f"当前 run 不可校准：{calibration.get('reason') or 'unknown'}")
        return
    st.json(calibration.get("config") or {}, expanded=False)
    report_bytes = export_experiment_feedback_report(str(payload.get("run_dir") or ""))
    if report_bytes:
        st.download_button(
            "下载 prediction-vs-experiment 报告",
            data=report_bytes,
            file_name=f"{payload.get('run_name', 'run')}_prediction_experiment_report.md",
            mime="text/markdown",
        )
    records = list(calibration.get("records") or [])
    target_rows = {str(item.get("target_id")): item for item in calibration.get("targets") or []}
    for target_id in ("hLF", "OPN"):
        target = target_rows.get(target_id)
        st.markdown(f"### {target_id}")
        if not target:
            st.info(f"{target_id} 尚无 calibration 记录。")
            continue
        metric_cols = st.columns(3)
        metric_cols[0].metric("Eligible", int(target.get("eligible_count", 0)))
        metric_cols[1].metric("不可校准", int(target.get("ineligible_count", 0)))
        consistency = target.get("direction_consistency_rate")
        metric_cols[2].metric(
            "方向一致率",
            "N/A" if consistency is None else f"{float(consistency):.1%}",
        )
        st.dataframe(target.get("top_k_metrics") or [], use_container_width=True, hide_index=True)
        st.dataframe(target.get("evidence_tier_metrics") or [], use_container_width=True, hide_index=True)
        target_records = [record for record in records if record.get("target_id") == target_id]
        st.dataframe(target_records, use_container_width=True, hide_index=True)


__all__ = ["render_experiment_feedback"]
