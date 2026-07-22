from __future__ import annotations

from collections import Counter
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
TARGET_METADATA_KEY = "experiment_feedback_target_metadata"
BATCH_METADATA_KEY = "experiment_feedback_batch_metadata"
IMPORT_FORM_STATE_KEY = "experiment_feedback_import_form_state"
PREDICTION_PATH_KEY = "experiment_feedback_prediction_path"
INCREASE_THRESHOLD_KEY = "experiment_feedback_increase_threshold"
DECREASE_THRESHOLD_KEY = "experiment_feedback_decrease_threshold"
BASELINE_HIT_RATE_KEY = "experiment_feedback_baseline_hit_rate"
TOP_K_KEY = "experiment_feedback_top_k"


# --- 中文人类可读映射：内部枚举值/字段名 -> 展示文案 -------------------------

_INTERVENTION_TYPE_LABELS = {"control": "对照", "KO": "敲除", "OE": "过表达"}
_DIRECTION_LABELS = {"increase": "提升", "decrease": "降低", "neutral": "无明显变化", "": "未预测方向"}
_LINK_STATUS_LABELS = {
    "matched": "已匹配",
    "ambiguous": "有歧义",
    "missing_prediction": "无对应预测",
    "context_mismatch": "条件不匹配",
}
_ELIGIBILITY_LABELS = {"eligible": "可核对", "ineligible": "不可核对"}
_QUALITY_STATUS_LABELS = {"valid": "有效", "warning": "有警告", "invalid": "无效", "excluded": "已排除"}
_FERMENTATION_STATUS_LABELS = {
    "normal": "正常",
    "contamination": "污染",
    "culture_failed": "培养失败",
    "assay_failed": "检测失败",
    "other_excluded": "其他排除",
}
_MEASUREMENT_STATUS_LABELS = {
    "valid": "正常",
    "below_lod": "低于检测限",
    "below_loq": "低于定量限",
    "above_range": "超出标准曲线范围",
    "missing": "缺失",
    "assay_failed": "检测失败",
    "excluded": "已排除",
}
_RANKING_ASSESSMENT_LABELS = {
    "insufficient_evidence": "样本量不足，仅供参考",
    "descriptive_evidence_available": "样本量足够，可参考排名相关性",
}
_VALIDATION_CODE_LABELS = {
    "schema_validation_error": "数据格式不符合规则",
    "unit_validation_error": "检测单位不符合规则",
    "record_id_conflict": "同一记录出现冲突内容",
    "condition_missing": "发酵条件信息缺失",
    "missing_experiment_reference": "引用了不存在的实验记录",
    "missing_intervention_reference": "引用了不存在的改造记录",
    "import_warning": "导入提示",
}
_RECORD_TYPE_LABELS = {
    "experiment": "实验记录",
    "intervention": "改造记录",
    "measurement": "检测记录",
    "prediction_link": "预测匹配记录",
    "bundle": "整批数据",
}
_LINK_REASON_LABELS = {
    "gene_id_missing": "缺少基因编号",
    "common_name_only": "只有基因别名，无法唯一确认对应基因",
}
# These prefixes/codes mirror calibration.py's `_ineligible()` reason strings. The UI
# can't import pcsec_pichia directly (facade boundary, enforced by
# test_pichia_experiment_feedback_ui_contract.py), so this list is kept in sync by hand —
# if calibration.py's reason format changes, update this too; an unmatched code just falls
# back to displaying the raw string, it won't raise.
_INELIGIBILITY_PREFIX_LABELS = {
    "fermentation_data_status": "发酵状态异常",
    "experiment_quality_status": "数据质量",
    "prediction_link": "与模型预测的匹配",
}
_INELIGIBILITY_EXACT_LABELS = {
    "experiment_context_incomplete": "发酵条件信息不完整",
    "combination_intervention_not_attributable": "同一实验包含多个改造，无法归因到单个基因",
    "candidate_measurement_not_evaluable": "候选检测值不可用（缺失/失败/排除/超出标准曲线范围）",
    "candidate_measurement_context_ambiguous": "候选检测信号来源不唯一",
    "control_match_missing": "没有找到匹配的同批次对照",
    "control_value_zero": "对照检测值为零，无法计算比值",
    "prediction_link_missing": "缺少预测匹配记录",
}


def _humanize_reason_code(code: str) -> str:
    if code in _INELIGIBILITY_EXACT_LABELS:
        return _INELIGIBILITY_EXACT_LABELS[code]
    if ":" in code:
        prefix, _, suffix = code.partition(":")
        prefix_label = _INELIGIBILITY_PREFIX_LABELS.get(prefix)
        if prefix_label is not None:
            suffix_label = (
                _FERMENTATION_STATUS_LABELS.get(suffix)
                or _QUALITY_STATUS_LABELS.get(suffix)
                or _LINK_STATUS_LABELS.get(suffix)
                or _INELIGIBILITY_EXACT_LABELS.get(suffix)
                or suffix
            )
            return f"{prefix_label}：{suffix_label}"
    return code


def _humanize_list(values: Any, label_map: dict[str, str]) -> str:
    items = list(values) if isinstance(values, (list, tuple)) else ([values] if values else [])
    return "、".join(label_map.get(str(item), str(item)) for item in items) or "—"


def _humanize_reason_list(values: Any) -> str:
    items = list(values) if isinstance(values, (list, tuple)) else ([values] if values else [])
    return "；".join(_humanize_reason_code(str(item)) for item in items) or "—"


def _translate_rows(
    rows: list[dict[str, Any]],
    *,
    column_labels: dict[str, str],
    value_maps: dict[str, dict[str, str]] | None = None,
    list_columns: dict[str, dict[str, str]] | None = None,
    reason_columns: tuple[str, ...] = (),
) -> list[dict[str, Any]]:
    value_maps = value_maps or {}
    list_columns = list_columns or {}
    translated: list[dict[str, Any]] = []
    for row in rows:
        new_row: dict[str, Any] = {}
        for key, label in column_labels.items():
            if key not in row:
                continue
            value = row[key]
            if key in reason_columns:
                new_row[label] = _humanize_reason_list(value)
            elif key in list_columns:
                new_row[label] = _humanize_list(value, list_columns[key])
            elif key in value_maps:
                new_row[label] = value_maps[key].get(str(value), value)
            elif isinstance(value, (list, tuple)):
                new_row[label] = "、".join(str(item) for item in value) or "—"
            elif value is None:
                new_row[label] = "—"
            else:
                new_row[label] = value
        translated.append(new_row)
    return translated


def _restore_import_form_state() -> None:
    defaults: dict[str, Any] = {
        PREDICTION_PATH_KEY: "",
        TARGET_METADATA_KEY: "",
        BATCH_METADATA_KEY: "",
        RUN_NAME_KEY: f"feedback-{datetime.now().strftime('%Y%m%d-%H%M%S')}",
        INCREASE_THRESHOLD_KEY: 1.05,
        DECREASE_THRESHOLD_KEY: 0.95,
        BASELINE_HIT_RATE_KEY: 0.10,
        TOP_K_KEY: "5,10",
    }
    saved = dict(st.session_state.get(IMPORT_FORM_STATE_KEY) or {})
    for key, default in defaults.items():
        if key not in st.session_state:
            st.session_state[key] = saved.get(key, default)


def _sync_import_form_field(key: str) -> None:
    saved = dict(st.session_state.get(IMPORT_FORM_STATE_KEY) or {})
    saved[key] = st.session_state.get(key)
    st.session_state[IMPORT_FORM_STATE_KEY] = saved


def _normalize_import_form_options(prediction_paths: list[str]) -> None:
    valid_options = {
        PREDICTION_PATH_KEY: {"", *prediction_paths},
        TARGET_METADATA_KEY: {"", "hLF", "OPN"},
    }
    for key, options in valid_options.items():
        if st.session_state.get(key) not in options:
            st.session_state[key] = ""
            _sync_import_form_field(key)


def render_experiment_feedback() -> None:
    st.header("实验反馈闭环")
    st.caption(
        "导入实验记录不代表已经完成核对；只有通过数据校验、找到匹配对照、且和模型预测唯一对应的记录，"
        "才会进入历史数据核对的统计。"
    )
    st.info("所有上传文件、冲突记录、匹配结果和核对结果只保存在本地文件（不进入版本库）；原始实验记录不会发送给大模型。")

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

    if selected_payload:
        _render_summary(selected_payload)

    import_tab, validation_tab, linkage_tab, calibration_tab = st.tabs(
        ["导入 / 修正", "数据校验", "预测匹配", "历史数据核对"]
    )
    with import_tab:
        _render_import()
    with validation_tab:
        _render_validation(selected_payload)
    with linkage_tab:
        _render_linkage(selected_payload)
    with calibration_tab:
        _render_calibration(selected_payload)


def _render_summary(payload: dict[str, Any]) -> None:
    validation = payload.get("validation") or {}
    linkage = payload.get("linkage") or {}
    calibration = payload.get("calibration") or {}
    errors = list(validation.get("errors") or [])

    st.subheader("这批数据说明了什么")
    if errors:
        st.error(
            f"这批数据里有 {len(errors)} 处问题会阻止后续核对，请到"
            "「数据校验」标签页查看具体问题，修正后重新上传为新的 run。"
        )
        return

    records = list(calibration.get("records") or [])
    control_count = int(linkage.get("control_count", 0))
    if not records:
        st.info(f"这批数据包含 {control_count} 条对照记录，没有候选改造记录可核对。")
        return

    eligible = [record for record in records if record.get("eligibility_status") == "eligible"]
    hits = [record for record in eligible if record.get("hit") is True]
    reason_counts: Counter[str] = Counter()
    for record in records:
        if record.get("eligibility_status") != "eligible":
            for reason in record.get("ineligibility_reasons") or []:
                reason_counts[_humanize_reason_code(str(reason))] += 1

    lines = [
        f"这批数据包含 {control_count} 条对照、{len(records)} 条候选改造记录。",
        f"其中 {len(eligible)} 条候选完成了历史数据核对，{len(records) - len(eligible)} 条因故未能核对。",
    ]
    if eligible:
        lines.append(f"完成核对的候选里，有 {len(hits)} 条达到了命中阈值（候选/对照比值超过设定门槛）。")
    if reason_counts:
        top_reasons = "，".join(f"{reason}（{count} 条）" for reason, count in reason_counts.most_common(3))
        lines.append(f"未能核对的主要原因：{top_reasons}。")
    for target in calibration.get("targets") or []:
        if target.get("ranking_assessment") == "insufficient_evidence":
            lines.append(
                f"{target.get('target_id')} 的可比较排名对数量为 "
                f"{target.get('comparable_rank_pair_count', 0)}，样本量不足，排名相关性仅供参考，不构成结论。"
            )
    st.info("\n\n".join(lines))
    st.caption("以上摘要由现有校验/匹配/核对结果自动组装，不改变任何判定逻辑，也不会写回模型或推荐等级。")


def _render_import() -> None:
    _restore_import_form_state()
    st.subheader("导入脱敏实验 bundle")
    st.caption("支持内部通用格式与研发发酵宽表。XLSX 可使用 metadata 工作表；目标蛋白和批次也可由下方表单补齐。修正问题后重新上传为新的 run，原 run 不覆盖。")
    experiment_upload = st.file_uploader(
        "实验文件（CSV / XLSX / JSONL）",
        type=["csv", "xlsx", "jsonl"],
        key=EXPERIMENT_UPLOAD_KEY,
    )
    prediction_paths = list_prediction_fact_packs("local_runs")
    _normalize_import_form_options(prediction_paths)
    prediction_path = st.selectbox(
        "已有的模型预测结果（可选）",
        options=[""] + prediction_paths,
        format_func=lambda value: Path(value).parent.name if value else "不选择已有预测结果",
        key=PREDICTION_PATH_KEY,
        on_change=_sync_import_form_field,
        args=(PREDICTION_PATH_KEY,),
    )
    prediction_upload = st.file_uploader(
        "或上传模型预测结果 fact_pack.json（可选，优先于已有路径）",
        type=["json"],
        key=PREDICTION_UPLOAD_KEY,
    )
    metadata_target, metadata_batch = st.columns(2)
    metadata_target.selectbox(
        "模板目标蛋白（文件未填写时补齐）",
        options=["", "hLF", "OPN"],
        format_func=lambda value: value or "从文件读取",
        key=TARGET_METADATA_KEY,
        on_change=_sync_import_form_field,
        args=(TARGET_METADATA_KEY,),
    )
    metadata_batch.text_input(
        "模板批次（文件未填写时补齐）",
        key=BATCH_METADATA_KEY,
        placeholder="例如 B01",
        on_change=_sync_import_form_field,
        args=(BATCH_METADATA_KEY,),
    )
    st.text_input(
        "Run 名称",
        key=RUN_NAME_KEY,
        on_change=_sync_import_form_field,
        args=(RUN_NAME_KEY,),
    )
    with st.expander("历史核对参数设置", expanded=False):
        increase_threshold = st.number_input(
            "命中阈值（候选/对照检测值比值，超过算命中）",
            min_value=1.0001,
            step=0.01,
            key=INCREASE_THRESHOLD_KEY,
            on_change=_sync_import_form_field,
            args=(INCREASE_THRESHOLD_KEY,),
        )
        decrease_threshold = st.number_input(
            "降低方向阈值（候选/对照比值低于此值算明显降低）",
            min_value=0.0001,
            max_value=0.9999,
            step=0.01,
            key=DECREASE_THRESHOLD_KEY,
            on_change=_sync_import_form_field,
            args=(DECREASE_THRESHOLD_KEY,),
        )
        baseline_hit_rate = st.number_input(
            "基线命中率（用于计算相对富集倍数）",
            min_value=0.0001,
            max_value=1.0,
            step=0.01,
            key=BASELINE_HIT_RATE_KEY,
            on_change=_sync_import_form_field,
            args=(BASELINE_HIT_RATE_KEY,),
        )
        top_k_text = st.text_input(
            "Top-K 命中率统计范围（逗号分隔，如 5,10）",
            key=TOP_K_KEY,
            on_change=_sync_import_form_field,
            args=(TOP_K_KEY,),
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
                experiment_metadata={
                    key: value
                    for key, value in {
                        "target_id": st.session_state.get(TARGET_METADATA_KEY),
                        "batch_id": st.session_state.get(BATCH_METADATA_KEY),
                    }.items()
                    if str(value or "").strip()
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
    st.subheader("数据校验")
    if not payload:
        st.info("请选择或导入一个 run。")
        return
    validation = payload.get("validation") or {}
    errors = list(validation.get("errors") or [])
    warnings = list(validation.get("warnings") or [])
    col_valid, col_errors, col_warnings = st.columns(3)
    col_valid.metric("数据格式校验通过", "是" if validation.get("is_valid") else "否")
    col_errors.metric("错误数", len(errors))
    col_warnings.metric("警告数", len(warnings))
    issue_columns = {
        "code": "问题类型",
        "message": "详细说明",
        "record_type": "记录类型",
        "record_id": "记录编号",
    }
    issue_value_maps = {"code": _VALIDATION_CODE_LABELS, "record_type": _RECORD_TYPE_LABELS}
    if errors:
        st.error("存在阻止正式核对的问题；修正后请作为新 run 重新上传。")
        st.dataframe(
            _translate_rows(errors, column_labels=issue_columns, value_maps=issue_value_maps),
            width='stretch',
            hide_index=True,
        )
    if warnings:
        st.warning("以下记录已保留，但需要复核。")
        st.dataframe(
            _translate_rows(warnings, column_labels=issue_columns, value_maps=issue_value_maps),
            width='stretch',
            hide_index=True,
        )
    run_dir = Path(str(payload.get("run_dir") or ""))
    conflict_bytes = export_experiment_feedback_issues(run_dir, issue_kind="conflicts")
    warning_bytes = export_experiment_feedback_issues(run_dir, issue_kind="warnings")
    export_conflicts, export_warnings = st.columns(2)
    export_conflicts.download_button(
        "导出冲突记录（conflicts.jsonl）",
        data=conflict_bytes,
        file_name=f"{payload.get('run_name', 'run')}_conflicts.jsonl",
        mime="application/x-ndjson",
        disabled=not conflict_bytes,
    )
    export_warnings.download_button(
        "导出警告记录（warnings.jsonl）",
        data=warning_bytes,
        file_name=f"{payload.get('run_name', 'run')}_warnings.jsonl",
        mime="application/x-ndjson",
        disabled=not warning_bytes,
    )


def _render_linkage(payload: dict[str, Any] | None) -> None:
    st.subheader("预测与实验的匹配情况")
    st.caption("匹配状态包括已匹配、有歧义、无对应预测、条件不匹配；只有基因别名、没有基因编号时不会单独形成已匹配。")
    if not payload:
        st.info("请选择或导入一个 run。")
        return
    linkage = payload.get("linkage") or {}
    columns = st.columns(4)
    columns[0].metric("已匹配", int(linkage.get("matched_count", 0)))
    columns[1].metric("有歧义", int(linkage.get("ambiguous_count", 0)))
    columns[2].metric("无对应预测", int(linkage.get("missing_prediction_count", 0)))
    columns[3].metric("条件不匹配", int(linkage.get("context_mismatch_count", 0)))
    links = list(linkage.get("links") or [])
    if links:
        column_labels = {
            "experiment_id": "实验编号",
            "target_id": "目标蛋白",
            "gene_id": "基因",
            "common_name": "基因别名",
            "intervention_type": "改造类型",
            "status": "匹配状态",
            "reason": "说明",
            "prediction_rank": "预测排名",
            "predicted_direction": "预测方向",
            "evidence_tier": "证据等级",
            "recommendation_tier": "推荐等级",
            "prediction_run_id": "预测批次",
            "evidence_id": "证据编号",
        }
        value_maps = {
            "intervention_type": _INTERVENTION_TYPE_LABELS,
            "status": _LINK_STATUS_LABELS,
            "predicted_direction": _DIRECTION_LABELS,
            "reason": _LINK_REASON_LABELS,
        }
        st.dataframe(
            _translate_rows(links, column_labels=column_labels, value_maps=value_maps),
            width='stretch',
            hide_index=True,
        )
    else:
        st.info("当前 run 没有非对照的匹配记录。")


def _render_calibration(payload: dict[str, Any] | None) -> None:
    st.subheader("hLF / OPN 历史数据核对")
    st.caption("不可核对的记录不会进入指标分母；阴性结果、失败、缺对照和有歧义匹配仍完整保留在记录列表中，不会被隐藏。")
    if not payload:
        st.info("请选择或导入一个 run。")
        return
    calibration = payload.get("calibration") or {}
    if not calibration.get("available"):
        st.warning(f"当前 run 不可核对：{calibration.get('reason') or '未知原因'}")
        return
    config = calibration.get("config") or {}
    st.markdown(
        f"**核对参数**：命中阈值 {config.get('increase_threshold_ratio', '—')} 、"
        f"降低方向阈值 {config.get('decrease_threshold_ratio', '—')} 、"
        f"基线命中率 {config.get('baseline_hit_rate', '—')} 、"
        f"Top-K {list(config.get('top_k') or [])} 、"
        f"主要检测类型/区室：{config.get('primary_assay_type', '—')}/{config.get('primary_compartment', '—')}"
    )
    report_bytes = export_experiment_feedback_report(str(payload.get("run_dir") or ""))
    if report_bytes:
        st.download_button(
            "下载预测 vs 实验核对报告",
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
            st.info(f"{target_id} 尚无核对记录。")
            continue
        comparable_rank_pairs = int(target.get("comparable_rank_pair_count", 0))
        minimum_rank_pairs = int((calibration.get("config") or {}).get("minimum_rank_pairs", 2))
        ranking_assessment = str(target.get("ranking_assessment") or "insufficient_evidence")
        metric_cols = st.columns(4)
        metric_cols[0].metric("可核对", int(target.get("eligible_count", 0)))
        metric_cols[1].metric("不可核对", int(target.get("ineligible_count", 0)))
        consistency = target.get("direction_consistency_rate")
        metric_cols[2].metric(
            "方向一致率",
            "无数据" if consistency is None else f"{float(consistency):.1%}",
        )
        metric_cols[3].metric("排名可比较对", f"{comparable_rank_pairs}/{minimum_rank_pairs}")
        assessment_label = _RANKING_ASSESSMENT_LABELS.get(ranking_assessment, ranking_assessment)
        if ranking_assessment == "insufficient_evidence":
            st.warning(
                f"排序证据状态：{assessment_label}——仅有 {comparable_rank_pairs} 个"
                f"预测排名+观测比值可比较对，低于最低要求 {minimum_rank_pairs}。"
                "下方 Top-K、命中率与富集倍数仅作描述性展示，不构成排名可靠性结论。"
            )
        else:
            st.info(f"排序证据状态：{assessment_label}；可比较排名对 {comparable_rank_pairs}/{minimum_rank_pairs}。")
        st.caption("Top-K 与富集倍数不会自动修改推荐等级或模型约束。")
        top_k_rows = target.get("top_k_metrics") or []
        if top_k_rows:
            st.dataframe(
                _translate_rows(
                    top_k_rows,
                    column_labels={
                        "k": "Top-K",
                        "tested_count": "参与统计的候选数",
                        "hit_count": "命中数",
                        "hit_rate": "命中率",
                        "relative_baseline_enrichment": "相对基线富集倍数",
                    },
                ),
                width='stretch',
                hide_index=True,
            )
        tier_rows = target.get("evidence_tier_metrics") or []
        if tier_rows:
            st.dataframe(
                _translate_rows(
                    tier_rows,
                    column_labels={
                        "evidence_tier": "证据等级",
                        "tested_count": "参与统计的候选数",
                        "hit_count": "命中数",
                        "hit_rate": "命中率",
                    },
                ),
                width='stretch',
                hide_index=True,
            )
        target_records = [record for record in records if record.get("target_id") == target_id]
        if target_records:
            column_labels = {
                "gene_id": "基因",
                "intervention_type": "改造类型",
                "evidence_id": "证据编号",
                "prediction_rank": "预测排名",
                "evidence_tier": "证据等级",
                "recommendation_tier": "推荐等级",
                "predicted_direction": "预测方向",
                "fermentation_data_status": "发酵状态",
                "eligibility_status": "是否可核对",
                "ineligibility_reasons": "不可核对原因",
                "candidate_value": "候选检测值",
                "control_value": "对照检测值",
                "observed_ratio": "观测比值",
                "observed_direction": "观测方向",
                "direction_consistent": "方向是否一致",
                "hit": "是否命中",
                "measurement_statuses": "检测状态",
                "experiment_id": "实验编号",
            }
            value_maps = {
                "intervention_type": _INTERVENTION_TYPE_LABELS,
                "predicted_direction": _DIRECTION_LABELS,
                "observed_direction": _DIRECTION_LABELS,
                "fermentation_data_status": _FERMENTATION_STATUS_LABELS,
                "eligibility_status": _ELIGIBILITY_LABELS,
            }
            list_columns = {"measurement_statuses": _MEASUREMENT_STATUS_LABELS}
            st.dataframe(
                _translate_rows(
                    target_records,
                    column_labels=column_labels,
                    value_maps=value_maps,
                    list_columns=list_columns,
                    reason_columns=("ineligibility_reasons",),
                ),
                width='stretch',
                hide_index=True,
            )


__all__ = ["render_experiment_feedback"]
