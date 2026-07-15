from __future__ import annotations

from datetime import datetime
from typing import Any

import streamlit as st

from app.services.pichia_oe_capacity_service import (
    DEFAULT_OE_CAPACITY_CANDIDATE_ROOT,
    export_oe_capacity_report,
    list_oe_capacity_runs,
    list_oe_capacity_targets,
    load_oe_capacity_candidate_review,
    preview_oe_capacity_candidate,
    preview_oe_capacity_promotion,
    promote_oe_capacity_candidate_selection,
    submit_oe_capacity_screen,
)


TARGET_KEY = "oe_capacity_target_id"
ACTIVE_TARGET_KEY = "oe_capacity_active_form_target"
FORM_STATE_KEY = "oe_capacity_form_state_by_target"
LAST_PREVIEWS_KEY = "oe_capacity_last_previews_by_target"
LAST_RUNS_KEY = "oe_capacity_last_runs_by_target"
CANDIDATE_SELECTIONS_KEY = "oe_capacity_candidate_selection_by_target"
CANDIDATE_PREVIEWS_KEY = "oe_capacity_promotion_preview_by_target"
CANDIDATE_WIDGET_KEY = "oe_capacity_widget_candidate_ids"
CANDIDATE_ACTIVE_TARGET_KEY = "oe_capacity_candidate_active_target"

FORM_WIDGET_KEYS = {
    "gene_id": "oe_capacity_widget_gene_id",
    "dose_mode": "oe_capacity_widget_dose_mode",
    "product_mode": "oe_capacity_widget_product_mode",
    "multiplier": "oe_capacity_widget_multiplier",
    "scenarios": "oe_capacity_widget_scenarios",
    "compare_proxy": "oe_capacity_widget_compare_proxy",
    "run_name": "oe_capacity_widget_run_name",
}


def _default_form_state(target_id: str) -> dict[str, Any]:
    return {
        "gene_id": "PAS_chr2-1_0047",
        "dose_mode": "明确倍数",
        "product_mode": "相对未校准",
        "multiplier": 2.0,
        "scenarios": ["low", "nominal", "high"],
        "compare_proxy": True,
        "run_name": f"oe-{target_id}-{datetime.now().strftime('%Y%m%d-%H%M%S')}",
    }


def _ensure_target_form_state(target_id: str) -> dict[str, Any]:
    states = dict(st.session_state.get(FORM_STATE_KEY) or {})
    state = dict(states.get(target_id) or _default_form_state(target_id))
    states[target_id] = state
    st.session_state[FORM_STATE_KEY] = states
    return state


def _save_widget_form(target_id: str) -> None:
    state = _ensure_target_form_state(target_id)
    for field, widget_key in FORM_WIDGET_KEYS.items():
        if widget_key in st.session_state:
            value = st.session_state[widget_key]
            state[field] = list(value) if field == "scenarios" else value
    states = dict(st.session_state.get(FORM_STATE_KEY) or {})
    states[target_id] = state
    st.session_state[FORM_STATE_KEY] = states


def _load_widget_form(target_id: str) -> None:
    state = _ensure_target_form_state(target_id)
    for field, widget_key in FORM_WIDGET_KEYS.items():
        value = state[field]
        st.session_state[widget_key] = list(value) if field == "scenarios" else value
    st.session_state[ACTIVE_TARGET_KEY] = target_id


def _sync_form_field(target_id: str, field: str) -> None:
    widget_key = FORM_WIDGET_KEYS[field]
    states = dict(st.session_state.get(FORM_STATE_KEY) or {})
    state = dict(states.get(target_id) or _default_form_state(target_id))
    value = st.session_state.get(widget_key, state[field])
    state[field] = list(value) if field == "scenarios" else value
    states[target_id] = state
    st.session_state[FORM_STATE_KEY] = states


def _switch_target_form() -> None:
    previous_target = st.session_state.get(ACTIVE_TARGET_KEY)
    if previous_target:
        _save_widget_form(str(previous_target))
    target_id = str(st.session_state[TARGET_KEY])
    _load_widget_form(target_id)
    _load_candidate_selection(target_id)


def _load_candidate_selection(target_id: str) -> None:
    selections = dict(st.session_state.get(CANDIDATE_SELECTIONS_KEY) or {})
    st.session_state[CANDIDATE_WIDGET_KEY] = list(selections.get(target_id) or [])
    st.session_state[CANDIDATE_ACTIVE_TARGET_KEY] = target_id


def _sync_candidate_selection(target_id: str) -> None:
    selections = dict(st.session_state.get(CANDIDATE_SELECTIONS_KEY) or {})
    selections[target_id] = list(st.session_state.get(CANDIDATE_WIDGET_KEY) or [])
    st.session_state[CANDIDATE_SELECTIONS_KEY] = selections


def render_oe_capacity() -> None:
    st.subheader("基因级 OE 容量对照")
    st.caption(
        "明确区分 reaction proxy、相对未校准 OE 与绝对容量可用性。"
        "结果是模型内相对比较，不预测 mg/L、真实表达倍数或实验成功率。"
    )
    targets = list_oe_capacity_targets()
    if not targets:
        st.error("未找到可用的 hLF / OPN 内置目标。")
        return
    labels = {row["target_id"]: f"{row['label']} · {row['target_id']}" for row in targets}
    target_ids = list(labels)
    if st.session_state.get(TARGET_KEY) not in target_ids:
        st.session_state[TARGET_KEY] = target_ids[0]
    if st.session_state.get(ACTIVE_TARGET_KEY) not in target_ids:
        _load_widget_form(str(st.session_state[TARGET_KEY]))
    if st.session_state.get(CANDIDATE_ACTIVE_TARGET_KEY) not in target_ids:
        _load_candidate_selection(str(st.session_state[TARGET_KEY]))
    target_id = st.selectbox(
        "目标蛋白",
        target_ids,
        format_func=lambda value: labels[value],
        key=TARGET_KEY,
        on_change=_switch_target_form,
    )
    st.info(
        "同工酶不会放宽整个反应；单个复合体亚基不会自动提高完整复合体容量；"
        "外部证据不会覆盖当前模型映射。"
    )
    st.caption(
        "正式 gene-capacity 需要与当前 target/context/model 精确匹配、使用 model_flux 单位且"
        "带来源和审核信息的 baseline capacity。缺少 reviewed_baseline_capacity 时只展示"
        "边界与旧 reaction proxy，不会从最优 flux 或通用上界推断容量。"
    )

    _render_external_candidate_review(target_id)

    gene_id = st.text_input(
        "模型 gene ID",
        key=FORM_WIDGET_KEYS["gene_id"],
        on_change=_sync_form_field,
        args=(target_id, "gene_id"),
        help="必须使用当前模型 gene ID；common name 不会自动关联。",
    ).strip()
    product_mode_label = st.radio(
        "产品层级",
        ("相对未校准", "Reaction proxy", "绝对容量研究"),
        horizontal=True,
        key=FORM_WIDGET_KEYS["product_mode"],
        on_change=_sync_form_field,
        args=(target_id, "product_mode"),
        help="绝对容量研究在缺少审核 baseline anchor 时会返回 unavailable，且不会调用求解器。",
    )
    left, right = st.columns(2)
    with left:
        dose_mode = st.radio(
            "OE 剂量输入",
            ("明确倍数", "类别输入（仅解释）"),
            horizontal=True,
            key=FORM_WIDGET_KEYS["dose_mode"],
            on_change=_sync_form_field,
            args=(target_id, "dose_mode"),
        )
        multiplier = st.number_input(
            "表达容量倍数",
            min_value=0.1,
            max_value=20.0,
            step=0.1,
            disabled=dose_mode != "明确倍数",
            key=FORM_WIDGET_KEYS["multiplier"],
            on_change=_sync_form_field,
            args=(target_id, "multiplier"),
        )
    with right:
        scenarios = st.multiselect(
            "参数不确定性场景",
            ("low", "nominal", "high"),
            key=FORM_WIDGET_KEYS["scenarios"],
            on_change=_sync_form_field,
            args=(target_id, "scenarios"),
        )
        compare_proxy = st.checkbox(
            "同时运行旧 reaction proxy 对照",
            key=FORM_WIDGET_KEYS["compare_proxy"],
            on_change=_sync_form_field,
            args=(target_id, "compare_proxy"),
        )
    if not scenarios:
        st.warning("至少选择一个参数场景。")

    preview_col, run_col = st.columns(2)
    with preview_col:
        if st.button("预览 mapping 与参数", use_container_width=True):
            if not gene_id:
                st.error("请填写模型 gene ID。")
            else:
                with st.spinner("正在建立当前模型 mapping catalog…"):
                    try:
                        preview = preview_oe_capacity_candidate(
                            target_id=target_id,
                            gene_id=gene_id,
                            dose_payload=(
                                {
                                    "dose_id": f"{float(multiplier):g}x-preview",
                                    "dose_mode": "explicit_multiplier",
                                    "expression_multiplier": float(multiplier),
                                }
                                if dose_mode == "明确倍数"
                                else {
                                    "dose_id": "categorical_oe_preview",
                                    "dose_mode": "categorical_only",
                                    "promoter": "unspecified",
                                }
                            ),
                            product_mode={
                                "相对未校准": "relative_uncalibrated",
                                "Reaction proxy": "reaction_proxy",
                                "绝对容量研究": "absolute_capacity",
                            }[product_mode_label],
                        )
                    except Exception as exc:
                        st.error(f"预览失败：{exc}")
                    else:
                        previews = dict(st.session_state.get(LAST_PREVIEWS_KEY) or {})
                        previews[target_id] = preview
                        st.session_state[LAST_PREVIEWS_KEY] = previews
    with run_col:
        run_name = st.text_input(
            "运行名称",
            key=FORM_WIDGET_KEYS["run_name"],
            on_change=_sync_form_field,
            args=(target_id, "run_name"),
        )
        if st.button(
            "运行所选 OE 产品模式",
            type="primary",
            use_container_width=True,
            disabled=not gene_id or not scenarios,
        ):
            dose_payload = (
                {
                    "dose_id": f"{float(multiplier):g}x",
                    "dose_mode": "explicit_multiplier",
                    "expression_multiplier": float(multiplier),
                }
                if dose_mode == "明确倍数"
                else {
                    "dose_id": "categorical_oe",
                    "dose_mode": "categorical_only",
                    "promoter": "unspecified",
                }
            )
            with st.spinner("正在运行 pcSec baseline 与容量对照；大模型可能需要数十秒…"):
                try:
                    result = submit_oe_capacity_screen(
                        target_id=target_id,
                        gene_ids=(gene_id,),
                        dose_payload=dose_payload,
                        parameter_scenarios=tuple(scenarios),
                        execution_mode="comparison",
                        product_mode={
                            "相对未校准": "relative_uncalibrated",
                            "Reaction proxy": "reaction_proxy",
                            "绝对容量研究": "absolute_capacity",
                        }[product_mode_label],
                        feature_enabled=True,
                        compare_proxy=compare_proxy,
                        run_name=run_name,
                    )
                except FileExistsError as exc:
                    st.warning(f"运行名称冲突：{exc}")
                except Exception as exc:
                    st.error(f"运行失败：{exc}")
                else:
                    runs = dict(st.session_state.get(LAST_RUNS_KEY) or {})
                    runs[target_id] = result
                    st.session_state[LAST_RUNS_KEY] = runs
                    st.success("OE 产品分层运行已完成，结果已写入 local_runs/oe_capacity/ui_runs。")

    preview = (st.session_state.get(LAST_PREVIEWS_KEY) or {}).get(target_id)
    if preview:
        _render_preview(preview)
    result = (st.session_state.get(LAST_RUNS_KEY) or {}).get(target_id)
    if result:
        _render_result(result)
    _render_history(target_id)


def _render_external_candidate_review(target_id: str) -> None:
    review = load_oe_capacity_candidate_review(
        DEFAULT_OE_CAPACITY_CANDIDATE_ROOT,
        target_id=target_id,
    )
    with st.expander("Round 6A 外部 baseline capacity 候选审核", expanded=False):
        st.caption(str(review.get("message") or ""))
        st.code(str(review.get("candidate_root") or DEFAULT_OE_CAPACITY_CANDIDATE_ROOT))
        if not review.get("available"):
            st.info("候选 cache 缺失或校验失败时，当前正式容量资产不会改变。")
            return
        st.caption(
            f"candidate manifest SHA-256: {review.get('candidate_manifest_sha256')} · "
            f"formal asset SHA-256: {review.get('formal_asset_sha256')}"
        )
        candidates = list(review.get("candidates") or [])
        if not candidates:
            st.warning("当前 target/context 没有适用候选；不会用 1000、最优 flux、固定 1.0 或 fixture 补齐。")
            return
        rows = [
            {
                "candidate_id": item.get("candidate_id"),
                "scope": item.get("applicability_scope"),
                "status": item.get("status"),
                "confidence": item.get("confidence"),
                "nominal_capacity": item.get("nominal_capacity"),
                "unit": item.get("unit"),
                "gene_id": ((item.get("model_bindings") or [{}])[0]).get("gene_id"),
                "formation_handle": ((item.get("model_bindings") or [{}])[0]).get(
                    "formation_or_dilution_reaction_id"
                ),
                "conflicts": ", ".join(item.get("conflicts") or []),
                "missing_information": ", ".join(item.get("missing_information") or []),
                "promotion_eligible": item.get("promotion_eligible"),
            }
            for item in candidates
        ]
        st.dataframe(rows, use_container_width=True, hide_index=True)
        with st.expander("来源、原始单位与转换链", expanded=False):
            st.json(
                {
                    "sources": review.get("sources") or [],
                    "candidates": candidates,
                }
            )
        eligible_ids = [
            str(item.get("candidate_id"))
            for item in candidates
            if item.get("promotion_eligible")
        ]
        previous = [
            value
            for value in st.session_state.get(CANDIDATE_WIDGET_KEY, [])
            if value in eligible_ids
        ]
        if previous != st.session_state.get(CANDIDATE_WIDGET_KEY, []):
            st.session_state[CANDIDATE_WIDGET_KEY] = previous
            _sync_candidate_selection(target_id)
        selected = st.multiselect(
            "选择进入 promotion 预览的候选",
            eligible_ids,
            key=CANDIDATE_WIDGET_KEY,
            on_change=_sync_candidate_selection,
            args=(target_id,),
            help="候选仍不是 reviewed anchor；此处只选择审核对象。",
        )
        if st.button("生成只读 promotion 预览", disabled=not selected):
            try:
                preview = preview_oe_capacity_promotion(
                    candidate_root=str(review["candidate_root"]),
                    candidate_ids=tuple(selected),
                    target_id=target_id,
                    expected_candidate_manifest_sha256=str(
                        review["candidate_manifest_sha256"]
                    ),
                    expected_asset_sha256=str(review["formal_asset_sha256"]),
                )
            except Exception as exc:
                st.error(f"promotion 预览失败：{exc}")
            else:
                previews = dict(st.session_state.get(CANDIDATE_PREVIEWS_KEY) or {})
                previews[target_id] = preview
                st.session_state[CANDIDATE_PREVIEWS_KEY] = previews
        preview = (st.session_state.get(CANDIDATE_PREVIEWS_KEY) or {}).get(target_id)
        if not preview:
            return
        st.json(preview)
        reviewer = st.text_input(
            "审核人",
            key=f"oe_capacity_reviewer_{target_id}",
            help="正式 promotion 会记录审核人和时间。",
        )
        explicit = st.checkbox(
            "我明确批准以上候选写入正式容量资产",
            key=f"oe_capacity_explicit_approval_{target_id}",
        )
        if st.button(
            "正式提升为 reviewed capacity anchor",
            type="primary",
            disabled=not explicit or not reviewer.strip() or not bool(preview.get("eligible")),
        ):
            try:
                promoted = promote_oe_capacity_candidate_selection(
                    candidate_root=str(review["candidate_root"]),
                    candidate_ids=tuple(preview.get("candidate_ids") or ()),
                    reviewer=reviewer,
                    expected_candidate_manifest_sha256=str(
                        preview["candidate_manifest_sha256"]
                    ),
                    expected_asset_sha256=str(preview["formal_asset_sha256"]),
                    explicit_approval=True,
                )
            except Exception as exc:
                st.error(f"正式 promotion 失败：{exc}")
            else:
                st.success(
                    "容量资产已原子更新；正式 acceptance 尚未自动启动，请按 Round 6B runner 重验收。"
                )
                st.json(promoted)


def _render_preview(preview: dict[str, Any]) -> None:
    st.markdown("### Mapping 与参数预览")
    product = preview.get("product") or {}
    if product:
        state_cols = st.columns(4)
        state_cols[0].metric("产品状态", str(product.get("product_state") or ""))
        state_cols[1].metric("校准状态", str(product.get("calibration_status") or ""))
        state_cols[2].metric(
            "绝对容量",
            str(product.get("absolute_capacity_availability") or ""),
        )
        state_cols[3].metric(
            "绝对求解门禁",
            "允许" if product.get("absolute_solver_allowed") else "禁止",
        )
        for warning in product.get("warnings") or []:
            st.warning(str(warning))
        st.caption(
            "模型指纹："
            + str(product.get("model_fingerprint") or preview.get("model_fingerprint") or "")
        )
    metrics = st.columns(3)
    metrics[0].metric("当前基因 mapping", int(preview.get("mapping_count") or 0))
    metrics[1].metric("可用参数集", int(preview.get("parameter_set_count") or 0))
    metrics[2].metric(
        "全 catalog 可执行 mapping",
        int(preview.get("executable_mapping_count") or 0),
    )
    mappings = preview.get("mappings") or []
    if mappings:
        st.dataframe(mappings, use_container_width=True, hide_index=True)
    else:
        st.warning("当前 gene ID 没有模型内 mapping；不会自动使用 common name 或外部同源覆盖。")
    with st.expander("参数来源与 low / nominal / high 区间", expanded=bool(mappings)):
        parameter_sets = preview.get("parameter_sets") or []
        if parameter_sets:
            st.json(parameter_sets)
        else:
            st.caption("没有可用参数；核心层会返回明确的 unavailable/not-executable 状态。")


def _render_result(result: dict[str, Any]) -> None:
    st.markdown("### 最近一次当前目标运行")
    top = st.columns(4)
    top[0].metric("完成", int(result.get("completed_count") or 0))
    top[1].metric("失败 / 不可执行", int(result.get("failure_count") or 0))
    top[2].metric("目标", str(result.get("target_id") or ""))
    top[3].metric("运行状态", str(result.get("status") or "completed"))
    st.caption(
        "产品状态："
        + ", ".join(str(item) for item in result.get("product_states") or [])
        + " · 模型指纹："
        + str(result.get("model_fingerprint") or "")
    )
    rows = result.get("rows") or []
    failures = result.get("failures") or []
    if rows:
        st.dataframe(rows, use_container_width=True, hide_index=True)
    if failures:
        st.warning("以下候选未完成所选产品模式；请查看产品状态、missing_information 和 warnings。")
        st.dataframe(failures, use_container_width=True, hide_index=True)
    solver_evidence = _solver_evidence_rows(result)
    if solver_evidence:
        with st.expander(
            "逐场景与 proxy solver 证据",
            expanded=any(not bool(row.get("success")) for row in solver_evidence),
        ):
            st.dataframe(solver_evidence, use_container_width=True, hide_index=True)
    for warning in result.get("warnings") or []:
        st.warning(str(warning))
    report = export_oe_capacity_report(str(result.get("run_dir") or ""))
    if report:
        st.download_button(
            "下载 OE capacity Markdown 报告",
            data=report,
            file_name=f"{result.get('run_name', 'oe-capacity')}.md",
            mime="text/markdown",
        )
    st.caption("本页不会自动修改 recommendation tier、模型资产或实验结论。")


def _solver_evidence_rows(result: dict[str, Any]) -> list[dict[str, Any]]:
    evidence: list[dict[str, Any]] = []
    for row in (*(result.get("rows") or []), *(result.get("failures") or [])):
        identity = {
            "gene_id": row.get("gene_id"),
            "screen_status": row.get("screen_status"),
        }
        for scenario in row.get("scenario_results") or []:
            scenario_id = scenario.get("parameter_scenario")
            for phase in ("baseline", "perturbed"):
                snapshot = scenario.get(phase) or {}
                evidence.append(
                    {
                        **identity,
                        "result_type": f"scenario_{phase}",
                        "candidate": scenario_id,
                        "success": snapshot.get("success"),
                        "solver_status": snapshot.get("solver_status"),
                        "objective": snapshot.get("secretion_objective"),
                        "message": snapshot.get("message"),
                        "failure_reason": scenario.get("failure_reason"),
                    }
                )
        for snapshot in row.get("proxy_attempts") or []:
            evidence.append(
                {
                    **identity,
                    "result_type": "proxy_attempt",
                    "candidate": snapshot.get("attempt_id"),
                    "success": snapshot.get("success"),
                    "solver_status": snapshot.get("solver_status"),
                    "objective": snapshot.get("secretion_objective"),
                    "message": snapshot.get("message"),
                    "failure_reason": (
                        "" if snapshot.get("success") else snapshot.get("message")
                    ),
                }
            )
    return evidence


def _render_history(target_id: str) -> None:
    st.markdown("### 当前目标历史运行")
    runs = list_oe_capacity_runs(target_id=target_id)
    if not runs:
        st.caption("尚无当前目标的 UI 运行记录。")
        return
    st.dataframe(
        [
            {
                "run_name": row.get("run_name"),
                "target_id": row.get("target_id"),
                "context_id": row.get("context_id"),
                "model_fingerprint": row.get("model_fingerprint"),
                "product_states": ", ".join(row.get("product_states") or []),
                "absolute_capacity_available": row.get("absolute_capacity_available"),
                "completed_count": row.get("completed_count"),
                "failure_count": row.get("failure_count"),
                "status": row.get("status"),
                "error_message": row.get("error_message"),
                "run_dir": row.get("run_dir"),
            }
            for row in runs
        ],
        use_container_width=True,
        hide_index=True,
    )


__all__ = ["render_oe_capacity"]
