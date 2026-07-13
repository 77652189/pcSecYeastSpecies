from __future__ import annotations

from datetime import datetime
from typing import Any

import streamlit as st

from app.services.pichia_oe_capacity_service import (
    export_oe_capacity_report,
    list_oe_capacity_runs,
    list_oe_capacity_targets,
    preview_oe_capacity_candidate,
    submit_oe_capacity_screen,
)


TARGET_KEY = "oe_capacity_target_id"
LAST_PREVIEWS_KEY = "oe_capacity_last_previews_by_target"
LAST_RUNS_KEY = "oe_capacity_last_runs_by_target"


def render_oe_capacity() -> None:
    st.subheader("基因级 OE 容量对照")
    st.caption(
        "并列比较 baseline、旧 reaction proxy 与 gene-enzyme capacity。"
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
    target_id = st.selectbox(
        "目标蛋白",
        target_ids,
        format_func=lambda value: labels[value],
        key=TARGET_KEY,
    )
    st.info(
        "同工酶不会放宽整个反应；单个复合体亚基不会自动提高完整复合体容量；"
        "外部证据不会覆盖当前模型映射。"
    )
    st.caption(
        "若出现 nonzero_baseline_formation_flux，表示该 formation handle 的基线通量为 0；"
        "相对倍数不能凭空创建非零容量，该结果只能作为边界说明。"
    )

    gene_key = f"oe_capacity_gene_id_{target_id}"
    dose_mode_key = f"oe_capacity_dose_mode_{target_id}"
    multiplier_key = f"oe_capacity_multiplier_{target_id}"
    scenarios_key = f"oe_capacity_scenarios_{target_id}"
    compare_key = f"oe_capacity_compare_proxy_{target_id}"
    run_name_key = f"oe_capacity_run_name_{target_id}"
    default_gene = "PAS_chr2-1_0047"
    if gene_key not in st.session_state:
        st.session_state[gene_key] = default_gene
    gene_id = st.text_input(
        "模型 gene ID",
        key=gene_key,
        help="必须使用当前模型 gene ID；common name 不会自动关联。",
    ).strip()
    left, right = st.columns(2)
    with left:
        dose_mode = st.radio(
            "OE 剂量输入",
            ("明确倍数", "类别输入（仅解释）"),
            horizontal=True,
            key=dose_mode_key,
        )
        if multiplier_key not in st.session_state:
            st.session_state[multiplier_key] = 2.0
        multiplier = st.number_input(
            "表达容量倍数",
            min_value=0.1,
            max_value=20.0,
            step=0.1,
            disabled=dose_mode != "明确倍数",
            key=multiplier_key,
        )
    with right:
        if scenarios_key not in st.session_state:
            st.session_state[scenarios_key] = ["low", "nominal", "high"]
        scenarios = st.multiselect(
            "参数不确定性场景",
            ("low", "nominal", "high"),
            key=scenarios_key,
        )
        if compare_key not in st.session_state:
            st.session_state[compare_key] = True
        compare_proxy = st.checkbox(
            "同时运行旧 reaction proxy 对照",
            key=compare_key,
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
                        )
                    except Exception as exc:
                        st.error(f"预览失败：{exc}")
                    else:
                        previews = dict(st.session_state.get(LAST_PREVIEWS_KEY) or {})
                        previews[target_id] = preview
                        st.session_state[LAST_PREVIEWS_KEY] = previews
    with run_col:
        if run_name_key not in st.session_state:
            st.session_state[run_name_key] = (
                f"oe-{target_id}-{datetime.now().strftime('%Y%m%d-%H%M%S')}"
            )
        run_name = st.text_input("运行名称", key=run_name_key)
        if st.button(
            "运行 baseline / proxy / gene-capacity",
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
                    st.success("OE 容量对照已完成，结果已写入 local_runs/oe_capacity/ui_runs。")

    preview = (st.session_state.get(LAST_PREVIEWS_KEY) or {}).get(target_id)
    if preview:
        _render_preview(preview)
    result = (st.session_state.get(LAST_RUNS_KEY) or {}).get(target_id)
    if result:
        _render_result(result)
    _render_history(target_id)


def _render_preview(preview: dict[str, Any]) -> None:
    st.markdown("### Mapping 与参数预览")
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
            st.caption("没有完整的 kcat、分子量和基线容量参数，gene-capacity 将不可执行。")


def _render_result(result: dict[str, Any]) -> None:
    st.markdown("### 最近一次当前目标运行")
    top = st.columns(3)
    top[0].metric("完成", int(result.get("completed_count") or 0))
    top[1].metric("失败 / 不可执行", int(result.get("failure_count") or 0))
    top[2].metric("目标", str(result.get("target_id") or ""))
    rows = result.get("rows") or []
    failures = result.get("failures") or []
    if rows:
        st.dataframe(rows, use_container_width=True, hide_index=True)
    if failures:
        st.warning("以下候选未完成 gene-capacity 对照；请查看 missing_information 和 warnings。")
        st.dataframe(failures, use_container_width=True, hide_index=True)
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
                "completed_count": row.get("completed_count"),
                "failure_count": row.get("failure_count"),
                "run_dir": row.get("run_dir"),
            }
            for row in runs
        ],
        use_container_width=True,
        hide_index=True,
    )


__all__ = ["render_oe_capacity"]
