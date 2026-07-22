from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

import streamlit as st

from app.services.pichia_shadow_cross_check_service import (
    load_pichia_shadow_cross_check_manifest,
    run_pichia_shadow_cross_check,
)


SHADOW_CROSS_CHECK_STATE_KEY = "shadow_lp_cross_check_last_response"
SHADOW_CROSS_CHECK_MANIFEST_KEY = "shadow_lp_cross_check_manifest_payload"

_SHADOW_STATUS_LABELS = {"ok": "一致", "review_required": "需人工复核"}


def render_shadow_cross_check() -> None:
    st.header("Shadow LP一致性验证")
    _ensure_state()
    run_tab, manifest_tab = st.tabs(["运行", "查看"])
    with run_tab:
        target_id = st.selectbox("目标", ["hLF", "OPN_ALPHA_FULL_PROJECT"], key="shadow_lp_cross_check_target")
        screen_run_id = st.text_input("筛查运行 ID", value="", key="shadow_lp_cross_check_screen_run_id")
        saved_result_path = st.text_input("已保存结果路径", value="", key="shadow_lp_cross_check_saved_result_path")
        output_dir = st.text_input("输出目录", value="", key="shadow_lp_cross_check_output_dir")
        if st.button("运行 Shadow LP cross-check", key="shadow_lp_cross_check_run"):
            try:
                response = run_pichia_shadow_cross_check(
                    target_id=target_id,
                    screen_run_id=screen_run_id,
                    saved_result_path=_optional_path(saved_result_path),
                    output_dir=_optional_path(output_dir),
                )
                st.session_state[SHADOW_CROSS_CHECK_STATE_KEY] = response
                _render_status(response)
            except Exception as exc:  # pragma: no cover - Streamlit surface should show service failures.
                st.error(f"{type(exc).__name__}: {exc}")
        _render_last_response(st.session_state.get(SHADOW_CROSS_CHECK_STATE_KEY))
    with manifest_tab:
        manifest_path = st.text_input("manifest 路径", value="", key="shadow_lp_cross_check_manifest_path")
        if st.button("读取 manifest", key="shadow_lp_cross_check_load_manifest"):
            try:
                payload = load_pichia_shadow_cross_check_manifest(Path(manifest_path))
                st.session_state[SHADOW_CROSS_CHECK_MANIFEST_KEY] = payload
            except Exception as exc:  # pragma: no cover - Streamlit surface should show service failures.
                st.error(f"{type(exc).__name__}: {exc}")
        _render_manifest(st.session_state.get(SHADOW_CROSS_CHECK_MANIFEST_KEY))


def _ensure_state() -> None:
    st.session_state.setdefault(SHADOW_CROSS_CHECK_STATE_KEY, None)
    st.session_state.setdefault(SHADOW_CROSS_CHECK_MANIFEST_KEY, None)


def _render_status(response: Mapping[str, Any]) -> None:
    if response.get("status") == "ok":
        st.success("一致（shadow LP 与参考结果对齐）")
    else:
        status = str(response.get("status", "review_required"))
        st.warning(_SHADOW_STATUS_LABELS.get(status, status))


def _render_last_response(response: Mapping[str, Any] | None) -> None:
    if not response:
        return
    st.subheader("最近一次结果")
    cols = st.columns(3)
    status = str(response.get("status", ""))
    cols[0].metric("状态", _SHADOW_STATUS_LABELS.get(status, status))
    cols[1].metric("在容差内", "是" if response.get("within_tolerance") else "否")
    cols[2].metric("相对差异", _display_value(response.get("relative_diff")))
    st.code(str(response.get("report_path", "")))
    warnings = response.get("warnings") or []
    if warnings:
        st.warning("; ".join(str(item) for item in warnings))


def _render_manifest(payload: Mapping[str, Any] | None) -> None:
    if not payload:
        return
    result = payload.get("result", {}) if isinstance(payload, Mapping) else {}
    st.subheader("清单内容")
    st.json(result)
    statement = payload.get("no_experimental_claim_statement", "")
    if statement:
        st.markdown(str(statement))


def _optional_path(value: str) -> Path | None:
    text = str(value or "").strip()
    return Path(text) if text else None


def _display_value(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)
