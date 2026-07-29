"""自建序列库管理：信号肽 / 引导肽 / 成熟蛋白 / 组合模板的增删改查。

界面原则（用户 2026-07-28 反馈"信息不要堆叠"）：一次只呈现一件事——
先选类别（标签页）→ 看清单（表格）→ 要改再展开表单，而不是把清单和四个表单一次性铺满屏幕。
内置条目标灰、只读；自建条目可改可删。
"""

from __future__ import annotations

from typing import Any

import pandas as pd
import streamlit as st

from app.services import custom_target_library_service as library
from app.services.pichia_target_catalog_service import (
    known_leaders,
    known_mature_proteins,
    known_signal_peptides,
)


_BUILTIN_LOADERS = {
    library.KIND_SIGNAL_PEPTIDE: known_signal_peptides,
    library.KIND_LEADER: known_leaders,
    library.KIND_MATURE: known_mature_proteins,
    library.KIND_TEMPLATE: dict,  # 组合模板没有内置来源（内置目标走"快速选择"）
}

_PTM_FIELDS = (
    ("disulfide_sites", "二硫键数"),
    ("n_glycosylation_sites", "N-糖基化位点数"),
    ("o_glycosylation_sites", "O-糖基化位点数"),
)


def merged_entries(kind: str) -> dict[str, dict[str, Any]]:
    return library.merge_with_builtin(_BUILTIN_LOADERS[kind](), kind)


def _overview_frame(entries: dict[str, dict[str, Any]], kind: str) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for entry_id, entry in entries.items():
        sequence = entry.get("sequence") or entry.get("mature_sequence") or ""
        rows.append(
            {
                "来源": "内置（只读）" if entry.get("source") == "builtin" else "自建",
                "编号": entry_id,
                "名称": entry.get("label", entry_id),
                "长度(aa)": len(str(sequence)),
                "序列预览": (str(sequence)[:28] + "…") if len(str(sequence)) > 28 else str(sequence),
            }
        )
    frame = pd.DataFrame(rows)
    if frame.empty:
        return frame
    return frame.sort_values(["来源", "编号"], ascending=[False, True]).reset_index(drop=True)


def _entry_form(kind: str, entries: dict[str, dict[str, Any]]) -> None:
    editable_ids = [key for key, value in entries.items() if value.get("editable")]
    choice = st.selectbox(
        "要做什么",
        ["➕ 新建", *[f"✏️ 修改 {entry_id}" for entry_id in editable_ids]],
        key=f"library_action_{kind}",
    )
    editing_id = choice.replace("✏️ 修改 ", "") if choice.startswith("✏️") else ""
    current = entries.get(editing_id, {}) if editing_id else {}

    with st.form(key=f"library_form_{kind}"):
        col_id, col_label = st.columns([1, 2])
        entry_id = col_id.text_input("编号", value=editing_id, disabled=bool(editing_id), help="字母/数字/下划线")
        label = col_label.text_input("名称", value=str(current.get("label", "")), help="下拉框里显示的名字")

        if kind == library.KIND_TEMPLATE:
            signal_sequence = st.text_area("信号肽序列（可空）", value=str(current.get("signal_peptide_sequence", "")), height=68)
            leader_sequence = st.text_area("引导肽序列（可空）", value=str(current.get("leader_sequence", "")), height=68)
            main_sequence = st.text_area("成熟蛋白序列", value=str(current.get("mature_sequence", "")), height=110)
        else:
            signal_sequence = leader_sequence = ""
            main_sequence = st.text_area("序列", value=str(current.get("sequence", "")), height=110,
                                         help="可直接粘贴，空格和换行会自动清理")

        ptm_values: dict[str, int] = {}
        if kind in {library.KIND_MATURE, library.KIND_TEMPLATE}:
            ptm_columns = st.columns(3)
            for column, (field, field_label) in zip(ptm_columns, _PTM_FIELDS):
                ptm_values[field] = int(
                    column.number_input(field_label, min_value=0, value=int(current.get(field, 0) or 0), step=1,
                                        key=f"library_{kind}_{field}")
                )

        submitted = st.form_submit_button("保存", type="primary")

    if not submitted:
        return

    entry: dict[str, Any] = {"id": entry_id or editing_id, "label": label, **ptm_values}
    if kind == library.KIND_TEMPLATE:
        entry.update(
            signal_peptide_sequence=signal_sequence,
            leader_sequence=leader_sequence,
            mature_sequence=main_sequence,
        )
    else:
        entry["sequence"] = main_sequence

    builtin_ids = {key for key, value in entries.items() if value.get("source") == "builtin"}
    ok, problems = library.save_entry(kind, entry, builtin_ids=builtin_ids)
    if ok:
        st.success(f"已保存 {library.KIND_LABELS[kind]}「{entry['id']}」。")
        st.rerun()
    for problem in problems:
        st.error(problem)


def _delete_control(kind: str, entries: dict[str, dict[str, Any]]) -> None:
    editable_ids = [key for key, value in entries.items() if value.get("editable")]
    if not editable_ids:
        return
    col_pick, col_confirm = st.columns([2, 1])
    target = col_pick.selectbox("删除自建条目", editable_ids, key=f"library_delete_pick_{kind}")
    confirmed = col_confirm.checkbox("确认删除", key=f"library_delete_confirm_{kind}")
    if st.button("删除", key=f"library_delete_btn_{kind}", disabled=not confirmed):
        if library.delete_entry(kind, target):
            st.success(f"已删除「{target}」。")
            st.rerun()
        st.error("删除失败：该条目可能已不存在。")


def render_target_library_manager() -> None:
    with st.expander("管理序列库（新建 / 修改 / 删除自己的信号肽、引导肽、蛋白与模板）", expanded=False):
        st.caption(
            "内置条目来自项目正式资产，**只读**；这里新增的条目存在本地运行目录，"
            "会立刻出现在上面的下拉框里。"
        )
        tabs = st.tabs([library.KIND_LABELS[kind] for kind in library.KINDS])
        for tab, kind in zip(tabs, library.KINDS):
            with tab:
                entries = merged_entries(kind)
                frame = _overview_frame(entries, kind)
                if frame.empty:
                    st.caption(f"还没有{library.KIND_LABELS[kind]}，用下面的表单新建一个。")
                else:
                    st.dataframe(frame, width="stretch", hide_index=True)
                st.divider()
                _entry_form(kind, entries)
                custom_count = sum(1 for value in entries.values() if value.get("editable"))
                if custom_count:
                    st.divider()
                    _delete_control(kind, entries)


__all__ = ["merged_entries", "render_target_library_manager"]
