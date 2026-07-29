"""E3 映射的**策展回填面板**（ADR-007）。

给研究人员一个应用内入口：自动起草的映射摆在这里，逐条打勾/否决/补角色，改完导出成策展文件。
把"从零查 78 个复合体"降成"审 59 条勾选题"。

**为什么是导出而不是直接写盘**：策展映射是长期科学资产，落在 `Data/` 下受保护目录，按数据治理
必须由人显式提交、并声明为科学资产变更；应用运行时一律不写受保护目录（写入只落 `local_runs/`）。
所以这里给下载按钮 + 明确的放置路径，最后一步由人来做。
"""

from __future__ import annotations

import json
from typing import Any

import pandas as pd
import streamlit as st

from app.services.gene_complex_mapping_service import (
    CURATED_MAPPING_RELATIVE_PATH,
    build_draft_mapping_rows,
    load_gene_complex_mapping,
)


_ROLE_HELP = {
    "required_subunit": "必需亚基：没有它复合体就不成立",
    "replaceable_isoenzyme": "可替换同工酶：有别的基因能顶替",
    "auxiliary": "辅助：参与但不是限速/必需",
}


def _draft_frame() -> pd.DataFrame:
    """草稿由服务层提供（UI 不得直接 import 引擎——test_streamlit_ui_does_not_import_engine_directly）。"""
    return pd.DataFrame(
        [
            {
                "复合体反应": str(row.get("complex_reaction_id") or ""),
                "候选基因": str(row.get("pichia_gene_id") or ""),
                "来源俗名": str(row.get("evidence_citation") or "").replace("策展俗名 ", ""),
                "复核结论": str(row.get("review_status") or ""),
                "亚基角色": str(row.get("subunit_role") or ""),
                "化学计量": str(row.get("stoichiometry_status") or ""),
                "备注": str(row.get("note") or ""),
            }
            for row in build_draft_mapping_rows()
        ]
    )


def _rows_from_edited(edited: pd.DataFrame) -> list[dict[str, Any]]:
    return [
        {
            "pichia_gene_id": str(row["候选基因"]),
            "complex_reaction_id": str(row["复合体反应"]),
            "subunit_role": str(row["亚基角色"]),
            "stoichiometry_status": str(row["化学计量"]),
            "review_status": str(row["复核结论"]),
            "evidence_source": "curated_review",
            "evidence_citation": f"策展俗名 {row['来源俗名']}" if str(row.get("来源俗名") or "") else "",
            "note": str(row.get("备注") or ""),
        }
        for _, row in edited.iterrows()
    ]


def render_gene_complex_mapping_review() -> None:
    """策展回填入口。默认折叠——日常用不到，只有做映射策展时才展开。"""
    existing, notes = load_gene_complex_mapping()
    label = (
        f"补全「实验时对应基因」映射 · 策展回填（当前已复核 {len(existing)} 条）"
        if existing
        else "补全「实验时对应基因」映射 · 策展回填（尚无策展数据）"
    )
    with st.expander(label, expanded=False):
        st.markdown(
            "**这是做什么的**：模型能说“过表达某个复合体有效”，但答不了“实验室该动哪个基因”——"
            "因为分泌机器在模型里没有基因关联。这里把同源比对的候选自动列好，"
            "**你只需逐条判断对不对**，改完导出给负责人提交。"
        )
        st.caption(
            "草稿一律标为“待复核 / 辅助 / 计量未知”这一最保守组合，**不会自己生效**；"
            "只有被改成“已复核 + 必需亚基 + 计量已知”的条目，才允许用于"
            "“过表达这个基因＝提升该复合体容量”的判断。"
        )

        try:
            draft = _draft_frame()
        except Exception as exc:  # noqa: BLE001 - 策展面板失败不该拖垮主流程
            st.caption(f"草稿生成失败：{exc}")
            return
        if draft.empty:
            st.caption("没有可起草的映射（同源比对未给出候选基因）。")
            return

        st.markdown(f"**自动起草 {len(draft)} 条**（覆盖 {draft['复合体反应'].nunique()} 个复合体反应）")
        edited = st.data_editor(
            draft,
            width="stretch",
            hide_index=True,
            disabled=["复合体反应", "候选基因", "来源俗名"],
            column_config={
                "复核结论": st.column_config.SelectboxColumn(
                    "复核结论",
                    options=["pending_review", "reviewed", "rejected"],
                    help="reviewed=确认这个基因确实参与该复合体；rejected=同源猜错了。",
                ),
                "亚基角色": st.column_config.SelectboxColumn(
                    "亚基角色",
                    options=list(_ROLE_HELP),
                    help=" / ".join(f"{key}：{value}" for key, value in _ROLE_HELP.items()),
                ),
                "化学计量": st.column_config.SelectboxColumn(
                    "化学计量",
                    options=["unknown", "known"],
                    help="亚基配比是否已知。未知时不得声称单基因过表达能提升复合体容量。",
                ),
                "备注": st.column_config.TextColumn("备注", help="写下判断依据 / 文献出处。"),
            },
            key="gene_complex_mapping_review_editor",
        )

        reviewed = int((edited["复核结论"] == "reviewed").sum())
        rejected = int((edited["复核结论"] == "rejected").sum())
        st.caption(f"当前：已复核 {reviewed} 条、否决 {rejected} 条、待复核 {len(edited) - reviewed - rejected} 条。")

        keep = edited[edited["复核结论"] != "rejected"]
        payload = {
            "schema_version": 1,
            "mappings": _rows_from_edited(keep),
        }
        st.download_button(
            "⬇️ 导出策展映射（JSON）",
            data=json.dumps(payload, ensure_ascii=False, indent=2).encode("utf-8"),
            file_name="gene_complex_mapping.json",
            mime="application/json",
            key="gene_complex_mapping_review_download",
            type="primary",
        )
        st.caption(
            f"导出后把文件放到仓库的 `{str(CURATED_MAPPING_RELATIVE_PATH).replace(chr(92), '/')}`，"
            "重启应用即生效。（这一步由人显式提交：策展映射是长期科学资产，应用不自动写入受保护目录。）"
        )
        for note in notes:
            st.caption(note)


__all__ = ["render_gene_complex_mapping_review"]
