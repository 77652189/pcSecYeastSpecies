"""E3 映射的**策展回填面板**（ADR-007）。

给研究人员一个应用内入口：同源比对的候选自动列好，逐条判断对不对，**点确认即保存生效**。

设计要点（2026-07-28 按用户反馈重做）：
- **不做"导出再导入"**——那是多余摩擦。复核结果直接存 `local_runs/` 工作副本并立即生效，
  沿用项目既有模式（实验反馈也是先落 local_runs、人工确认后才提升到 `Data/`）。
  应用不自动写受保护的 `Data/`，但完全可以写工作区。
- **按模型预测效应排序**：先审最值钱的（+8.15% 的 PDI1/ERO1 轴排最前），而不是从字母序第一条开始。
- **写入前二次确认**：勾选确认才落盘，避免误点。
"""

from __future__ import annotations

from typing import Any

import pandas as pd
import streamlit as st

from app.services.gene_complex_mapping_service import (
    CURATED_MAPPING_RELATIVE_PATH,
    build_draft_mapping_rows,
    load_gene_complex_mapping,
    save_reviewed_mappings,
    working_gene_complex_mapping_path,
)

REVIEW_CHOICES = ["待复核", "确认参与", "否决（同源猜错）"]
ROLE_CHOICES = ["辅助/不确定", "必需亚基", "可替换同工酶"]
STOICHIOMETRY_CHOICES = ["未知", "已知"]

_REVIEW_TO_CONTRACT = {"待复核": "pending_review", "确认参与": "reviewed", "否决（同源猜错）": "rejected"}
_ROLE_TO_CONTRACT = {"辅助/不确定": "auxiliary", "必需亚基": "required_subunit", "可替换同工酶": "replaceable_isoenzyme"}
_STOICHIOMETRY_TO_CONTRACT = {"未知": "unknown", "已知": "known"}
_CONTRACT_TO_REVIEW = {value: key for key, value in _REVIEW_TO_CONTRACT.items()}
_CONTRACT_TO_ROLE = {value: key for key, value in _ROLE_TO_CONTRACT.items()}
_CONTRACT_TO_STOICHIOMETRY = {value: key for key, value in _STOICHIOMETRY_TO_CONTRACT.items()}


def _existing_by_key(paths: Any | None = None) -> dict[tuple[str, str], Any]:
    rows, _ = load_gene_complex_mapping(paths)
    return {(row.pichia_gene_id, row.complex_reaction_id): row for row in rows}


def build_review_frame(target_id: str = "hLF") -> pd.DataFrame:
    """待审队列：草稿 + 已保存的复核结果，**按模型预测效应降序**（先审最值钱的）。"""
    existing = _existing_by_key()
    try:
        from app.services.screen_effect_lookup import load_screen_effect_lookup

        effects = load_screen_effect_lookup(target_id)
    except Exception:  # noqa: BLE001 - 效应只是排序依据，读不到就按原顺序
        effects = {}

    rows: list[dict[str, Any]] = []
    for draft in build_draft_mapping_rows():
        gene_id = str(draft.get("pichia_gene_id") or "")
        reaction_id = str(draft.get("complex_reaction_id") or "")
        saved = existing.get((gene_id, reaction_id))
        effect = effects.get(("OE", reaction_id)) or effects.get(("KO", reaction_id))
        rows.append(
            {
                "模型预测提升(%)": (effect[0] * 100.0) if effect else None,
                "复合体反应": reaction_id,
                "候选基因": gene_id,
                "来源俗名": str(draft.get("evidence_citation") or "").replace("策展俗名 ", ""),
                "复核结论": _CONTRACT_TO_REVIEW.get(
                    saved.review_status if saved else "", REVIEW_CHOICES[0]
                ),
                "亚基角色": _CONTRACT_TO_ROLE.get(saved.subunit_role if saved else "", ROLE_CHOICES[0]),
                "化学计量": _CONTRACT_TO_STOICHIOMETRY.get(
                    saved.stoichiometry_status if saved else "", STOICHIOMETRY_CHOICES[0]
                ),
                "判断依据": (saved.note if saved else "") or "",
            }
        )
    frame = pd.DataFrame(rows)
    if frame.empty:
        return frame
    return frame.sort_values("模型预测提升(%)", ascending=False, na_position="last").reset_index(drop=True)


def rows_to_contract_payload(edited: pd.DataFrame) -> list[dict[str, Any]]:
    """界面中文选项 → 契约字段。仍是"待复核"的不落盘（没判断过的东西不该占位）。"""
    payload: list[dict[str, Any]] = []
    for _, row in edited.iterrows():
        review = _REVIEW_TO_CONTRACT.get(str(row.get("复核结论") or ""), "pending_review")
        if review == "pending_review":
            continue
        payload.append(
            {
                "pichia_gene_id": str(row.get("候选基因") or ""),
                "complex_reaction_id": str(row.get("复合体反应") or ""),
                "subunit_role": _ROLE_TO_CONTRACT.get(str(row.get("亚基角色") or ""), "auxiliary"),
                "stoichiometry_status": _STOICHIOMETRY_TO_CONTRACT.get(str(row.get("化学计量") or ""), "unknown"),
                "review_status": review,
                "evidence_source": "curated_review",
                "evidence_citation": f"策展俗名 {row.get('来源俗名')}" if str(row.get("来源俗名") or "") else "",
                "note": str(row.get("判断依据") or ""),
            }
        )
    return payload


def render_gene_complex_mapping_review(target_id: str = "hLF") -> None:
    existing, _ = load_gene_complex_mapping()
    reviewed_count = sum(1 for row in existing if row.is_reviewed)
    label = (
        f"补全「实验时对应基因」映射 · 策展复核（已确认 {reviewed_count} 条）"
        if reviewed_count
        else "补全「实验时对应基因」映射 · 策展复核"
    )
    with st.expander(label, expanded=False):
        st.markdown(
            "模型能说“过表达某个复合体有效”，但答不了“实验室该动哪个基因”。"
            "下面是同源比对给出的候选，**按模型预测提升排好序——先审最值钱的那几条**。"
            "改完点保存即刻生效，不用导来导去。"
        )
        st.caption(
            "只有改成「确认参与 + 必需亚基 + 计量已知」的条目，才允许用于"
            "“过表达这个基因＝提升该复合体容量”的判断；其余组合只作参考、不解锁该结论。"
        )

        try:
            frame = build_review_frame(target_id)
        except Exception as exc:  # noqa: BLE001 - 策展面板失败不该拖垮主流程
            st.caption(f"待审队列生成失败：{exc}")
            return
        if frame.empty:
            st.caption("没有可复核的映射（同源比对未给出候选基因）。")
            return

        edited = st.data_editor(
            frame,
            width="stretch",
            hide_index=True,
            disabled=["模型预测提升(%)", "复合体反应", "候选基因", "来源俗名"],
            column_config={
                "模型预测提升(%)": st.column_config.NumberColumn(
                    "模型预测提升(%)", format="%.3f", help="该复合体在模型里的相对效应，用来决定先审哪条。"
                ),
                "复核结论": st.column_config.SelectboxColumn("复核结论", options=REVIEW_CHOICES, width="medium"),
                "亚基角色": st.column_config.SelectboxColumn("亚基角色", options=ROLE_CHOICES, width="medium"),
                "化学计量": st.column_config.SelectboxColumn(
                    "化学计量", options=STOICHIOMETRY_CHOICES, help="亚基配比是否已知。"
                ),
                "判断依据": st.column_config.TextColumn("判断依据", help="文献 / 数据库出处，便于他人复查。"),
            },
            key="gene_complex_mapping_review_editor",
        )

        payload = rows_to_contract_payload(edited)
        decided = len(payload)
        confirmed = sum(1 for row in payload if row["review_status"] == "reviewed")
        st.caption(
            f"本次将保存 **{decided}** 条已判断的条目（确认参与 {confirmed}、否决 {decided - confirmed}）；"
            f"仍为“待复核”的 {len(edited) - decided} 条不写入。"
        )

        agree = st.checkbox(
            "我确认以上判断，写入并立即生效",
            key="gene_complex_mapping_review_confirm",
            disabled=decided == 0,
        )
        if st.button(
            "保存复核结果",
            key="gene_complex_mapping_review_save",
            type="primary",
            disabled=not agree or decided == 0,
        ):
            saved, problems = save_reviewed_mappings(payload)
            st.success(f"已保存 {saved} 条并立即生效。")
            for problem in problems:
                st.warning(problem)
            st.rerun()

        st.caption(
            f"保存位置：`{working_gene_complex_mapping_path()}`（运行工作区，随时可改）。"
            f"要沉淀为长期科学资产，再由人显式提交到 `{str(CURATED_MAPPING_RELATIVE_PATH).replace(chr(92), '/')}`——"
            "正式资产优先于工作副本。"
        )


__all__ = ["build_review_frame", "render_gene_complex_mapping_review", "rows_to_contract_payload"]
