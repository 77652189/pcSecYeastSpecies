"""hLF/OPN 候选面板：把"不能按基因跑"的候选按复合体反应跑通。

背景（2026-07-28 查实）：69 条策展候选里只有 13 条能按基因跑，其余 56 条被显示成"仅复核/模型外"，
读起来像死胡同——但它们**全部**带着模型里真实存在的复合体反应（项目头号 hLF 杠杆
`sec_Pdi1p_complex_formation` 就在这 56 里）。不能按基因跑的原因是模型对分泌机器没有 GPR。
此前面板既不显示策展俗名、也不显示反应 ID，等于把唯一可用入口藏了起来。
"""

from __future__ import annotations

from typing import Any

from app.ui.views import hlf_opn_candidate_panel as view


def _row(**overrides: Any) -> dict[str, Any]:
    row = {
        "target_context": "hLF",
        "gene_id": "PAS_chr4_0844",
        "source_common_name": "PDI1（单独）",
        "source_category": "二硫键 (DSB)",
        "operability_status": "not_in_model",
        "recommended_intervention": "OE",
        "review_reactions": ["sec_Pdi1p_complex_formation"],
    }
    row.update(overrides)
    return row


def test_candidate_frame_surfaces_curated_common_name_and_runnable_reaction() -> None:
    frame = view._candidate_frame([_row()])

    assert frame.iloc[0]["俗名"] == "PDI1（单独）", "研究员认得的是俗名，不是位点号"
    assert frame.iloc[0]["可跑的复合体反应"] == "sec_Pdi1p_complex_formation"
    assert frame.iloc[0]["分泌环节"] == "二硫键 (DSB)"
    # 俗名必须排在位点号之前
    columns = list(frame.columns)
    assert columns.index("俗名") < columns.index("基因")


def test_not_in_model_label_says_reaction_is_runnable_not_dead_end() -> None:
    assert "可跑" in view._label("not_in_model")
    assert "可跑" in view._label("unresolved_name")
    # 不能再是那种读着像死胡同的措辞
    assert view._label("not_in_model") != "不在模型内"
    assert view._label("unresolved_name") != "命名未解析"


def test_review_reactions_group_by_recommended_intervention() -> None:
    rows = [
        _row(recommended_intervention="OE", review_reactions=["sec_Pdi1p_complex_formation"]),
        _row(recommended_intervention="KO", review_reactions=["sec_Och1p_complex_formation"]),
        _row(recommended_intervention="review_only", review_reactions=["sec_Ignored_complex_formation"]),
        _row(recommended_intervention="OE", review_reactions=["sec_Pdi1p_complex_formation"]),  # 去重
    ]

    grouped = view._review_reactions_by_intervention(rows)

    assert grouped["oe"] == ["sec_Pdi1p_complex_formation"]
    assert grouped["ko"] == ["sec_Och1p_complex_formation"]
    assert "sec_Ignored_complex_formation" not in grouped["ko"] + grouped["oe"], "review_only 不该进反应框"


def test_apply_review_reactions_writes_reaction_boxes_only_and_merges(monkeypatch) -> None:
    """只碰反应框（基因框由另一个按钮负责，有测试锁着），且不覆盖已有内容。"""

    class _FakeSt:
        def __init__(self) -> None:
            self.session_state: dict[str, Any] = {}

    fake = _FakeSt()
    fake.session_state["pichia_draft_oe_reactions"] = "sec_Existing_complex_formation"
    monkeypatch.setattr(view, "st", fake)

    added = view._apply_review_reaction_inputs(
        {"ko": ["sec_Och1p_complex_formation"], "oe": ["sec_Pdi1p_complex_formation"]}
    )

    assert added == {"ko": ["sec_Och1p_complex_formation"], "oe": ["sec_Pdi1p_complex_formation"]}
    assert fake.session_state["pichia_draft_ko_reactions"] == "sec_Och1p_complex_formation"
    assert (
        fake.session_state["pichia_draft_oe_reactions"]
        == "sec_Existing_complex_formation\nsec_Pdi1p_complex_formation"
    ), "必须 merge，不能覆盖用户已填内容"
    # 基因框一律不碰
    assert "pichia_draft_ko_genes" not in fake.session_state
    assert "pichia_draft_oe_genes" not in fake.session_state


def test_no_reactions_means_no_writes(monkeypatch) -> None:
    class _FakeSt:
        def __init__(self) -> None:
            self.session_state: dict[str, Any] = {}

    fake = _FakeSt()
    monkeypatch.setattr(view, "st", fake)

    assert view._apply_review_reaction_inputs({"ko": [], "oe": []}) == {"ko": [], "oe": []}
    assert fake.session_state == {}
