"""统一命名入口：反应 id → 研究员看得懂的名字。

背景：命名逻辑此前散在多处（基因名一处、复合体俗名一处、各页面各自截断 id 一处），
而模型自带的 `rxnNames`（29026 条）**根本没被加载器读进来**，所以图表上只能看到
PROTEINS_glyc / METtm_no_1_fwd 这类 id。这里锁住合并后的解析优先级与诚实标注。
"""

from __future__ import annotations

import pytest

from app.services import display_naming


def test_model_reaction_names_are_loaded_and_win_over_other_sources() -> None:
    """.mat 里的名称最权威（代谢反应几乎全覆盖）。"""
    assert display_naming.reaction_display_name("METtm_no_1_fwd") == "L-methionine mitochondrial permease"
    assert display_naming.reaction_display_name("LPPROTL_no_1_fwd") == "lipoate protein ligase"
    assert "Glycerol" in display_naming.reaction_display_name("PROTEINS_glyc")


def test_unclear_reactions_are_flagged_not_passed_through_as_a_name() -> None:
    """模型作者把一部分反应直接标成 unclear reaction——原样显示会让人以为那是个正经靶点。"""
    label = display_naming.reaction_display_name("4HPHACS")

    assert label == display_naming.UNCLEAR_REACTION_LABEL
    assert "unclear" not in label.lower()


def test_secretion_complexes_fall_back_to_curated_common_names() -> None:
    """分泌机器复合体在 .mat 里**没有**名称，可读名只存在于策展库。"""
    name = display_naming.reaction_display_name("sec_Pdi1p_complex_formation")

    assert name, "复合体反应必须能拿到策展俗名"
    assert "PDI1" in name


def test_unknown_reaction_returns_empty_not_a_guess() -> None:
    assert display_naming.reaction_display_name("NOT_A_REAL_REACTION_XYZ") == ""
    assert display_naming.reaction_display_name("") == ""


def test_label_always_keeps_the_model_id_visible() -> None:
    """id 才是模型实际算的对象，研究员要靠它跟仿真输入对上——不能只给名字。"""
    label = display_naming.reaction_display_label("METtm_no_1_fwd")

    assert "METtm_no_1_fwd" in label
    assert "L-methionine" in label


def test_label_of_unnamed_reaction_is_just_the_id() -> None:
    assert display_naming.reaction_display_label("NOT_A_REAL_REACTION_XYZ") == "NOT_A_REAL_REACTION_XYZ"


def test_long_names_are_truncated_for_charts() -> None:
    long_name = "x" * 200
    assert len(display_naming._shorten(long_name, 40)) == 40
    assert display_naming._shorten("short", 40) == "short"


@pytest.mark.parametrize(
    "view_module, function_name",
    [("app.ui.views.simulation_results", "_short_reaction_label"), ("app.ui.views.genome_wide_screen", "_short_reaction")],
)
def test_both_chart_label_helpers_delegate_to_the_single_resolver(view_module: str, function_name: str) -> None:
    """两个页面此前各有一份截断逻辑，同一个反应显示不一致。锁住它们都走统一入口。"""
    import importlib
    import inspect

    source = inspect.getsource(getattr(importlib.import_module(view_module), function_name))

    assert "display_naming" in source, f"{view_module}.{function_name} 必须走统一命名入口"
    assert "reaction_display_name" in source
