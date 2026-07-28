"""仿真验证「基因扰动」段的信息层级。

用户 2026-07-28 反馈这段"比较混乱"：辅助面板（候选库 / 基因查找 / 输入说明）挡在输入框前面，
四个输入框平权并列，看不出该填哪个。这里锁住重排后的层级，防回退。
"""

from __future__ import annotations

import inspect

from app.ui.views import simulation_gene_inputs as gene_inputs


def _source() -> str:
    return inspect.getsource(gene_inputs.render_gene_perturbation_form)


def test_helper_panels_are_collapsed_not_dumped_inline_above_the_inputs() -> None:
    """帮手面板必须收进**默认折叠**的入口，只占一行；此前两个面板默认展开、把输入框挤下屏。

    注意它们仍必须位于输入框**之前**（不是 UX 偏好而是功能约束：面板的"加入 KO/OE 输入"要写
    pichia_draft_* 的 session_state，Streamlit 不允许在同名 key 控件实例化后再改）——
    见 test_hlf_opn_candidate_panel_is_embedded_before_gene_textareas。所以这里锁"折叠"，不锁"在后"。
    """
    source = _source()

    helper_expander = source.index('st.expander("不知道基因 ID？在这里查找 / 从候选库挑", expanded=False)')
    candidate_panel = source.index("render_hlf_opn_candidate_panel(target_id)")
    ko_gene_widget = source.index('key="pichia_draft_ko_genes"')

    assert helper_expander < candidate_panel, "两个帮手面板必须包在折叠入口里"
    assert candidate_panel < ko_gene_widget, "帮手必须仍在输入框之前（session_state 写入约束）"


def test_reaction_level_inputs_are_demoted_below_gene_inputs() -> None:
    """反应级是高级诊断入口，不该与基因输入平权并列。

    按**输入框标签**定位，不按 session key——key 名在顶部的 _has_text() 预填探测里先出现过一次。
    """
    source = _source()

    assert source.index("敲除基因（KO gene）") < source.index("敲除反应（复合体级 KO）")
    assert source.index("过表达基因（OE gene proxy）") < source.index("过表达反应（高级 / OE reaction）")
    assert "复合体 / 反应级扰动（高级）" in source


def test_prefilled_reaction_candidate_auto_expands_advanced_section() -> None:
    """关键正确性：策展库 / 复合体候选会被预填进高级区。若不自动展开，用户从筛查页点
    "在仿真验证中核实"跳来只看到两个空的基因框，会误判跳转失效。"""
    source = _source()

    assert "expanded=reaction_prefilled" in source
    assert '_has_text("pichia_draft_ko_reactions")' in source
    assert '_has_text("pichia_draft_oe_reactions")' in source


def test_has_text_detects_only_real_content() -> None:
    gene_inputs.st.session_state.clear()
    gene_inputs.st.session_state["blank"] = "   "
    gene_inputs.st.session_state["filled"] = "PAS_chr2-2_0107"

    assert gene_inputs._has_text("missing") is False
    assert gene_inputs._has_text("blank") is False
    assert gene_inputs._has_text("filled") is True


def test_all_four_prefill_session_keys_are_still_present() -> None:
    """重排不得改动 session_state key——apply_simulation_prefill 依赖它们。"""
    source = _source()

    for key in (
        "pichia_draft_ko_genes",
        "pichia_draft_ko_reactions",
        "pichia_draft_oe_genes",
        "pichia_draft_oe_reactions",
    ):
        assert key in source, key
