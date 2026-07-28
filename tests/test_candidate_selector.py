"""E2 统一候选选择器（ADR-007）：基因与复合体并列可选，勾选后由系统路由到正确输入框。

核心不变量：用户**不需要知道**自己要动的东西在模型里算基因还是算反应——路由必须自动且正确。
路由错了会把复合体反应 id 填进基因框（GPR 解析不出来 → 静默无效）。
"""

from __future__ import annotations

from typing import Any

import pandas as pd

from app.ui.views import candidate_selector as selector


def _candidate(**overrides: Any) -> dict[str, Any]:
    row: dict[str, Any] = {
        "source_common_name": "PDI1（单独）",
        "source_category": "二硫键 (DSB)",
        "gene_id": "PAS_chr4_0844",
        "operability_status": "not_in_model",
        "recommended_intervention": "OE",
        "review_reactions": ["sec_Pdi1p_complex_formation"],
    }
    row.update(overrides)
    return row


def test_gene_executable_candidate_becomes_a_gene_row() -> None:
    frame = selector.build_unified_candidate_rows(
        [
            _candidate(
                source_common_name="PEP4",
                operability_status="model_ko_executable",
                recommended_intervention="KO",
                gene_id="PAS_PEP4",
                executable_ko_reactions=["some_reaction"],
                review_reactions=[],
            )
        ]
    )

    assert frame.iloc[0]["作用对象"] == selector.KIND_GENE
    assert frame.iloc[0]["模型对象"] == "PAS_PEP4", "能按基因跑就该用基因 id，不用反应 id"
    assert frame.iloc[0]["把握"] == "基因可直接跑"


def test_gene_unusable_candidate_falls_back_to_its_complex_reaction() -> None:
    """PDI1 这类：基因不在模型基因集，但复合体反应真实可跑（项目头号 hLF 杠杆）。"""
    frame = selector.build_unified_candidate_rows([_candidate()])

    assert frame.iloc[0]["作用对象"] == selector.KIND_COMPLEX
    assert frame.iloc[0]["模型对象"] == "sec_Pdi1p_complex_formation"
    assert "待复核" in frame.iloc[0]["把握"], "必须点明基因归属还没坐实"
    assert frame.iloc[0]["候选"] == "PDI1（单独）", "显示研究员认得的俗名"


def test_routing_sends_genes_and_complexes_to_different_boxes() -> None:
    selected = pd.DataFrame(
        [
            {"改造方式": "KO", "作用对象": selector.KIND_GENE, "模型对象": "PAS_KO_GENE"},
            {"改造方式": "OE", "作用对象": selector.KIND_GENE, "模型对象": "PAS_OE_GENE"},
            {"改造方式": "KO", "作用对象": selector.KIND_COMPLEX, "模型对象": "sec_KO_complex_formation"},
            {"改造方式": "OE", "作用对象": selector.KIND_COMPLEX, "模型对象": "sec_OE_complex_formation"},
        ]
    )

    routed = selector.route_selected_candidates(selected)

    assert routed["pichia_draft_ko_genes"] == ["PAS_KO_GENE"]
    assert routed["pichia_draft_oe_genes"] == ["PAS_OE_GENE"]
    assert routed["pichia_draft_ko_reactions"] == ["sec_KO_complex_formation"]
    assert routed["pichia_draft_oe_reactions"] == ["sec_OE_complex_formation"]


def test_review_only_and_dead_candidates_are_not_selectable() -> None:
    """选了也跑不了的不该出现在列表里。"""
    frame = selector.build_unified_candidate_rows(
        [
            _candidate(recommended_intervention="review_only"),
            _candidate(source_common_name="无入口", review_reactions=[], gene_id=""),
        ]
    )

    assert frame.empty


def test_candidates_sharing_one_complex_reaction_collapse_to_one_row() -> None:
    """PDI1/ERO1/ERV2 共用一个复合体反应——不去重会出现几行完全等效的候选。"""
    shared = "sec_PDI1_ERV2_Ero1p_complex_formation"
    frame = selector.build_unified_candidate_rows(
        [
            _candidate(source_common_name="PDI1", review_reactions=[shared]),
            _candidate(source_common_name="ERO1", review_reactions=[shared]),
            _candidate(source_common_name="ERV2", review_reactions=[shared]),
        ]
    )

    assert len(frame) == 1
    assert frame.iloc[0]["模型对象"] == shared


def test_apply_merges_into_session_state_without_overwriting(monkeypatch) -> None:
    class _FakeSt:
        def __init__(self) -> None:
            self.session_state: dict[str, Any] = {}

    fake = _FakeSt()
    fake.session_state["pichia_draft_oe_genes"] = "PAS_ALREADY_THERE"
    monkeypatch.setattr(selector, "st", fake)

    added = selector.apply_routed_candidates(
        {
            "pichia_draft_ko_genes": [],
            "pichia_draft_oe_genes": ["PAS_NEW"],
            "pichia_draft_ko_reactions": [],
            "pichia_draft_oe_reactions": ["sec_New_complex_formation"],
        }
    )

    assert added == 2
    assert fake.session_state["pichia_draft_oe_genes"] == "PAS_ALREADY_THERE\nPAS_NEW"
    assert fake.session_state["pichia_draft_oe_reactions"] == "sec_New_complex_formation"
    assert "pichia_draft_ko_genes" not in fake.session_state, "空的不要写"


def test_empty_selection_routes_to_nothing() -> None:
    routed = selector.route_selected_candidates(pd.DataFrame())

    assert all(not items for items in routed.values())


def test_selector_renders_before_text_areas_in_the_form() -> None:
    """session_state 约束：选择器要写 pichia_draft_*，必须在同名 key 控件实例化之前。"""
    import inspect

    from app.ui.views import simulation_gene_inputs

    source = inspect.getsource(simulation_gene_inputs.render_gene_perturbation_form)

    assert source.index("render_candidate_selector") < source.index('key="pichia_draft_ko_genes"')
