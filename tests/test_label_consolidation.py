"""枚举释义的单一来源约束。

背景（2026-07-29 收拢前的实测）：14 个本地枚举字典散在各视图里，28 个键重复定义，
其中 `manual_review_required` 在三处是**三种说法**（"需要人工检查" / "需人工确认：有部分数据库证据"
/ "需人工复核"）。扁平字典装不下同码多义，所以收拢方案是"中央字典 + 显式按域覆盖"，
而不是把所有条目倒进一个 dict——后者会产生错译。
"""

from __future__ import annotations

import pytest

from app.core.i18n import DOMAIN_VALUE_LABELS, SIMULATION_RESULT_VALUE_LABELS, sim_result_value_label
from app.ui.views import simulation_display as display
from app.ui.views import simulation_gene_inputs as gene_inputs


_DERIVED_DICTS = {
    "TARGET_SEMANTICS_LABELS": display.TARGET_SEMANTICS_LABELS,
    "MAPPING_LEVEL_LABELS": display.MAPPING_LEVEL_LABELS,
    "MAPPING_CONFIDENCE_LABELS": display.MAPPING_CONFIDENCE_LABELS,
    "GPR_ROLE_LABELS": display.GPR_ROLE_LABELS,
    "CAPACITY_EFFECT_LABELS": display.CAPACITY_EFFECT_LABELS,
    "SIMULATION_BASIS_LABELS": display.SIMULATION_BASIS_LABELS,
    "KO_SUPPORT_STATUS_LABELS": display.KO_SUPPORT_STATUS_LABELS,
    "OE_SUPPORT_STATUS_LABELS": display.OE_SUPPORT_STATUS_LABELS,
    "WET_LAB_READINESS_LABELS": display.WET_LAB_READINESS_LABELS,
}


@pytest.mark.parametrize("name", sorted(_DERIVED_DICTS))
def test_every_derived_label_resolves_to_real_chinese(name: str) -> None:
    """派生字典的值不能悄悄回退成原始英文码——那说明中央字典漏登记了这个键。"""
    fell_back = [code for code, label in _DERIVED_DICTS[name].items() if label == code]

    assert fell_back == [], f"{name} 里这些码在中央字典查不到：{fell_back}"


def test_context_dependent_codes_keep_their_distinct_meanings() -> None:
    """一码多义必须靠按域覆盖保住，不能被扁平合并抹平。"""
    assert display.CAPACITY_EFFECT_LABELS["manual_review_required"] == "需要人工检查"
    assert "数据库证据" in display.WET_LAB_READINESS_LABELS["manual_review_required"]
    assert display.GPR_ROLE_LABELS["single_gene"] == "单基因酶"
    # 不带域时给通用释义
    assert sim_result_value_label("manual_review_required") == "需人工复核"
    assert sim_result_value_label("single_gene") == "单基因"


def test_mapping_labels_are_no_longer_defined_twice() -> None:
    """这两个字典曾在 simulation_display 与 simulation_gene_inputs 各写一份（内容当时相同、迟早分叉）。"""
    assert gene_inputs.MAPPING_LEVEL_LABELS is display.MAPPING_LEVEL_LABELS
    assert gene_inputs.MAPPING_CONFIDENCE_LABELS is display.MAPPING_CONFIDENCE_LABELS


def test_domain_overrides_are_registered_centrally_not_in_views() -> None:
    """按域覆盖只能登记在 app.core.i18n，视图里不得再私开字典。"""
    assert set(DOMAIN_VALUE_LABELS) == {"gpr_role", "capacity_effect", "wet_lab_readiness"}
    for domain, mapping in DOMAIN_VALUE_LABELS.items():
        assert mapping, f"{domain} 覆盖表不该为空"
        for code, label in mapping.items():
            assert label != code, f"{domain}.{code} 覆盖值不能等于原码"


@pytest.mark.parametrize("name", sorted(_DERIVED_DICTS))
def test_enum_dicts_are_derived_not_hardcoded(name: str) -> None:
    """防回退：这些**枚举释义**字典必须由中央字典派生，不能又写回本地字面量。

    只管枚举值字典；DataFrame 的列名映射（input_gene_id -> 输入基因）是另一回事，
    它们描述的是表格结构而非 payload 枚举，留在视图里是合理的。
    """
    import inspect
    import re

    source = inspect.getsource(display)
    match = re.search(rf"^{name} = \{{(.*?)^\}}", source, re.MULTILINE | re.DOTALL)

    assert match, f"没找到 {name} 的定义"
    body = match.group(1)
    assert "sim_result_value_label" in body, f"{name} 必须从中央字典派生"
    hardcoded = re.findall(r':\s*"[^"]*[一-鿿]', body)
    assert hardcoded == [], f"{name} 里又出现了硬编码中文：{hardcoded[:3]}"


def test_central_dictionary_absorbed_the_scattered_entries() -> None:
    """收拢后中央字典应显著变大（收编前 136 条，收编了约 79 个新键）。"""
    assert len(SIMULATION_RESULT_VALUE_LABELS) >= 210
