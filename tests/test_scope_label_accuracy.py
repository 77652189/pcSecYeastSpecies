"""筛查范围标签上的数字必须与真实策展库对得上。

2026-07-28 用户发现：界面写"约30个反应"，实际是 **61 个唯一反应 × KO/OE 两向 = 122 个候选**——
少报一半，还连带把耗时预期说小了。数字写死在文案里就会随策展库增删悄悄变陈旧，故用测试锁住。
"""

from __future__ import annotations

import re

from app.ui.views.genome_wide_screen import SCOPE_LABELS


def _curated_reaction_count() -> int:
    from pcsec_pichia.screens.genome_wide_tradeoff import catalog_reaction_candidates

    return len({candidate["reaction_id"] for candidate in catalog_reaction_candidates()})


def test_scope_label_reaction_count_matches_curated_catalog() -> None:
    label = SCOPE_LABELS["catalog"]
    numbers = {int(match) for match in re.findall(r"\d+", label)}

    assert _curated_reaction_count() in numbers, (
        f"策展 scope 标签 {label!r} 里的数字与真实唯一反应数 {_curated_reaction_count()} 对不上——"
        "策展库变了就要同步改文案。"
    )


def test_catalog_scope_tests_both_directions_per_reaction() -> None:
    """每个反应测 KO+OE 两向，所以候选行数是唯一反应数的两倍——文案不能只提反应数不提方向。"""
    from pcsec_pichia.screens.genome_wide_tradeoff import catalog_reaction_candidates

    candidates = catalog_reaction_candidates()

    assert len(candidates) == 2 * _curated_reaction_count()
    assert {candidate["intervention_type"] for candidate in candidates} == {"KO", "OE"}
    assert "KO/OE" in SCOPE_LABELS["catalog"]


def test_gene_scope_label_still_flags_the_hour_scale_cost() -> None:
    assert "小时级" in SCOPE_LABELS["gene"]
    assert "1025" in SCOPE_LABELS["gene"]
