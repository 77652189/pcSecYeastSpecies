"""KO 候选可视化的数据准备（与 OE 短名单对等 + KO 特有的分泌↔生长两维）。

重点：分泌提升与生长代价必须分开表达——只按分泌排名会把"提升一点、生长掉一半"的陷阱
candidate 排到前面；以及既有两张 KO 明细表漏掉的"轻微生长代价"档必须能被看到。
"""

from __future__ import annotations

import pandas as pd

from app.services import genome_wide_screen_analysis as analysis


def _ko(gene_id: str, ratio: float, growth: float = 1.0, process: str = "metabolic_or_other") -> dict[str, object]:
    return {
        "target_id": "hLF",
        "gene_id": gene_id,
        "intervention_type": "KO",
        "candidate_kind": "gene",
        "secretion_ratio_vs_wildtype": ratio,
        "growth_retention_ratio": growth,
        "secretory_process": process,
        "skipped_reason": None,
        "max_feasible_mu": 0.1,
    }


def test_ranked_keeps_only_secretion_wins_and_sorts_by_effect() -> None:
    frame = pd.DataFrame(
        [
            _ko("G_SMALL", 1.02),
            _ko("G_BIG", 1.20),
            _ko("G_FLAT", 1.00),  # 无提升 -> 不进排序图
            _ko("G_WORSE", 0.80),  # 变差 -> 不进排序图
        ]
    )

    view = analysis.build_ko_candidate_view(frame, "hLF")

    assert view["ranked"]["gene_id"].tolist() == ["G_BIG", "G_SMALL"]
    assert view["counts"]["secretion_up"] == 2
    assert view["counts"]["secretion_down"] == 1


def test_growth_impact_buckets_separate_clean_wins_from_traps() -> None:
    frame = pd.DataFrame(
        [
            _ko("G_CLEAN", 1.05, growth=1.0),
            _ko("G_TRAP", 1.30, growth=0.50),  # 提升最大但生长掉一半 = 陷阱，必须可区分
            _ko("G_SLIGHT", 1.04, growth=0.995),
        ]
    )

    view = analysis.build_ko_candidate_view(frame, "hLF")
    by_gene = dict(zip(view["ranked"]["gene_id"], view["ranked"]["growth_impact"]))

    assert by_gene["G_CLEAN"] == "growth_fully_retained"
    assert by_gene["G_TRAP"] == "growth_cost"
    assert by_gene["G_SLIGHT"] == "growth_slight_cost"
    assert view["counts"]["clean_wins"] == 1
    assert view["counts"]["growth_cost"] == 1
    # 陷阱 candidate 仍排在最前（效应最大）——所以着色/分档是唯一能防误读的手段
    assert view["ranked"].iloc[0]["gene_id"] == "G_TRAP"


def test_slight_growth_cost_band_is_visible_here_though_detail_tables_drop_it(tmp_path) -> None:
    """既有维度表一张要求 >=0.999、另一张要求 <0.99，(0.99, 0.999) 这档两边都不落。

    经真实 CSV 装载路径（会补齐各默认列），与页面实际走的一致。
    """
    csv_path = tmp_path / "gene_tradeoff_rows.csv"
    pd.DataFrame([{**_ko("G_GAP", 1.05, growth=0.995), "gpr_role": "single_gene", "affected_reactions": ""}]).to_csv(
        csv_path, index=False
    )
    frame = analysis.load_gene_tradeoff_csv(str(csv_path))

    view = analysis.build_ko_candidate_view(frame, "hLF")
    result = analysis.analyze_single_target(frame, "hLF")

    assert view["counts"]["slight_growth_cost"] == 1
    # 证明它确实在两张既有明细表里都不出现
    assert result.ko_clean_wins.empty
    assert result.ko_yield_up_growth_cost.empty


def test_essential_genes_without_ratio_are_excluded_not_plotted_as_zero() -> None:
    """必需基因没有比值（KO 不可行），不能被当成 0 效应画进图里。"""
    frame = pd.DataFrame(
        [
            _ko("G_OK", 1.05),
            {**_ko("G_ESSENTIAL", 1.0), "secretion_ratio_vs_wildtype": None, "max_feasible_mu": None},
        ]
    )

    view = analysis.build_ko_candidate_view(frame, "hLF")

    assert view["counts"]["ko_rows"] == 1
    assert "G_ESSENTIAL" not in view["scatter"]["gene_id"].tolist()


def test_scatter_caps_points_and_reports_truncation() -> None:
    frame = pd.DataFrame([_ko(f"G{i}", 1.0 + i / 10000.0) for i in range(50)])

    view = analysis.build_ko_candidate_view(frame, "hLF", scatter_max_points=10)

    assert len(view["scatter"]) == 10
    assert view["counts"]["scatter_truncated"] == 40
    # 保留的是效应最大的那些（截掉的都是接近 0 的中性 KO）
    assert view["scatter"].iloc[0]["gene_id"] == "G49"


def test_oe_rows_and_other_targets_are_not_mixed_in() -> None:
    frame = pd.DataFrame(
        [
            _ko("G_KO", 1.05),
            {**_ko("G_OE", 1.50), "intervention_type": "OE"},
            {**_ko("G_OTHER", 1.40), "target_id": "OPN"},
        ]
    )

    view = analysis.build_ko_candidate_view(frame, "hLF")

    assert view["ranked"]["gene_id"].tolist() == ["G_KO"]


def test_empty_and_no_win_cases_degrade_without_error() -> None:
    assert analysis.build_ko_candidate_view(pd.DataFrame(), "hLF")["available"] is False

    no_wins = analysis.build_ko_candidate_view(pd.DataFrame([_ko("G", 0.90)]), "hLF")
    assert no_wins["ranked"].empty
    assert no_wins["counts"]["secretion_up"] == 0
