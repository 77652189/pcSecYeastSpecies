"""把筛查真实效应接进候选选择器：挑候选从"按机制猜"变成"按模型算的效应排"。

诚实边界：效应是模型内部相对量、不是产量预测；KO 必须连生长代价一起看（实测常见
"提升几个点、生长掉一半"的陷阱候选）。没跑过筛查时必须优雅降级、不能报错。
"""

from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd

from app.services import screen_effect_lookup as lookup
from app.ui.views import candidate_selector as selector


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    fields = ["target_id", "gene_id", "intervention_type", "secretion_ratio_vs_wildtype", "growth_retention_ratio"]
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def test_effects_are_keyed_by_intervention_and_model_object(tmp_path) -> None:
    """同一个对象的 KO 和 OE 是两条不同结果，不能混为一谈。"""
    csv_path = tmp_path / "gene_tradeoff_rows.csv"
    _write_csv(
        csv_path,
        [
            {"target_id": "hLF", "gene_id": "sec_Pdi1p_complex_formation", "intervention_type": "OE",
             "secretion_ratio_vs_wildtype": 1.0815, "growth_retention_ratio": 1.0},
            {"target_id": "hLF", "gene_id": "sec_Pdi1p_complex_formation", "intervention_type": "KO",
             "secretion_ratio_vs_wildtype": 0.4, "growth_retention_ratio": 0.9},
            {"target_id": "OPN", "gene_id": "sec_Pdi1p_complex_formation", "intervention_type": "OE",
             "secretion_ratio_vs_wildtype": 1.001, "growth_retention_ratio": 1.0},
        ],
    )

    effects = lookup._read_effects(str(csv_path), csv_path.stat().st_mtime, "hLF")

    assert effects[("OE", "sec_Pdi1p_complex_formation")][0] == pytest_approx(0.0815)
    assert effects[("KO", "sec_Pdi1p_complex_formation")][0] == pytest_approx(-0.6)
    # 别的靶点的行不能串进来
    assert all(value[0] != pytest_approx(0.001) for value in effects.values())


def pytest_approx(value: float, tol: float = 1e-9):
    class _Approx:
        def __eq__(self, other: object) -> bool:
            return isinstance(other, (int, float)) and abs(float(other) - value) < 1e-6

        def __ne__(self, other: object) -> bool:
            return not self.__eq__(other)

        def __repr__(self) -> str:  # pragma: no cover - 仅用于失败信息
            return f"~{value}"

    return _Approx()


def test_malformed_rows_are_skipped_not_crashing(tmp_path) -> None:
    csv_path = tmp_path / "gene_tradeoff_rows.csv"
    _write_csv(
        csv_path,
        [
            {"target_id": "hLF", "gene_id": "", "intervention_type": "OE",
             "secretion_ratio_vs_wildtype": 1.1, "growth_retention_ratio": 1.0},
            {"target_id": "hLF", "gene_id": "R1", "intervention_type": "不认识",
             "secretion_ratio_vs_wildtype": 1.1, "growth_retention_ratio": 1.0},
            {"target_id": "hLF", "gene_id": "R2", "intervention_type": "OE",
             "secretion_ratio_vs_wildtype": "", "growth_retention_ratio": 1.0},
            {"target_id": "hLF", "gene_id": "R3", "intervention_type": "OE",
             "secretion_ratio_vs_wildtype": 1.05, "growth_retention_ratio": ""},
        ],
    )

    effects = lookup._read_effects(str(csv_path), csv_path.stat().st_mtime, "hLF")

    assert set(effects) == {("OE", "R3")}, "只有可解析的行进来"
    assert effects[("OE", "R3")][1] == 1.0, "生长缺失时按不影响生长处理"


def test_missing_file_degrades_to_empty(tmp_path) -> None:
    assert lookup._read_effects(str(tmp_path / "nope.csv"), 0.0, "hLF") == {}


def test_composite_target_id_resolves_to_the_screened_target(monkeypatch) -> None:
    """回归：把默认构建模式改成三段式后，target_id 变成拼接串
    （alpha-factor_MRFPS_OPN_alpha-pro_OPN_ALPHA_FULL_PROJECT），与筛查按规范靶点存的结果对不上，
    效应列整列消失。且筛查里的 OPN 靶点实际叫 OPN_ALPHA_FULL_PROJECT 而非 OPN——精确匹配也不够。"""
    monkeypatch.setattr(lookup, "available_screen_targets", lambda paths=None: ["hLF", "OPN_ALPHA_FULL_PROJECT"])

    composite = "alpha-factor_MRFPS_OPN_alpha-pro_OPN_ALPHA_FULL_PROJECT"
    assert lookup.resolve_effect_source(composite, "OPN") == "OPN_ALPHA_FULL_PROJECT"
    # 规范靶点原样命中
    assert lookup.resolve_effect_source("hLF", "hLF") == "hLF"
    # 完全无关的上下文不该硬凑一个
    assert lookup.resolve_effect_source("something_else", "KLM") == ""


def test_effect_source_is_reported_so_borrowed_numbers_are_disclosed(monkeypatch) -> None:
    """自定义组合借用规范靶点的效应时，必须把来源报出来——否则会被当成这个构建体自己的预测。"""
    monkeypatch.setattr(
        "app.services.screen_effect_lookup.resolve_effect_source", lambda *a, **k: "hLF"
    )
    monkeypatch.setattr(
        "app.services.screen_effect_lookup.load_screen_effect_lookup",
        lambda *a, **k: {("OE", "R1"): (0.08, 1.0)},
    )
    frame = pd.DataFrame([{"改造方式": "OE", "模型对象": "R1", "作用对象": selector.KIND_COMPLEX}])

    _, has_effects, source = selector._attach_screen_effects(frame, "my_custom_build", target_context="hLF")

    assert has_effects is True
    assert source == "hLF", "来源靶点要回传给界面提示"


def test_selector_shows_dash_when_no_screen_results_exist(monkeypatch) -> None:
    """没跑过筛查是**预期状态**——选择器仍要能用，效应列为空。"""
    monkeypatch.setattr(
        "app.services.screen_effect_lookup.load_screen_effect_lookup", lambda *a, **k: {}
    )
    frame = pd.DataFrame(
        [{"改造方式": "OE", "模型对象": "sec_Pdi1p_complex_formation", "作用对象": selector.KIND_COMPLEX}]
    )

    result, has_effects, _ = selector._attach_screen_effects(frame, "hLF")

    assert has_effects is False
    assert "模型预测提升(%)" not in result.columns


def test_selector_attaches_effect_and_growth_columns(monkeypatch) -> None:
    monkeypatch.setattr(
        "app.services.screen_effect_lookup.load_screen_effect_lookup",
        lambda *a, **k: {("OE", "sec_Pdi1p_complex_formation"): (0.0815, 1.0), ("KO", "PAS_TRAP"): (0.162, 0.1)},
    )
    frame = pd.DataFrame(
        [
            {"改造方式": "OE", "模型对象": "sec_Pdi1p_complex_formation", "作用对象": selector.KIND_COMPLEX},
            {"改造方式": "KO", "模型对象": "PAS_TRAP", "作用对象": selector.KIND_GENE},
            {"改造方式": "OE", "模型对象": "PAS_NOT_SCREENED", "作用对象": selector.KIND_GENE},
        ]
    )

    result, has_effects, _ = selector._attach_screen_effects(frame, "hLF")

    assert has_effects is True
    assert result.iloc[0]["模型预测提升(%)"] == pytest_approx(8.15)
    # 陷阱候选：提升可观但生长只剩 0.1，必须能一眼看见
    assert result.iloc[1]["生长保持"] == pytest_approx(0.1)
    # 没被筛查覆盖的留空（NaN），不能瞎填 0——0 会被读成"测过、没效果"
    assert pd.isna(result.iloc[2]["模型预测提升(%)"])


def test_effect_column_is_numeric_so_sorting_is_by_magnitude_not_alphabet(monkeypatch) -> None:
    """曾是格式化字符串列 → Streamlit 按字典序排，"+8.15%" 会排在 "+10.5%" 之前（'8'>'1'），
    研究员按提升排序时看到的名次是错的。必须存数值。"""
    monkeypatch.setattr(
        "app.services.screen_effect_lookup.load_screen_effect_lookup",
        lambda *a, **k: {("OE", "BIG"): (0.105, 1.0), ("OE", "SMALL"): (0.0815, 1.0)},
    )
    frame = pd.DataFrame(
        [
            {"改造方式": "OE", "模型对象": "SMALL", "作用对象": selector.KIND_GENE},
            {"改造方式": "OE", "模型对象": "BIG", "作用对象": selector.KIND_GENE},
        ]
    )

    result, _, _ = selector._attach_screen_effects(frame, "hLF")

    column = result["模型预测提升(%)"]
    assert pd.api.types.is_numeric_dtype(column), "必须是数值列，否则排序按字典序"
    assert list(column.sort_values(ascending=False).index) == [1, 0], "10.5% 必须排在 8.15% 之前"


def test_tiny_noise_level_effects_stay_numeric_not_scientific_notation(monkeypatch) -> None:
    """真实数据里有 9.27e-05% 这种量级；此前 .3g 直接把科学计数法显示给研究员，读不了。
    现在存数值、由 NumberColumn 定点格式化（显示 0.000）。"""
    monkeypatch.setattr(
        "app.services.screen_effect_lookup.load_screen_effect_lookup",
        lambda *a, **k: {("OE", "NOISE"): (9.27e-07, 1.0)},
    )
    frame = pd.DataFrame([{"改造方式": "OE", "模型对象": "NOISE", "作用对象": selector.KIND_COMPLEX}])

    result, _, _ = selector._attach_screen_effects(frame, "hLF")

    value = result.iloc[0]["模型预测提升(%)"]
    assert isinstance(value, float)
    assert value == pytest_approx(9.27e-05)
