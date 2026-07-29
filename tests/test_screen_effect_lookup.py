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


def test_selector_shows_dash_when_no_screen_results_exist(monkeypatch) -> None:
    """没跑过筛查是**预期状态**——选择器仍要能用，效应列为空。"""
    monkeypatch.setattr(
        "app.services.screen_effect_lookup.load_screen_effect_lookup", lambda *a, **k: {}
    )
    frame = pd.DataFrame(
        [{"改造方式": "OE", "模型对象": "sec_Pdi1p_complex_formation", "作用对象": selector.KIND_COMPLEX}]
    )

    result, has_effects = selector._attach_screen_effects(frame, "hLF")

    assert has_effects is False
    assert "模型预测提升" not in result.columns


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

    result, has_effects = selector._attach_screen_effects(frame, "hLF")

    assert has_effects is True
    assert result.iloc[0]["模型预测提升"].startswith("+8.15")
    # 陷阱候选：提升可观但生长只剩 0.1，必须能一眼看见
    assert result.iloc[1]["生长保持"] == "0.100"
    assert result.iloc[2]["模型预测提升"] == "—", "没被筛查覆盖的显示为 —，不能瞎填 0"
