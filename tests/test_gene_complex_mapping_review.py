"""策展复核面板：界面里改 → 二次确认 → 直接保存生效（不做导出再导入）。

用户 2026-07-28 判定"导出再导入毫无必要"。约束是**不自动写受保护的 Data/**，不是不能保存——
沿用项目既有模式：复核结果落 local_runs/ 工作区并立即生效，提升为正式资产是另一个人工步骤。
"""

from __future__ import annotations

import json
from typing import Any

import pandas as pd

from app.services import gene_complex_mapping_service as service
from app.ui.views import gene_complex_mapping_review as review


def _edited(**overrides: Any) -> pd.DataFrame:
    row: dict[str, Any] = {
        "模型预测提升(%)": 8.152,
        "复合体反应": "sec_Pdi1p_complex_formation",
        "候选基因": "PAS_chr4_0844",
        "来源俗名": "PDI1（单独）",
        "复核结论": "确认参与",
        "亚基角色": "必需亚基",
        "化学计量": "已知",
        "判断依据": "文献 X",
    }
    row.update(overrides)
    return pd.DataFrame([row])


def test_chinese_choices_map_onto_the_contract_vocabulary() -> None:
    payload = review.rows_to_contract_payload(_edited())

    assert payload == [
        {
            "pichia_gene_id": "PAS_chr4_0844",
            "complex_reaction_id": "sec_Pdi1p_complex_formation",
            "subunit_role": "required_subunit",
            "stoichiometry_status": "known",
            "review_status": "reviewed",
            "evidence_source": "curated_review",
            "evidence_citation": "策展俗名 PDI1（单独）",
            "note": "文献 X",
        }
    ]


def test_untouched_rows_are_not_written() -> None:
    """没判断过的东西不该落盘占位。"""
    assert review.rows_to_contract_payload(_edited(复核结论="待复核")) == []


def test_rejected_rows_are_written_so_the_judgement_is_not_lost() -> None:
    """否决也是判断——要记下来，否则下次又冒出来让人重审。"""
    payload = review.rows_to_contract_payload(_edited(复核结论="否决（同源猜错）"))

    assert len(payload) == 1
    assert payload[0]["review_status"] == "rejected"


def test_save_writes_working_copy_and_takes_effect_immediately(tmp_path, monkeypatch) -> None:
    class _Paths:
        repo_root = tmp_path
        local_runs_dir = tmp_path / "local_runs"

    monkeypatch.setattr(service, "PATHS", _Paths)
    service._load_cached.cache_clear()

    saved, problems = service.save_reviewed_mappings(
        review.rows_to_contract_payload(_edited()), paths=_Paths
    )

    assert saved == 1
    assert problems == ()
    written = service.working_gene_complex_mapping_path(_Paths)
    assert written.exists(), "必须真的落盘，不能只是给个下载按钮"
    payload = json.loads(written.read_text(encoding="utf-8"))
    assert payload["mappings"][0]["pichia_gene_id"] == "PAS_chr4_0844"

    # 立即生效：保存后再读就能查到（缓存已失效）
    rows, _ = service.load_gene_complex_mapping(_Paths)
    assert [row.pichia_gene_id for row in rows] == ["PAS_chr4_0844"]
    assert rows[0].may_enter_executable_single_gene_oe is True


def test_contract_violations_are_reported_and_not_persisted(tmp_path, monkeypatch) -> None:
    class _Paths:
        repo_root = tmp_path
        local_runs_dir = tmp_path / "local_runs"

    monkeypatch.setattr(service, "PATHS", _Paths)
    service._load_cached.cache_clear()

    saved, problems = service.save_reviewed_mappings(
        [{"pichia_gene_id": "", "complex_reaction_id": "R", "subunit_role": "auxiliary",
          "stoichiometry_status": "unknown", "review_status": "reviewed", "evidence_source": "x"}],
        paths=_Paths,
    )

    assert saved == 0
    assert problems, "不合规的条目要报出来，不能静默丢弃"


def test_formal_asset_wins_over_working_copy(tmp_path, monkeypatch) -> None:
    """正式资产经过显式提交，优先于尚未提升的工作副本。"""
    class _Paths:
        repo_root = tmp_path
        local_runs_dir = tmp_path / "local_runs"

    monkeypatch.setattr(service, "PATHS", _Paths)
    service._load_cached.cache_clear()

    def _entry(note: str) -> dict[str, Any]:
        return {
            "pichia_gene_id": "PAS_A", "complex_reaction_id": "sec_X_complex_formation",
            "subunit_role": "required_subunit", "stoichiometry_status": "known",
            "review_status": "reviewed", "evidence_source": "s", "note": note,
        }

    working = service.working_gene_complex_mapping_path(_Paths)
    working.parent.mkdir(parents=True, exist_ok=True)
    working.write_text(json.dumps({"schema_version": 1, "mappings": [_entry("工作副本")]}), encoding="utf-8")
    formal = service.curated_gene_complex_mapping_path(_Paths)
    formal.parent.mkdir(parents=True, exist_ok=True)
    formal.write_text(json.dumps({"schema_version": 1, "mappings": [_entry("正式资产")]}), encoding="utf-8")

    rows, _ = service.load_gene_complex_mapping(_Paths)

    assert len(rows) == 1, "同一 (基因, 复合体) 不得重复计入"
    assert rows[0].note == "正式资产"
