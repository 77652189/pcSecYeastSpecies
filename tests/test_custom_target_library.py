"""自建序列库的增删改查 + 校验。

关键：一条含非法字符的序列会静默产出无意义的仿真结果，所以非法条目**拒绝保存并说明原因**，
而不是"尽力而为"地存进去。内置条目来自正式科学资产，不可被覆盖/删除。
"""

from __future__ import annotations

import json
from typing import Any

import pytest

from app.services import custom_target_library_service as lib


@pytest.fixture()
def paths(tmp_path):
    class _Paths:
        repo_root = tmp_path
        local_runs_dir = tmp_path / "local_runs"

    return _Paths


def _entry(**overrides: Any) -> dict[str, Any]:
    row: dict[str, Any] = {"id": "my_sp", "label": "我的信号肽", "sequence": "MKLVFLVLLFLGALG"}
    row.update(overrides)
    return row


def test_create_read_update_delete_round_trip(paths) -> None:
    ok, problems = lib.save_entry(lib.KIND_SIGNAL_PEPTIDE, _entry(), paths=paths)
    assert (ok, problems) == (True, [])
    assert [row["id"] for row in lib.load_library(paths)[lib.KIND_SIGNAL_PEPTIDE]] == ["my_sp"]

    ok, _ = lib.save_entry(lib.KIND_SIGNAL_PEPTIDE, _entry(label="改过的名字"), paths=paths)
    stored = lib.load_library(paths)[lib.KIND_SIGNAL_PEPTIDE]
    assert ok and len(stored) == 1, "同 id 应该是更新而不是追加"
    assert stored[0]["label"] == "改过的名字"

    assert lib.delete_entry(lib.KIND_SIGNAL_PEPTIDE, "my_sp", paths=paths) is True
    assert lib.load_library(paths)[lib.KIND_SIGNAL_PEPTIDE] == []
    assert lib.delete_entry(lib.KIND_SIGNAL_PEPTIDE, "my_sp", paths=paths) is False


def test_illegal_amino_acids_are_rejected_not_stored(paths) -> None:
    ok, problems = lib.save_entry(lib.KIND_SIGNAL_PEPTIDE, _entry(sequence="MKLV123XZ"), paths=paths)

    assert ok is False
    assert any("非氨基酸字符" in problem for problem in problems)
    assert lib.load_library(paths)[lib.KIND_SIGNAL_PEPTIDE] == [], "非法条目绝不能落盘"


def test_sequence_whitespace_and_case_are_normalised(paths) -> None:
    """粘贴 FASTA 片段带换行空格很常见，不该因此报错。"""
    lib.save_entry(lib.KIND_SIGNAL_PEPTIDE, _entry(sequence=" mklv flvl\nlflg alg "), paths=paths)

    assert lib.load_library(paths)[lib.KIND_SIGNAL_PEPTIDE][0]["sequence"] == "MKLVFLVLLFLGALG"


def test_builtin_ids_cannot_be_shadowed(paths) -> None:
    ok, problems = lib.save_entry(
        lib.KIND_SIGNAL_PEPTIDE, _entry(id="alpha_factor"), builtin_ids={"alpha_factor"}, paths=paths
    )

    assert ok is False
    assert any("已存在" in problem for problem in problems)


def test_missing_id_label_or_sequence_is_reported(paths) -> None:
    ok, problems = lib.save_entry(lib.KIND_SIGNAL_PEPTIDE, {"id": "", "label": "", "sequence": ""}, paths=paths)

    assert ok is False
    assert len(problems) >= 3


def test_id_charset_is_restricted(paths) -> None:
    ok, problems = lib.save_entry(lib.KIND_SIGNAL_PEPTIDE, _entry(id="有 空格/斜杠"), paths=paths)

    assert ok is False
    assert any("编号" in problem for problem in problems)


def test_template_requires_mature_sequence_and_validates_optional_segments(paths) -> None:
    ok, problems = lib.save_entry(
        lib.KIND_TEMPLATE,
        {"id": "t1", "label": "模板", "mature_sequence": "", "signal_peptide_sequence": "MKL"},
        paths=paths,
    )
    assert ok is False and any("成熟蛋白序列" in problem for problem in problems)

    ok, problems = lib.save_entry(
        lib.KIND_TEMPLATE,
        {"id": "t1", "label": "模板", "mature_sequence": "MKLV", "leader_sequence": "MK123"},
        paths=paths,
    )
    assert ok is False and any("引导肽" in problem for problem in problems)


def test_negative_ptm_counts_are_rejected(paths) -> None:
    ok, problems = lib.save_entry(
        lib.KIND_MATURE,
        {"id": "m1", "label": "蛋白", "sequence": "MKLV", "disulfide_sites": -1},
        paths=paths,
    )

    assert ok is False and any("不能为负" in problem for problem in problems)


def test_merge_marks_builtin_readonly_and_custom_editable(paths) -> None:
    lib.save_entry(lib.KIND_SIGNAL_PEPTIDE, _entry(), paths=paths)

    merged = lib.merge_with_builtin({"builtin_sp": {"label": "内置", "sequence": "MK"}}, lib.KIND_SIGNAL_PEPTIDE, paths=paths)

    assert merged["builtin_sp"]["editable"] is False
    assert merged["builtin_sp"]["source"] == "builtin"
    assert merged["my_sp"]["editable"] is True
    assert merged["my_sp"]["source"] == "custom"


def test_copying_a_builtin_entry_under_a_new_id_is_how_editing_starts(paths) -> None:
    """内置只读 ⇒ 没有自建条目时"修改"无从下手。起步方式是"复制内置 → 改 → 存新编号"，
    这条路径必须走得通（用户 2026-07-28 反馈"看上去只有新建"就是因为它缺失）。"""
    builtin = {"alpha_factor": {"label": "内置 alpha-factor", "sequence": "MRFPSIFTAVLFAASSALA"}}

    copied = dict(lib.merge_with_builtin(builtin, lib.KIND_SIGNAL_PEPTIDE, paths=paths)["alpha_factor"])
    copied["id"] = "alpha_factor_my_variant"
    copied["label"] = "我的改版"
    copied["sequence"] = copied["sequence"] + "GG"

    ok, problems = lib.save_entry(
        lib.KIND_SIGNAL_PEPTIDE, copied, builtin_ids=set(builtin), paths=paths
    )

    assert (ok, problems) == (True, [])
    merged = lib.merge_with_builtin(builtin, lib.KIND_SIGNAL_PEPTIDE, paths=paths)
    assert merged["alpha_factor"]["sequence"] == "MRFPSIFTAVLFAASSALA", "内置条目必须原样不动"
    assert merged["alpha_factor_my_variant"]["editable"] is True
    assert merged["alpha_factor_my_variant"]["sequence"].endswith("GG")


def test_corrupt_or_stale_library_degrades_to_empty(paths) -> None:
    path = lib.library_path(paths)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("{not json", encoding="utf-8")
    assert lib.load_library(paths)[lib.KIND_SIGNAL_PEPTIDE] == []

    path.write_text(json.dumps({"schema_version": 999, "signal_peptides": [_entry()]}), encoding="utf-8")
    assert lib.load_library(paths)[lib.KIND_SIGNAL_PEPTIDE] == []
