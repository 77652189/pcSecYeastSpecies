"""自建序列库：信号肽 / 引导肽 / 成熟蛋白 / 组合模板的增删改查。

**分层原则**：内置库来自正式科学资产（`Data/`、目标 spec），**只读**——应用不修改受保护资产。
用户自建的条目落 `local_runs/custom_target_library/library.json`（运行工作区），与内置库合并后
供构建页选择。内置条目不可改不可删，用户条目可增可改可删；两者在界面上明确区分来源。

序列校验不是走过场：一条含非法字符的序列会静默产出无意义的仿真结果，所以非法条目**拒绝保存**
并说明原因，而不是"尽力而为"地存进去。
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

from app.ui.common import PATHS


LIBRARY_RELATIVE_PATH = Path("custom_target_library") / "library.json"
LIBRARY_SCHEMA_VERSION = 1

KIND_SIGNAL_PEPTIDE = "signal_peptides"
KIND_LEADER = "leaders"
KIND_MATURE = "mature_proteins"
KIND_TEMPLATE = "templates"
KINDS = (KIND_SIGNAL_PEPTIDE, KIND_LEADER, KIND_MATURE, KIND_TEMPLATE)

KIND_LABELS = {
    KIND_SIGNAL_PEPTIDE: "信号肽",
    KIND_LEADER: "引导肽",
    KIND_MATURE: "成熟蛋白",
    KIND_TEMPLATE: "组合模板",
}

# 20 种标准氨基酸。X/B/Z 这类模糊符号一律不收——模型按残基组成算成本，模糊符号会让结果失真。
AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")
_ID_PATTERN = re.compile(r"^[A-Za-z0-9_\-]+$")


def library_path(paths: Any | None = None) -> Path:
    resolved = paths or PATHS
    return resolved.local_runs_dir / LIBRARY_RELATIVE_PATH


def clean_sequence(raw: str) -> str:
    """去掉空白与换行（粘贴 FASTA 片段很常见），统一大写。"""
    return re.sub(r"\s+", "", str(raw or "")).upper()


def validate_entry(kind: str, entry: dict[str, Any], *, existing_ids: set[str] | None = None) -> list[str]:
    """返回问题列表；空列表＝可保存。"""
    problems: list[str] = []
    entry_id = str(entry.get("id") or "").strip()
    if not entry_id:
        problems.append("缺少编号（id）。")
    elif not _ID_PATTERN.match(entry_id):
        problems.append("编号只能用字母、数字、下划线和连字符。")
    elif existing_ids and entry_id in existing_ids:
        problems.append(f"编号 {entry_id} 已存在（内置条目不可覆盖，自建条目请用“修改”）。")

    if not str(entry.get("label") or "").strip():
        problems.append("缺少名称。")

    sequence_fields = (
        [("sequence", "序列")]
        if kind != KIND_TEMPLATE
        else [("mature_sequence", "成熟蛋白序列")]
    )
    for field, label in sequence_fields:
        sequence = clean_sequence(entry.get(field, ""))
        if not sequence:
            problems.append(f"{label}不能为空。")
            continue
        illegal = sorted(set(sequence) - AMINO_ACIDS)
        if illegal:
            problems.append(f"{label}含非氨基酸字符：{''.join(illegal)}（只接受 20 种标准氨基酸单字母）。")

    if kind == KIND_TEMPLATE:
        for field, label in (("signal_peptide_sequence", "信号肽"), ("leader_sequence", "引导肽")):
            sequence = clean_sequence(entry.get(field, ""))
            if sequence:
                illegal = sorted(set(sequence) - AMINO_ACIDS)
                if illegal:
                    problems.append(f"{label}含非氨基酸字符：{''.join(illegal)}。")

    if kind in {KIND_MATURE, KIND_TEMPLATE}:
        for field, label in (
            ("disulfide_sites", "二硫键数"),
            ("n_glycosylation_sites", "N-糖基化位点数"),
            ("o_glycosylation_sites", "O-糖基化位点数"),
        ):
            value = entry.get(field, 0)
            try:
                if int(value) < 0:
                    problems.append(f"{label}不能为负。")
            except (TypeError, ValueError):
                problems.append(f"{label}必须是整数。")
    return problems


def load_library(paths: Any | None = None) -> dict[str, list[dict[str, Any]]]:
    """读自建库。文件缺失 / 损坏 / schema 不符 → 返回空库（不抛异常）。"""
    empty: dict[str, list[dict[str, Any]]] = {kind: [] for kind in KINDS}
    path = library_path(paths)
    if not path.exists():
        return empty
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError, UnicodeDecodeError):
        return empty
    if not isinstance(payload, dict) or payload.get("schema_version") != LIBRARY_SCHEMA_VERSION:
        return empty
    return {kind: [row for row in payload.get(kind, []) if isinstance(row, dict)] for kind in KINDS}


def _write_library(library: dict[str, list[dict[str, Any]]], paths: Any | None = None) -> None:
    path = library_path(paths)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {"schema_version": LIBRARY_SCHEMA_VERSION, **{kind: library.get(kind, []) for kind in KINDS}}
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")


def save_entry(
    kind: str,
    entry: dict[str, Any],
    *,
    builtin_ids: set[str] | None = None,
    paths: Any | None = None,
) -> tuple[bool, list[str]]:
    """新增或更新一条自建条目。返回 (是否成功, 问题列表)。

    同 id 视为更新；与内置 id 冲突则拒绝——内置来自正式科学资产，不允许被同名条目遮蔽。
    """
    if kind not in KINDS:
        return False, [f"未知类别：{kind}"]
    library = load_library(paths)
    entry_id = str(entry.get("id") or "").strip()
    own_ids = {str(row.get("id") or "") for row in library[kind]}
    # 更新自己时不算重复；与内置冲突永远算重复。
    conflict_ids = set(builtin_ids or set()) | (own_ids - {entry_id})
    problems = validate_entry(kind, entry, existing_ids=conflict_ids)
    if problems:
        return False, problems

    normalized = dict(entry)
    normalized["id"] = entry_id
    for field in ("sequence", "signal_peptide_sequence", "leader_sequence", "mature_sequence"):
        if field in normalized:
            normalized[field] = clean_sequence(normalized[field])
    for field in ("disulfide_sites", "n_glycosylation_sites", "o_glycosylation_sites"):
        if field in normalized:
            normalized[field] = int(normalized[field] or 0)

    library[kind] = [row for row in library[kind] if str(row.get("id") or "") != entry_id] + [normalized]
    _write_library(library, paths)
    return True, []


def delete_entry(kind: str, entry_id: str, *, paths: Any | None = None) -> bool:
    """删除自建条目；内置条目不在此库中，天然删不掉。"""
    if kind not in KINDS:
        return False
    library = load_library(paths)
    remaining = [row for row in library[kind] if str(row.get("id") or "") != str(entry_id)]
    if len(remaining) == len(library[kind]):
        return False
    library[kind] = remaining
    _write_library(library, paths)
    return True


def merge_with_builtin(
    builtin: dict[str, dict[str, Any]],
    kind: str,
    *,
    paths: Any | None = None,
) -> dict[str, dict[str, Any]]:
    """内置库 + 自建库合并成一个下拉可选集合，并标注来源（内置只读）。"""
    merged: dict[str, dict[str, Any]] = {}
    for key, value in builtin.items():
        merged[key] = {**value, "source": "builtin", "editable": False}
    for row in load_library(paths).get(kind, []):
        entry_id = str(row.get("id") or "")
        if not entry_id:
            continue
        merged[entry_id] = {**row, "source": "custom", "editable": True}
    return merged


__all__ = [
    "AMINO_ACIDS",
    "KINDS",
    "KIND_LABELS",
    "KIND_LEADER",
    "KIND_MATURE",
    "KIND_SIGNAL_PEPTIDE",
    "KIND_TEMPLATE",
    "LIBRARY_SCHEMA_VERSION",
    "clean_sequence",
    "delete_entry",
    "library_path",
    "load_library",
    "merge_with_builtin",
    "save_entry",
    "validate_entry",
]
