"""E3 基因↔复合体映射的服务门面（ADR-007）。

职责只有三件：定位策展文件、读进来（缺席即优雅降级为空）、把"实验时该动哪几个基因"整理成
UI 能直接显示的一句话。判定与门禁（契约校验、`review_status` / `stoichiometry_status`）全在引擎
`pcsec_pichia.services.gene_complex_mapping`，本模块不做科学判断。

策展文件位置：默认 `Data/pcSecPichia/gene_complex_mapping.json`（人工策展的稳定资产，由策展者放置）。
本模块**只读**——受保护科学资产目录一律不写（写入一律落 `local_runs/`）。
"""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any

from app import ensure_python_pichia_on_path
from app.ui.common import PATHS


CURATED_MAPPING_RELATIVE_PATH = Path("Data") / "pcSecPichia" / "gene_complex_mapping.json"
# 应用内复核结果的工作副本。沿用项目既有模式（实验反馈也是先落 local_runs、人工确认后才提升到
# Data/）：应用**可以**写工作区，只是不自动写受保护的科学资产目录。研究员点确认即存这里、立即生效，
# 不需要"导出再导入"。要不要提升为正式资产是另一件事（人显式提交）。
WORKING_MAPPING_RELATIVE_PATH = Path("gene_complex_mapping") / "reviewed.json"


def curated_gene_complex_mapping_path(paths: Any | None = None) -> Path:
    resolved = paths or PATHS
    return resolved.repo_root / CURATED_MAPPING_RELATIVE_PATH


def working_gene_complex_mapping_path(paths: Any | None = None) -> Path:
    resolved = paths or PATHS
    return resolved.local_runs_dir / WORKING_MAPPING_RELATIVE_PATH


def save_reviewed_mappings(rows: list[dict[str, Any]], *, paths: Any | None = None) -> tuple[int, tuple[str, ...]]:
    """把界面上的复核结果存进工作副本，立即生效。返回 (写入条数, 被契约拒掉的问题)。

    先过引擎的契约校验——不合规的条目丢弃并报告，绝不落盘成"看起来像证据"的垃圾。
    """
    ensure_python_pichia_on_path()
    from pcsec_pichia.services.gene_complex_mapping import (
        serialize_gene_complex_mappings,
        validate_gene_complex_mapping_payloads,
    )

    valid, problems = validate_gene_complex_mapping_payloads(rows)
    payload = serialize_gene_complex_mappings(valid)
    output_path = working_gene_complex_mapping_path(paths)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    _load_cached.cache_clear()  # 让下一次读取立刻看到新结果
    return len(valid), problems


@lru_cache(maxsize=4)
def _load_cached(path_text: str, mtime: float) -> tuple[tuple[Any, ...], tuple[str, ...]]:
    ensure_python_pichia_on_path()
    from pcsec_pichia.services.gene_complex_mapping import load_gene_complex_mapping_file

    return load_gene_complex_mapping_file(Path(path_text))


def load_gene_complex_mapping(paths: Any | None = None) -> tuple[tuple[Any, ...], tuple[str, ...]]:
    """读映射：正式资产（`Data/`）+ 应用内复核的工作副本（`local_runs/`）。

    同一 (基因, 复合体) 若两处都有，**正式资产优先**——它经过显式提交，工作副本只是尚未提升的复核结果。
    两处都没有时返回 ((), (说明,))，全链路优雅降级。
    """
    merged: list[Any] = []
    notes: list[str] = []
    seen: set[tuple[str, str]] = set()
    for path in (curated_gene_complex_mapping_path(paths), working_gene_complex_mapping_path(paths)):
        # mtime 入 key：文件更新后自动失效；文件不存在时用 0.0，避免 stat 抛错。
        mtime = path.stat().st_mtime if path.exists() else 0.0
        rows, file_notes = _load_cached(str(path), mtime)
        for row in rows:
            key = (row.pichia_gene_id, row.complex_reaction_id)
            if key in seen:
                continue
            seen.add(key)
            merged.append(row)
        notes.extend(file_notes)
    return tuple(merged), tuple(notes)


def gene_complex_mapping_available(paths: Any | None = None) -> bool:
    rows, _ = load_gene_complex_mapping(paths)
    return bool(rows)


def lab_genes_for_complex(
    complex_reaction_id: str,
    *,
    paths: Any | None = None,
    include_pending_review: bool = True,
) -> list[dict[str, str]]:
    """模型说改这个复合体有效 → 实验室该动哪几个基因。无策展数据时返回空列表。"""
    rows, _ = load_gene_complex_mapping(paths)
    if not rows:
        return []
    ensure_python_pichia_on_path()
    from pcsec_pichia.services.gene_complex_mapping import genes_for_complex

    matches = genes_for_complex(rows, complex_reaction_id, reviewed_only=not include_pending_review)
    return [
        {
            "gene_id": row.pichia_gene_id,
            "subunit_role": row.subunit_role,
            "review_status": row.review_status,
            "stoichiometry_status": row.stoichiometry_status,
            "evidence_source": row.evidence_source,
            "single_gene_oe_eligible": "yes" if row.may_enter_executable_single_gene_oe else "no",
        }
        for row in matches
    ]


def lab_gene_hint_for_complex(complex_reaction_id: str, *, paths: Any | None = None) -> str:
    """给 UI 的一句话："实验时对应基因：…"。无策展数据 → 空串（UI 应显示为"—"）。"""
    matches = lab_genes_for_complex(complex_reaction_id, paths=paths)
    if not matches:
        return ""
    parts: list[str] = []
    for match in matches:
        suffix = "" if match["review_status"] == "reviewed" else "（待复核）"
        parts.append(f"{match['gene_id']}{suffix}")
    return "、".join(parts)


def build_draft_mapping_rows(paths: Any | None = None) -> list[dict[str, Any]]:
    """从策展候选自动起草映射，返回纯 dict 供 UI 渲染（UI 不得直接 import 引擎）。

    草稿一律 pending_review + auxiliary + unknown，**不会自己生效**——判定逻辑在引擎，这里只搬运。
    """
    ensure_python_pichia_on_path()
    from pcsec_pichia.services.gene_complex_mapping import build_draft_mappings_from_candidates

    from app.services.pichia_gene_catalog_service import load_hlf_opn_candidate_genes

    candidates = load_hlf_opn_candidate_genes(target_context=None, include_shared=True)
    return [row.to_dict() for row in build_draft_mappings_from_candidates(candidates)]


def gene_complex_mapping_summary(paths: Any | None = None) -> dict[str, Any]:
    rows, notes = load_gene_complex_mapping(paths)
    ensure_python_pichia_on_path()
    from pcsec_pichia.services.gene_complex_mapping import summarize_gene_complex_mapping_rows

    summary = dict(summarize_gene_complex_mapping_rows(rows))
    summary["notes"] = list(notes)
    summary["curated_file"] = str(CURATED_MAPPING_RELATIVE_PATH).replace("\\", "/")
    return summary


__all__ = [
    "CURATED_MAPPING_RELATIVE_PATH",
    "WORKING_MAPPING_RELATIVE_PATH",
    "build_draft_mapping_rows",
    "curated_gene_complex_mapping_path",
    "save_reviewed_mappings",
    "working_gene_complex_mapping_path",
    "gene_complex_mapping_available",
    "gene_complex_mapping_summary",
    "lab_gene_hint_for_complex",
    "lab_genes_for_complex",
    "load_gene_complex_mapping",
]
