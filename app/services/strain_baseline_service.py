"""改造后候选系统的"野生型基线读服务"（D4 复用地基 · 服务层 · 零新增求解）。

D2/D3 的分层复用要拿一份**同 run 口径**的野生型基线做锚点。本服务把 run 的口径
(目标·碳源·培养基·μ·折叠/翻译约束档·兼容口径) 映成 [strain_baseline_cache] 的口径指纹，
读缓存里蒸馏好的候选基线（OE+KO 逐候选层+相对效应）；命中即返回，未命中就诚实报"该口径下
还没有后台基线"，让面板给出触发后台构建的指引——本服务**从不求解、从不伪造**。

分层：UI（app/ui）禁止 import 引擎，故口径指纹 / 缓存读写这类引擎调用集中在本服务层，返回纯 dict。
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any, Sequence

from app import ensure_python_pichia_on_path

ensure_python_pichia_on_path()

from pcsec_pichia.strain_baseline_cache import (  # noqa: E402 - engine import after path bootstrap
    DEFAULT_STRAIN_BASELINE_CACHE_DIR,
    StrainBaselineCacheKey,
    cache_key_digest,
    distill_tradeoff_rows,
    load_cached_baseline,
    store_baseline,
    strain_fingerprint,
)
from pcsec_pichia.strain_modifications import StrainModifications  # noqa: E402


def build_baseline_cache_key(
    *,
    target_id: str,
    carbon_source_id: str,
    media_type: int,
    growth_rate: float,
    enable_ribosome_translation_constraint: bool,
    enable_misfolding_constraint: bool,
    compatibility_mode: str = "corrected",
    mode: str = "fast",
    solver_method: str = "default",
    modifications: StrainModifications | None = None,
) -> StrainBaselineCacheKey:
    """把 run 口径（+ 可选改造）映成基线缓存键。野生型（无改造）→ model_variant_fingerprint="wildtype"。"""
    return StrainBaselineCacheKey(
        target_id=str(target_id),
        carbon_source_id=str(carbon_source_id),
        media_type=int(media_type),
        growth_rate=float(growth_rate),
        write_ribosome_translation_constraint=bool(enable_ribosome_translation_constraint),
        write_misfolding_constraints=bool(enable_misfolding_constraint),
        mode=str(mode),
        compatibility_mode=str(compatibility_mode),
        solver_method=str(solver_method),
        model_variant_fingerprint=strain_fingerprint(modifications),
    )


def _partition(rows: Sequence[dict[str, Any]]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """按 intervention_type 分成 (OE, KO) 两组（不排序、不打标——排序/复用标注是 D3 的职责）。"""
    oe = [row for row in rows if str(row.get("intervention_type")).upper() == "OE"]
    ko = [row for row in rows if str(row.get("intervention_type")).upper() == "KO"]
    return oe, ko


def load_strain_baseline_readout(
    *,
    target_id: str,
    carbon_source_id: str,
    media_type: int,
    growth_rate: float,
    enable_ribosome_translation_constraint: bool,
    enable_misfolding_constraint: bool,
    compatibility_mode: str = "corrected",
    mode: str = "fast",
    solver_method: str = "default",
    modifications: StrainModifications | None = None,
    cache_dir: Path | str = DEFAULT_STRAIN_BASELINE_CACHE_DIR,
) -> dict[str, Any]:
    """读该口径的野生型（或指定改造）基线；命中返回候选行，未命中诚实报缺（供面板给触发指引）。"""
    key = build_baseline_cache_key(
        target_id=target_id,
        carbon_source_id=carbon_source_id,
        media_type=media_type,
        growth_rate=growth_rate,
        enable_ribosome_translation_constraint=enable_ribosome_translation_constraint,
        enable_misfolding_constraint=enable_misfolding_constraint,
        compatibility_mode=compatibility_mode,
        mode=mode,
        solver_method=solver_method,
        modifications=modifications,
    )
    caliber = {
        "target_id": key.target_id,
        "carbon_source_id": key.carbon_source_id,
        "media_type": key.media_type,
        "growth_rate": key.growth_rate,
        "write_ribosome_translation_constraint": key.write_ribosome_translation_constraint,
        "write_misfolding_constraints": key.write_misfolding_constraints,
        "mode": key.mode,
        "compatibility_mode": key.compatibility_mode,
    }
    payload = load_cached_baseline(key, cache_dir)
    if payload is None:
        return {
            "available": False,
            "caliber": caliber,
            "model_variant_fingerprint": key.model_variant_fingerprint,
            "cache_key_digest": cache_key_digest(key),
            "reason": "该口径下还没有后台基线缓存——请对此口径跑野生型全基因组后台构建后再读。",
            "oe_rows": [],
            "ko_rows": [],
            "candidate_count": 0,
        }
    rows = list(payload.get("rows") or [])
    oe_rows, ko_rows = _partition(rows)
    return {
        "available": True,
        "caliber": caliber,
        "model_variant_fingerprint": key.model_variant_fingerprint,
        "cache_key_digest": cache_key_digest(key),
        "built_at": payload.get("built_at"),
        "source_run": payload.get("source_run"),
        "candidate_count": payload.get("candidate_count", len(rows)),
        "rows": rows,
        "oe_rows": oe_rows,
        "ko_rows": ko_rows,
    }


def ingest_tradeoff_csv_into_baseline_cache(
    *,
    csv_path: Path | str,
    target_id: str,
    carbon_source_id: str,
    media_type: int,
    growth_rate: float,
    enable_ribosome_translation_constraint: bool,
    enable_misfolding_constraint: bool,
    compatibility_mode: str = "corrected",
    mode: str = "fast",
    solver_method: str = "default",
    modifications: StrainModifications | None = None,
    source_run: str | None = None,
    cache_dir: Path | str = DEFAULT_STRAIN_BASELINE_CACHE_DIR,
) -> dict[str, Any]:
    """把一次已完成的全基因组筛查 CSV 蒸馏进该口径的基线缓存（后台构建的收尾一步）。

    调用方（后台构建工具）负责保证 CSV 确实是在这份口径下跑出来的——本函数按 target 蒸馏 KO/OE 行、
    按口径指纹落缓存，供面板/ D2/D3 复用。
    """
    key = build_baseline_cache_key(
        target_id=target_id,
        carbon_source_id=carbon_source_id,
        media_type=media_type,
        growth_rate=growth_rate,
        enable_ribosome_translation_constraint=enable_ribosome_translation_constraint,
        enable_misfolding_constraint=enable_misfolding_constraint,
        compatibility_mode=compatibility_mode,
        mode=mode,
        solver_method=solver_method,
        modifications=modifications,
    )
    with Path(csv_path).open("r", newline="", encoding="utf-8") as handle:
        records = list(csv.DictReader(handle))
    rows = distill_tradeoff_rows(records, target_id=target_id)
    stored_path = store_baseline(key, rows, cache_dir, source_run=source_run)
    return {
        "stored_path": str(stored_path),
        "candidate_count": len(rows),
        "cache_key_digest": cache_key_digest(key),
        "model_variant_fingerprint": key.model_variant_fingerprint,
    }


__all__ = [
    "build_baseline_cache_key",
    "ingest_tradeoff_csv_into_baseline_cache",
    "load_strain_baseline_readout",
]
