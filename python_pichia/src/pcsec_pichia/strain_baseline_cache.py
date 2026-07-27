"""改造后候选系统的"口径指纹基线缓存"（D4 复用地基 · ADR-004 #1 迭代2）。

分层复用要在改造后菌株上**只重算受影响层、其余复用野生型结果**（D2/D3）。复用的前提是有一份
**同口径**的野生型基线可比：口径 = (目标 · 碳源 · 培养基 · μ · 折叠/翻译约束档 · 兼容口径 · solver)。
本模块按这份口径的内容哈希做键，把全基因组筛查蒸馏出的候选基线（OE+KO 逐候选：分泌层 + 相对效应
+ 生长保持）存本地缓存；D2/D3 的打标/复用直接读缓存、**零新增求解**。

设计边界（与 [solve_cache] 同源）：
- 缓存是本地派生产物（默认 ``local_runs/strain_baseline_cache/``，随 local_runs gitignored），
  可随时删重建；不入 git（也因此不会外泄任何改造规格）。
- 键必须涵盖决定基线的全部口径输入；分类口径 / 筛查语义变更必须 bump ``BASELINE_CACHE_SCHEMA_VERSION``。
  **本版本 post-date `87f99ac` secretory_process 分类修复**——旧的 `unknown` 分类基线键即失配、
  永不误命中，D2/D3 只吃新分类（"不碰旧 unknown 缓存"的硬保证）。
- 菌株指纹 ``model_variant_fingerprint``：野生型 = ``"wildtype"``（分层复用的可比锚点）；已改造 =
  KO/OE 规格的稳定哈希（复用 solve_cache 同一 `model_variant_fingerprint` 思路）。
- 复用是**近似**（LP 全局耦合），本模块只负责"存/取野生型基线"这一可比锚点；改造后哪些候选可复用
  由 D2 判定并**显式标注失效边界**，本模块不下结论、不做求解。
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Iterable, Sequence

from pcsec_pichia.strain_modifications import StrainModifications

# 分类口径 / 筛查语义变更时必须 bump。post-date 87f99ac：旧 unknown 分类基线永不命中。
BASELINE_CACHE_SCHEMA_VERSION = "2026-07-24-post-classification-fix"
DEFAULT_STRAIN_BASELINE_CACHE_DIR = Path("local_runs") / "strain_baseline_cache"

# 蒸馏进基线的逐候选字段（D2/D3 打标 / 复用排序只需这些；与 build_shortlist_readout 消费口径一致）。
# `secretory_process` 是筛查侧的中文展示标签（gene_perturbation_map 词表）；`affected_reactions` 是该候选
# 扰动的反应，D2 用它经 classify_secretory_process **重新分类**成与 LP 瓶颈同词表的层，做受影响层比较——
# 两套分类词表不同（展示标签 vs 归因键），不能直接字符串比对，故基线必须保留反应级信息。
BASELINE_ROW_FIELDS: tuple[str, ...] = (
    "target_id",
    "gene_id",
    "common_name",
    "candidate_kind",
    "intervention_type",  # KO / OE
    "secretory_process",  # 分泌层的中文展示标签（仅用于展示，不用于跨源层比较）
    "affected_reactions",  # 该候选扰动的反应（; 连接）；D2 据此重分类出可比层
    "secretion_ratio_vs_wildtype",
    "growth_retention_ratio",
    "mapping_confidence",
    "support_status",
)
_BASELINE_NUMERIC_FIELDS = ("secretion_ratio_vs_wildtype", "growth_retention_ratio")


@dataclass(frozen=True)
class StrainBaselineCacheKey:
    """决定一份候选基线的全部口径输入（内容寻址键）。

    ``growth_rate`` 是基线构建时的 μ（并行筛查工具的 ``reference_growth_rate``）；``mode`` 是 μ 扫描
    档（fast/precise，决定 max_feasible_mu 粒度）。面板"跟随 run 口径"读缓存时，用 run 的
    (目标·碳源·培养基·μ·约束档·兼容口径) 直接算出同一个键。
    """

    target_id: str
    carbon_source_id: str
    media_type: int
    growth_rate: float
    write_ribosome_translation_constraint: bool
    write_misfolding_constraints: bool
    mode: str = "fast"
    compatibility_mode: str = "corrected"
    solver_method: str = "default"
    # 菌株变体指纹：野生型 = "wildtype"。改造后必须相应设置，否则会误命中野生型基线。
    model_variant_fingerprint: str = "wildtype"
    entrypoint: str = "strain_baseline_shortlist"
    schema_version: str = BASELINE_CACHE_SCHEMA_VERSION


def _canonical_json(obj: Any) -> str:
    return json.dumps(obj, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def strain_fingerprint(modifications: StrainModifications | None) -> str:
    """菌株改造规格的稳定指纹（复用 solve_cache ``model_variant_fingerprint`` 思路）。

    空 / None → ``"wildtype"``（野生型基线，分层复用的可比锚点）。非空 → KO/OE 规格的 sha256 短哈希。
    KO/OE id 先排序再哈希：改造是集合语义，先后顺序不该改变指纹。
    """
    if modifications is None or modifications.is_empty():
        return "wildtype"
    payload = {
        "ko": sorted(str(r).strip() for r in modifications.ko_reaction_ids if str(r).strip()),
        "oe": sorted(str(r).strip() for r in modifications.oe_reaction_ids if str(r).strip()),
        "oe_factor": round(float(modifications.oe_factor), 6),
    }
    digest = hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()
    return f"strain-{digest[:16]}"


def cache_key_digest(key: StrainBaselineCacheKey) -> str:
    """键的稳定 sha256 摘要（内容寻址）。"""
    return hashlib.sha256(_canonical_json(asdict(key)).encode("utf-8")).hexdigest()


def cache_path(
    key: StrainBaselineCacheKey,
    cache_dir: Path | str = DEFAULT_STRAIN_BASELINE_CACHE_DIR,
) -> Path:
    return Path(cache_dir) / f"{key.entrypoint}_{cache_key_digest(key)}.json"


def _coerce_row(record: dict[str, Any]) -> dict[str, Any]:
    """把一条筛查记录裁到 BASELINE_ROW_FIELDS，数值字段安全转 float（不可解析 → None）。"""
    row: dict[str, Any] = {}
    for field_name in BASELINE_ROW_FIELDS:
        value = record.get(field_name)
        if field_name in _BASELINE_NUMERIC_FIELDS:
            try:
                row[field_name] = float(value) if value not in (None, "") else None
            except (TypeError, ValueError):
                row[field_name] = None
        else:
            row[field_name] = "" if value is None else str(value)
    return row


def distill_tradeoff_rows(
    records: Iterable[dict[str, Any]],
    target_id: str,
) -> list[dict[str, Any]]:
    """把全基因组筛查的 tradeoff 记录蒸馏成基线候选行（纯函数，无 pandas 依赖）。

    只保留该 target、且 intervention_type 为 KO/OE 的记录，裁到 BASELINE_ROW_FIELDS。CSV→dict 记录
    的转换（含 pandas）留在调用方（工具 / 服务层），本引擎函数只吃朴素 dict、方便单测。
    """
    rows: list[dict[str, Any]] = []
    for record in records:
        if str(record.get("target_id")) != str(target_id):
            continue
        if str(record.get("intervention_type")).upper() not in ("KO", "OE"):
            continue
        rows.append(_coerce_row(record))
    return rows


def baseline_to_payload(
    key: StrainBaselineCacheKey,
    rows: Sequence[dict[str, Any]],
    *,
    source_run: str | None = None,
    built_at: str | None = None,
    extra: dict[str, Any] | None = None,
) -> dict[str, Any]:
    payload = {
        "schema_version": key.schema_version,
        "cache_key": asdict(key),
        "cache_key_digest": cache_key_digest(key),
        "built_at": built_at or datetime.now(timezone.utc).isoformat(),
        "source_run": source_run,
        "candidate_count": len(rows),
        "rows": [dict(row) for row in rows],
    }
    if extra:
        payload.update(extra)
    return payload


def load_cached_baseline(
    key: StrainBaselineCacheKey,
    cache_dir: Path | str = DEFAULT_STRAIN_BASELINE_CACHE_DIR,
) -> dict[str, Any] | None:
    """命中则返回缓存载荷 dict；未命中 / 文件损坏 / 摘要不符 → None（当作未命中）。"""
    path = cache_path(key, cache_dir)
    if not path.exists():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError, UnicodeDecodeError):
        return None
    if not isinstance(payload, dict) or payload.get("cache_key_digest") != cache_key_digest(key):
        return None
    return payload


def store_baseline(
    key: StrainBaselineCacheKey,
    rows: Sequence[dict[str, Any]],
    cache_dir: Path | str = DEFAULT_STRAIN_BASELINE_CACHE_DIR,
    *,
    source_run: str | None = None,
    built_at: str | None = None,
    extra: dict[str, Any] | None = None,
) -> Path:
    path = cache_path(key, cache_dir)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = baseline_to_payload(key, rows, source_run=source_run, built_at=built_at, extra=extra)
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path


def cached_baseline_build(
    key: StrainBaselineCacheKey,
    compute: Callable[[], Sequence[dict[str, Any]]],
    *,
    cache_dir: Path | str = DEFAULT_STRAIN_BASELINE_CACHE_DIR,
    force: bool = False,
    source_run: str | None = None,
) -> tuple[dict[str, Any], bool]:
    """读穿缓存：命中返回 (载荷, True)；未命中 / force 则算并写缓存，返回 (载荷, False)。

    ``compute`` 只在未命中 / force 时调用，返回蒸馏好的候选行序列（全基因组筛查是 hour-scale，
    故几乎总在后台离线跑；本函数是编排薄壳，不含求解逻辑）。
    """
    if not force:
        cached = load_cached_baseline(key, cache_dir)
        if cached is not None:
            return cached, True
    rows = list(compute())
    store_baseline(key, rows, cache_dir, source_run=source_run)
    return load_cached_baseline(key, cache_dir) or baseline_to_payload(key, rows, source_run=source_run), False


__all__ = [
    "BASELINE_CACHE_SCHEMA_VERSION",
    "BASELINE_ROW_FIELDS",
    "DEFAULT_STRAIN_BASELINE_CACHE_DIR",
    "StrainBaselineCacheKey",
    "baseline_to_payload",
    "cache_key_digest",
    "cache_path",
    "cached_baseline_build",
    "distill_tradeoff_rows",
    "load_cached_baseline",
    "store_baseline",
    "strain_fingerprint",
]
