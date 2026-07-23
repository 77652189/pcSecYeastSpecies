"""分泌求解结果缓存（内容寻址、默认读缓存、显式重算）。

分泌 LP 求解慢（载模型 + 解 LP）且**确定**（`highs-ds` 确定性），所以对相同语义输入
反复求解是浪费。本模块按语义输入的内容哈希做键，把结果 JSON 存本地缓存目录：
默认命中即读、未命中才算并写；显式 `force` 时重算覆盖。

设计边界：
- **不改任何求解语义**——只是把确定性求解的结果记忆化。缓存是本地派生产物
  （默认 `local_runs/solve_cache/`，已 gitignored），可随时删除重建。
- **键必须涵盖决定结果的全部语义输入**，否则会返回过期/错误结果。求解逻辑或模型数据
  变更时 **bump `CACHE_SCHEMA_VERSION`**（旧键即失配，不会误命中）。
- 改了模型（KO/OE 等）必须相应设置 `model_variant_fingerprint`，否则改动后的求解会
  误命中野生型缓存条目。
- 只缓存 `success=True` 的结果——避免把瞬态/环境失败固化进缓存。
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable

if TYPE_CHECKING:
    from pcsec_pichia.simulation import SecretionSimulationResult

# 求解语义 / 模型数据变更时必须 bump（碳源标定接线后的当前版本）。
CACHE_SCHEMA_VERSION = "2026-07-23-carbon-source-calibrated"
DEFAULT_SOLVE_CACHE_DIR = Path("local_runs") / "solve_cache"

# SecretionSimulationResult 中 JSON 往返会退化成 list、需还原成 tuple 的字段。
_RESULT_TUPLE_FIELDS = ("open_growth_reaction_ids", "warnings")


@dataclass(frozen=True)
class SecretionSolveCacheKey:
    """决定一次分泌容量求解结果的全部语义输入（内容寻址键）。"""

    target_id: str
    carbon_source_id: str
    media_type: int
    growth_rate: float
    write_ribosome_translation_constraint: bool
    write_misfolding_constraints: bool
    solver_method: str
    # 自定义目标的参数指纹（序列 / 二硫键 / 糖基化位点等影响求解的）；内建目标可留空
    # （由 target_id + schema_version 唯一确定）。
    target_fingerprint: str = ""
    # 模型变体指纹：KO/OE 等对模型的改动。未改动 = "wildtype"。
    # 改了模型必须改这个，否则改动后的求解会误命中野生型缓存！
    model_variant_fingerprint: str = "wildtype"
    entrypoint: str = "secretion_capacity"
    schema_version: str = CACHE_SCHEMA_VERSION


def _canonical_json(obj: Any) -> str:
    return json.dumps(obj, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def cache_key_digest(key: SecretionSolveCacheKey) -> str:
    """键的稳定 sha256 摘要（内容寻址）。"""
    return hashlib.sha256(_canonical_json(asdict(key)).encode("utf-8")).hexdigest()


def cache_path(key: SecretionSolveCacheKey, cache_dir: Path | str = DEFAULT_SOLVE_CACHE_DIR) -> Path:
    return Path(cache_dir) / f"{key.entrypoint}_{cache_key_digest(key)}.json"


def secretion_result_to_payload(
    key: SecretionSolveCacheKey,
    result: "SecretionSimulationResult",
) -> dict[str, Any]:
    return {
        "schema_version": key.schema_version,
        "cache_key": asdict(key),
        "cache_key_digest": cache_key_digest(key),
        "result": asdict(result),
    }


def _coerce_result_fields(data: dict[str, Any]) -> dict[str, Any]:
    coerced = dict(data)
    for field_name in _RESULT_TUPLE_FIELDS:
        value = coerced.get(field_name)
        if isinstance(value, list):
            coerced[field_name] = tuple(value)
    lp_sensitivity = coerced.get("lp_sensitivity")
    if isinstance(lp_sensitivity, dict):
        coerced["lp_sensitivity"] = {
            name: tuple(value) if isinstance(value, list) else value
            for name, value in lp_sensitivity.items()
        }
    return coerced


def secretion_result_from_payload(payload: dict[str, Any]) -> "SecretionSimulationResult":
    from pcsec_pichia.simulation import SecretionSimulationResult

    return SecretionSimulationResult(**_coerce_result_fields(dict(payload["result"])))


def load_cached_secretion_result(
    key: SecretionSolveCacheKey,
    cache_dir: Path | str = DEFAULT_SOLVE_CACHE_DIR,
) -> "SecretionSimulationResult | None":
    """命中则返回缓存的结果；未命中 / 文件损坏 / 摘要不符 → None（当作未命中）。"""
    path = cache_path(key, cache_dir)
    if not path.exists():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return None
    if not isinstance(payload, dict) or payload.get("cache_key_digest") != cache_key_digest(key):
        return None
    return secretion_result_from_payload(payload)


def store_secretion_result(
    key: SecretionSolveCacheKey,
    result: "SecretionSimulationResult",
    cache_dir: Path | str = DEFAULT_SOLVE_CACHE_DIR,
) -> Path:
    path = cache_path(key, cache_dir)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = secretion_result_to_payload(key, result)
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path


def cached_secretion_solve(
    key: SecretionSolveCacheKey,
    compute: Callable[[], "SecretionSimulationResult"],
    *,
    cache_dir: Path | str = DEFAULT_SOLVE_CACHE_DIR,
    force: bool = False,
) -> tuple["SecretionSimulationResult", bool]:
    """读穿缓存：命中返回 (缓存结果, True)；未命中 / force 则算并（仅 success）写缓存，返回 (结果, False)。"""
    if not force:
        cached = load_cached_secretion_result(key, cache_dir)
        if cached is not None:
            return cached, True
    result = compute()
    if getattr(result, "success", False):
        store_secretion_result(key, result, cache_dir)
    return result, False


__all__ = [
    "CACHE_SCHEMA_VERSION",
    "DEFAULT_SOLVE_CACHE_DIR",
    "SecretionSolveCacheKey",
    "cache_key_digest",
    "cache_path",
    "cached_secretion_solve",
    "load_cached_secretion_result",
    "secretion_result_from_payload",
    "secretion_result_to_payload",
    "store_secretion_result",
]
