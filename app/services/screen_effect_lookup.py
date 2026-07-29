"""把已跑筛查结果里的**真实效应**接给候选选择器。

动机：选择器此前只按机制分组给候选，"该先试哪个"仍要靠猜。而策展 scope 的筛查是**分钟级**的，
跑完就有每个候选对每个靶点的真实相对效应——把这份数字接进选择器，挑候选就从"按机制猜"
变成"按模型算的效应排"。

诚实边界：效应是**模型内部相对量、不是 titer 预测**；KO 必须连生长代价一起看（实测常见
"提升几个点、生长掉一半"的陷阱候选）。没有跑过筛查时全部降级为空，选择器仍可用。
"""

from __future__ import annotations

import csv
from functools import lru_cache
from pathlib import Path
from typing import Any

from app.services.genome_wide_screen_registry import list_runs, run_scope_family
from app.ui.common import PATHS


def _candidate_csv_paths(target_id: str, paths: Any | None = None) -> list[str]:
    """能覆盖该靶点的筛查结果 CSV，策展 scope 优先（它才含复合体反应）、新的在前。"""
    resolved = paths or PATHS
    runs = [run for run in list_runs(resolved) if run.csv_path and target_id in (run.targets or [])]
    runs.sort(key=lambda run: (run_scope_family(run) != "catalog", run.updated_at or ""), reverse=False)
    # 策展在前；同类里新的在前
    catalog = [run for run in runs if run_scope_family(run) == "catalog"]
    others = [run for run in runs if run_scope_family(run) != "catalog"]
    catalog.sort(key=lambda run: run.updated_at or "", reverse=True)
    others.sort(key=lambda run: run.updated_at or "", reverse=True)
    return [str(run.csv_path) for run in (*catalog, *others)]


@lru_cache(maxsize=8)
def _read_effects(csv_path: str, mtime: float, target_id: str) -> dict[tuple[str, str], tuple[float, float]]:
    """{(KO/OE, 模型对象 id): (相对效应, 生长保持)}。策展行的 gene_id 列存的就是反应 id。"""
    effects: dict[tuple[str, str], tuple[float, float]] = {}
    try:
        with Path(csv_path).open(encoding="utf-8-sig", newline="") as handle:
            for row in csv.DictReader(handle):
                if str(row.get("target_id") or "") != target_id:
                    continue
                object_id = str(row.get("gene_id") or "").strip()
                intervention = str(row.get("intervention_type") or "").strip().upper()
                if not object_id or intervention not in {"KO", "OE"}:
                    continue
                try:
                    ratio = float(row.get("secretion_ratio_vs_wildtype") or "")
                except (TypeError, ValueError):
                    continue
                try:
                    growth = float(row.get("growth_retention_ratio") or "")
                except (TypeError, ValueError):
                    growth = 1.0
                key = (intervention, object_id)
                # 同一候选在多份结果里出现时保留第一份（调用方已按策展优先 / 新在前排序）。
                effects.setdefault(key, (ratio - 1.0, growth))
    except (OSError, UnicodeDecodeError, csv.Error):
        return {}
    return effects


def available_screen_targets(paths: Any | None = None) -> list[str]:
    """已有筛查结果覆盖了哪些靶点。"""
    resolved = paths or PATHS
    targets: set[str] = set()
    for run in list_runs(resolved):
        if run.csv_path:
            targets.update(str(item) for item in (run.targets or []) if str(item).strip())
    return sorted(targets)


def resolve_effect_source(
    target_id: str, target_context: str | None = None, *, paths: Any | None = None
) -> str:
    """把界面上的 target_id 解析成**筛查结果里真实存在的靶点**，解析不出则返回空串。

    为什么需要这层：三段式 / 自定义构建的 target_id 是拼出来的复合串
    （如 `alpha-factor_MRFPS_OPN_alpha-pro_OPN_ALPHA_FULL_PROJECT`），跟筛查按规范靶点存的结果对不上，
    直接查会得到 0 条、效应列整列消失（2026-07-28 实测到的回归）。而且筛查里的 OPN 靶点实际叫
    `OPN_ALPHA_FULL_PROJECT` 而非 `OPN`，所以只按上下文精确匹配也不够，要做子串匹配。

    解析出的靶点可能与当前正在构建的目标**不完全相同**（例如自定义三段式组合 vs 规范 OPN），
    调用方应把这一点如实告诉用户，不要让人以为效应就是这个自定义构建体算出来的。
    """
    targets = available_screen_targets(paths)
    if target_id in targets:
        return target_id
    context = str(target_context or "").strip().upper()
    if context:
        exact = [item for item in targets if item.upper() == context]
        if exact:
            return exact[0]
        partial = [item for item in targets if context in item.upper()]
        if partial:
            return partial[0]
    return ""


def load_screen_effect_lookup(target_id: str, *, paths: Any | None = None) -> dict[tuple[str, str], tuple[float, float]]:
    """合并该靶点所有可用筛查结果（策展优先）。没跑过筛查 → 返回空 dict。"""
    merged: dict[tuple[str, str], tuple[float, float]] = {}
    for csv_path in _candidate_csv_paths(target_id, paths):
        path = Path(csv_path)
        if not path.exists():
            continue
        for key, value in _read_effects(csv_path, path.stat().st_mtime, target_id).items():
            merged.setdefault(key, value)
    return merged


__all__ = ["available_screen_targets", "load_screen_effect_lookup", "resolve_effect_source"]
