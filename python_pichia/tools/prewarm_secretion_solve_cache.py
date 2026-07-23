"""B1 预热：把常见 (目标 × 碳源 × μ) 分泌容量求解结果预热进本地缓存。

分泌 LP 求解慢且确定（`highs-ds`）。本工具离线/后台跑一次，把常见条件解出来存进
`solve_cache`（默认 `local_runs/solve_cache/`，gitignored）；之后跨条件稳健性 / 面板等
读穿缓存即可零重复求解（见 `pcsec_pichia.solve_cache`）。

默认条件集 = {hLF, OPN} × {glucose, glycerol, methanol} × {μ=0.10} ∪ hLF 生产相锚点
（glucose, μ=0.013，见 `pcsec_pichia.process_anchors`）。不改任何求解语义；只记忆化确定性结果。

用法（从 python_pichia/，src/ 在 PYTHONPATH）:
    python tools/prewarm_secretion_solve_cache.py
    python tools/prewarm_secretion_solve_cache.py --cache-dir ../local_runs/solve_cache --force
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

_TOOLS_DIR = Path(__file__).resolve().parent
_SRC_DIR = _TOOLS_DIR.parent / "src"
_REPO_ROOT = _TOOLS_DIR.parents[1]
for _p in (str(_SRC_DIR), str(_TOOLS_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from pcsec_pichia.loading import load_pcsec_pichia_inputs  # noqa: E402
from pcsec_pichia.process_anchors import list_process_growth_anchors  # noqa: E402
from pcsec_pichia.simulation import DEFAULT_SOLVER_METHOD, solve_secretion_capacity  # noqa: E402
from pcsec_pichia.solve_cache import (  # noqa: E402
    DEFAULT_SOLVE_CACHE_DIR,
    SecretionSolveCacheKey,
    cached_secretion_solve,
)
from pcsec_pichia.targets import load_builtin_targets  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--targets", nargs="+", default=["hLF", "OPN_ALPHA_FULL_PROJECT"])
    parser.add_argument("--carbon-sources", nargs="+", default=["glucose", "glycerol", "methanol"])
    parser.add_argument("--mu", nargs="+", type=float, default=[0.10], help="生长速率网格")
    parser.add_argument("--media-type", type=int, default=4)
    parser.add_argument("--cache-dir", type=Path, default=None, help=f"默认 {DEFAULT_SOLVE_CACHE_DIR}")
    parser.add_argument("--include-anchors", action="store_true", default=True, help="附加 process_anchors 工艺操作点（默认开）")
    parser.add_argument("--no-anchors", dest="include_anchors", action="store_false")
    parser.add_argument("--force", action="store_true", help="忽略已有缓存、强制重算覆盖")
    return parser.parse_args()


def _conditions(args: argparse.Namespace) -> list[tuple[str, str, float]]:
    """去重的 (target, carbon_source, mu) 条件集。"""
    conditions: list[tuple[str, str, float]] = []
    for target_id in args.targets:
        for carbon_source_id in args.carbon_sources:
            for mu in args.mu:
                conditions.append((target_id, carbon_source_id, float(mu)))
    if args.include_anchors:
        for anchor in list_process_growth_anchors():
            if anchor.target_family in args.targets:
                conditions.append((anchor.target_family, anchor.carbon_source_id, float(anchor.growth_rate)))
    # 去重、保序
    seen: set[tuple[str, str, float]] = set()
    unique: list[tuple[str, str, float]] = []
    for cond in conditions:
        if cond not in seen:
            seen.add(cond)
            unique.append(cond)
    return unique


def main() -> int:
    args = parse_args()
    cache_dir = args.cache_dir or DEFAULT_SOLVE_CACHE_DIR
    conditions = _conditions(args)
    targets = {t.target_id: t for t in load_builtin_targets(_REPO_ROOT)}
    inputs_by_carbon: dict[str, object] = {}

    computed = 0
    hit = 0
    failed = 0
    print(f"预热 {len(conditions)} 个条件 → {cache_dir}（force={args.force}）")
    for target_id, carbon_source_id, mu in conditions:
        target = targets.get(target_id)
        if target is None:
            print(f"  跳过未知目标 {target_id}")
            continue
        if carbon_source_id not in inputs_by_carbon:
            inputs_by_carbon[carbon_source_id] = load_pcsec_pichia_inputs(
                _REPO_ROOT, media_type=args.media_type, carbon_source_id=carbon_source_id
            )
        inp = inputs_by_carbon[carbon_source_id]
        key = SecretionSolveCacheKey(
            target_id=target_id,
            carbon_source_id=carbon_source_id,
            media_type=int(args.media_type),
            growth_rate=float(mu),
            write_ribosome_translation_constraint=False,
            write_misfolding_constraints=False,
            solver_method=DEFAULT_SOLVER_METHOD,
        )

        def compute():
            return solve_secretion_capacity(
                inp.prepared_model, target, inp.amino_acids, inp.metabolic, inp.secretory, inp.combined, growth_rate=float(mu)
            )

        started = time.time()
        result, from_cache = cached_secretion_solve(key, compute, cache_dir=cache_dir, force=args.force)
        elapsed = time.time() - started
        if not result.success:
            failed += 1
            tag = "求解失败(不缓存)"
        elif from_cache:
            hit += 1
            tag = "命中缓存"
        else:
            computed += 1
            tag = "已求解并缓存"
        print(f"  {target_id:22s} {carbon_source_id:16s} μ={mu:<6.4g} {tag:14s} obj={result.objective_value!r} ({elapsed:.1f}s)")

    print(f"完成：已求解 {computed} / 命中 {hit} / 失败 {failed}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
