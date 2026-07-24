"""B2：把候选短名单跑过"真实工艺条件矩阵"，缓存每个候选在各条件下的相对分泌效应。

短名单是在**默认条件（葡萄糖 μ=0.10）**的筛查里选出来的 OE 提升候选。B2 问的是：这些 top
候选在 hLF/OPN 的**其它真实工艺碳源条件**下还成立吗？——即"跨条件稳健性"的原始数据。

条件矩阵（干净单碳源操作点，统一 μ=0.10；混合/过渡条件因单一生物量近似不进排序比较）:
    hLF  -> 甘油(生长相) / 葡萄糖(生产相)
    OPN  -> 甲醇
（μ 统一用默认 0.10 = hLF 甘油生长相实测值，见 ADR-006；生产相虽慢但按工艺决定不作模型 μ。）

复用引擎的 `enable_oe_dose_response` 管道（把短名单反应作为 `oe_dose_response_reactions` 传进去，
每个条件换 `carbon_source_id`），不写新引擎代码。离线/后台跑一次，缓存
`{target}_condition_matrix.json`；分类判定留 B3、面板展示留 B4（本工具**只产矩阵数据**）。

用法（从 python_pichia/，src/ 在 PYTHONPATH）:
    python tools/run_shortlist_condition_matrix.py \
        --screen-csv ../local_runs/catalog_reaction_screen_hlf_opn_expanded61_v2/catalog_reaction_tradeoff_rows.csv \
        --output-dir ../local_runs/candidate_shortlist_readout
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import pandas as pd

_TOOLS_DIR = Path(__file__).resolve().parent
_SRC_DIR = _TOOLS_DIR.parent / "src"
for _p in (str(_SRC_DIR), str(_TOOLS_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from build_candidate_shortlist_readout import _oe_shortlist  # noqa: E402 - sibling tool, shared shortlist logic
from pcsec_pichia.engines.base import PichiaSimulationRequest  # noqa: E402
from pcsec_pichia.pipeline import run_pichia_secretion_simulation  # noqa: E402

# 每个目标的干净单碳源真实工艺条件（μ 统一默认 0.10）。
DEFAULT_TARGET_CONDITIONS: dict[str, tuple[str, ...]] = {
    "hLF": ("glycerol", "glucose"),
    "OPN_ALPHA_FULL_PROJECT": ("methanol",),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--screen-csv", type=Path, required=True, help="筛查 tradeoff CSV（复用，短名单来源）")
    parser.add_argument("--output-dir", type=Path, required=True, help="写 {target}_condition_matrix.json 的目录")
    parser.add_argument("--targets", nargs="+", default=list(DEFAULT_TARGET_CONDITIONS))
    parser.add_argument("--top-n", type=int, default=8, help="短名单取前 N 个候选（有界求解，别太大）")
    parser.add_argument("--mu", type=float, default=0.10, help="固定生长速率（默认 0.10 = hLF 甘油生长相实测值）")
    parser.add_argument("--factors", nargs="+", type=float, default=None, help="OE 倍数网格（默认用引擎默认）")
    return parser.parse_args()


def _run_condition(target_id: str, reactions: list[str], carbon_source_id: str, mu: float, factors, output_dir: Path) -> dict:
    """在单个碳源条件下跑短名单 OE 剂量响应，返回该条件的原始效应数据。"""
    t0 = time.time()
    result = run_pichia_secretion_simulation(
        PichiaSimulationRequest(
            target_id=target_id,
            candidate_id=target_id,
            carbon_source_id=carbon_source_id,
            mu=float(mu),
            enable_oe_dose_response=True,
            oe_dose_response_reactions=tuple(reactions),
            oe_dose_response_factors=tuple(factors) if factors else (),
        ),
        output_dir=output_dir / "_engine_runs" / f"{target_id}_{carbon_source_id}",
    )
    elapsed = time.time() - t0
    if not result.success or result.summary_path is None:
        return {"carbon_source_id": carbon_source_id, "success": False, "elapsed_seconds": round(elapsed, 1)}
    summary = json.loads(result.summary_path.read_text(encoding="utf-8"))
    odr = (summary.get("protein_cost_analysis") or {}).get("oe_dose_response") or {}
    reaction_shapes = odr.get("reaction_shapes") or []
    return {
        "carbon_source_id": carbon_source_id,
        "mu": float(mu),
        "success": bool(odr.get("success")),
        "result_status": odr.get("result_status"),
        "baseline_objective": odr.get("baseline_objective"),
        "reaction_shapes": reaction_shapes,
        "shapes_by_reaction": {str(s.get("reaction_id")): s for s in reaction_shapes if isinstance(s, dict)},
        "warnings": odr.get("warnings"),
        "elapsed_seconds": round(elapsed, 1),
    }


def matrix_for_target(csv_path: Path, target_id: str, conditions: tuple[str, ...], top_n: int, mu: float, factors, output_dir: Path) -> None:
    frame = pd.read_csv(csv_path)
    shortlist = _oe_shortlist(frame, target_id, top_n)
    reactions = [str(row["reaction"]) for row in shortlist]
    if not reactions:
        print(f"[{target_id}] 短名单无 OE 提升候选反应，跳过。")
        return

    condition_rows = [_run_condition(target_id, reactions, cs, mu, factors, output_dir) for cs in conditions]

    # 跨条件的 per-reaction 视图（B4 面板方便消费：每个候选在各条件下的 shape/benefit）。
    per_reaction: dict[str, dict[str, object]] = {}
    for reaction in reactions:
        per_reaction[reaction] = {
            row["carbon_source_id"]: row.get("shapes_by_reaction", {}).get(reaction)
            for row in condition_rows
            if row.get("success")
        }

    payload = {
        "target_id": target_id,
        "mu": float(mu),
        "conditions": [str(c) for c in conditions],
        "shortlist": shortlist,
        "swept_reactions": reactions,
        "per_condition": condition_rows,
        "per_reaction_across_conditions": per_reaction,
        "note": "只产条件矩阵原始数据；跨条件稳健性分类见 B3、面板展示见 B4。μ 统一默认 0.10（见 ADR-006）。",
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    out_path = output_dir / f"{target_id}_condition_matrix.json"
    out_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2, default=str), encoding="utf-8")
    ok = sum(1 for row in condition_rows if row.get("success"))
    print(f"[{target_id}] {len(reactions)} 候选 × {len(conditions)} 条件（成功 {ok}/{len(conditions)}） -> {out_path.name}")


def main() -> int:
    args = parse_args()
    for target_id in args.targets:
        conditions = DEFAULT_TARGET_CONDITIONS.get(target_id, ("glucose",))
        matrix_for_target(args.screen_csv, target_id, conditions, args.top_n, args.mu, args.factors, args.output_dir)
    print("DONE ->", args.output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
