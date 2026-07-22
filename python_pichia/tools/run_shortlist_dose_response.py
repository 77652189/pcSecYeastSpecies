"""① R2 增强：对候选短名单跑 OE 剂量响应扫描，缓存每个候选的形状。

短名单本身零求解（复用筛查缓存），但"剂量响应形状"（过表达越多，分泌是持续上升＝线性，
还是很快到顶＝饱和）需要对短名单每个反应扫多个 OE 倍数、每个倍数重解目标 LP——**有界的
额外求解**（仅 top 短名单几个反应，不是全基因组 1000+）。因此这一步离线/后台跑一次，把
形状缓存成 `{target}_dose_response.json`，Streamlit 面板只读缓存、零求解显示。

复用引擎已有的 `enable_oe_dose_response` 管道（把短名单反应 ID 作为 `oe_dose_response_reactions`
传进去），不写新引擎代码。只开 dose-response、不开 cost_slope（protein_cost_analysis 只要
oe_dose_response 非 None 就会填充，省一轮 cost_slope 求解）。

用法（从 python_pichia/，src/ 在 PYTHONPATH）:
    python tools/run_shortlist_dose_response.py \
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--screen-csv", type=Path, required=True, help="筛查 tradeoff CSV（复用，已求解）")
    parser.add_argument("--output-dir", type=Path, required=True, help="写 {target}_dose_response.json 的目录")
    parser.add_argument("--targets", nargs="+", default=["hLF", "OPN_ALPHA_FULL_PROJECT"])
    parser.add_argument("--top-n", type=int, default=8, help="短名单取前 N 个候选去扫（有界求解，别太大）")
    parser.add_argument(
        "--factors",
        nargs="+",
        type=float,
        default=None,
        help="OE 倍数网格（默认用引擎 DEFAULT_OE_DOSE_RESPONSE_FACTORS）",
    )
    return parser.parse_args()


def sweep_target(csv_path: Path, target_id: str, top_n: int, factors, output_dir: Path) -> None:
    frame = pd.read_csv(csv_path)
    shortlist = _oe_shortlist(frame, target_id, top_n)
    reactions = [str(row["reaction"]) for row in shortlist]
    if not reactions:
        print(f"[{target_id}] 短名单无 OE 提升候选反应，跳过。")
        return

    t0 = time.time()
    result = run_pichia_secretion_simulation(
        PichiaSimulationRequest(
            target_id=target_id,
            candidate_id=target_id,
            enable_oe_dose_response=True,
            oe_dose_response_reactions=tuple(reactions),
            oe_dose_response_factors=tuple(factors) if factors else (),
        ),
        output_dir=output_dir / "_engine_runs" / target_id,
    )
    elapsed = time.time() - t0
    if not result.success or result.summary_path is None:
        print(f"[{target_id}] 管道未成功（success={result.success}），跳过缓存。")
        return

    summary = json.loads(result.summary_path.read_text(encoding="utf-8"))
    odr = (summary.get("protein_cost_analysis") or {}).get("oe_dose_response") or {}
    reaction_shapes = odr.get("reaction_shapes") or []
    shapes_by_reaction = {str(s.get("reaction_id")): s for s in reaction_shapes if isinstance(s, dict)}

    payload = {
        "target_id": target_id,
        "result_status": odr.get("result_status"),
        "success": odr.get("success"),
        "tested_factors": odr.get("tested_factors"),
        "baseline_objective": odr.get("baseline_objective"),
        "swept_reactions": reactions,
        "reaction_shapes": reaction_shapes,
        "shapes_by_reaction": shapes_by_reaction,
        "warnings": odr.get("warnings"),
        "elapsed_seconds": round(elapsed, 1),
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    out_path = output_dir / f"{target_id}_dose_response.json"
    out_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2, default=str), encoding="utf-8")
    shape_counts: dict[str, int] = {}
    for s in reaction_shapes:
        shape_counts[str(s.get("shape"))] = shape_counts.get(str(s.get("shape")), 0) + 1
    print(
        f"[{target_id}] 扫了 {len(reactions)} 个反应、{elapsed:.1f}s -> {out_path.name}；"
        f"形状分布={shape_counts}"
    )


def main() -> int:
    args = parse_args()
    for target_id in args.targets:
        sweep_target(args.screen_csv, target_id, args.top_n, args.factors, args.output_dir)
    print("DONE ->", args.output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
