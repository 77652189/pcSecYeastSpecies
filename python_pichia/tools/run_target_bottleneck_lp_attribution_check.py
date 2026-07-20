"""Direction-3-as-acceptance-dimension check (not a new screening product).

Question: at the current model's optimum, which constraint block is the
real binding bottleneck on hLF/OPN secretion - and does that bottleneck
shift to a different resource as forced target-expression level is swept
up/down? Both pieces already exist and are already tested
(analyze_target_protein_lp_attribution is called unconditionally in
pipeline.py's run_pichia_secretion_simulation; run_protein_cost_slope_compatibility
is the opt-in `enable_cost_slope_compatibility` sweep exercised by
test_pipeline_entrypoints.py's cost-slope test). This script just runs the
existing, already-shipped pipeline entrypoint with that flag on for both
targets and surfaces the LP-sensitivity-derived tables it already produces,
so we don't have to build a new bottleneck-discovery tool from scratch.

Deliberately reuses default cost_slope_growth_rates/cost_slope_capacity_fractions
(engines/base.py) rather than inventing new sweep points.

Usage (from python_pichia/, with src/ on PYTHONPATH):
    python tools/run_target_bottleneck_lp_attribution_check.py --target-id hLF --output-root ../local_runs/target_bottleneck_lp_attribution
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

_SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if str(_SRC_DIR) not in sys.path:
    sys.path.insert(0, str(_SRC_DIR))

from pcsec_pichia.engines.base import PichiaSimulationRequest  # noqa: E402
from pcsec_pichia.pipeline import run_pichia_secretion_simulation  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-id", default="hLF")
    parser.add_argument("--output-root", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    candidate_id = args.target_id
    output_dir = Path(args.output_root) / args.target_id

    t0 = time.time()
    result = run_pichia_secretion_simulation(
        PichiaSimulationRequest(
            target_id=args.target_id,
            candidate_id=candidate_id,
            enable_cost_slope_compatibility=True,
        ),
        output_dir=output_dir,
    )
    print(f"[{time.strftime('%H:%M:%S')}] {args.target_id}: solved in {time.time() - t0:.1f}s, success={result.success}")

    summary = json.loads(result.summary_path.read_text(encoding="utf-8"))
    protein_cost = summary.get("protein_cost_analysis") or {}
    lp_attribution = protein_cost.get("lp_attribution") or {}
    cost_slope = protein_cost.get("cost_slope_compatibility") or {}

    report = {
        "target_id": args.target_id,
        "objective_evidence": lp_attribution.get("objective_evidence"),
        "dominant_constraint_blocks": lp_attribution.get("dominant_constraint_blocks"),
        "top_constraint_marginals": lp_attribution.get("top_constraint_marginals"),
        "top_bound_marginals": lp_attribution.get("top_bound_marginals"),
        "lp_attribution_warnings": lp_attribution.get("warnings"),
        "cost_slope_secretion_ratio_policy": cost_slope.get("secretion_ratio_policy"),
        "cost_slope_capacity_reference": cost_slope.get("capacity_reference"),
        "cost_slope_rows": cost_slope.get("rows"),
        "glucose_cost_slopes": cost_slope.get("glucose_cost_slopes"),
        "ribosome_cost_slopes": cost_slope.get("ribosome_cost_slopes"),
        "cost_slope_warnings": cost_slope.get("warnings"),
    }
    text = json.dumps(report, ensure_ascii=False, indent=2, default=str)
    print(text)
    out_path = Path(args.output_root) / f"target_bottleneck_lp_attribution_{args.target_id}.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text, encoding="utf-8")
    print(f"[{time.strftime('%H:%M:%S')}] wrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
