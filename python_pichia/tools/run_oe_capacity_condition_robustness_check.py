from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
for path in (REPO_ROOT, REPO_ROOT / "python_pichia" / "src"):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from app.services.pichia_oe_capacity_service import submit_oe_capacity_screen
from pcsec_pichia.oe_capacity import evaluate_condition_robustness, ranking_case_from_screen_result


"""Direction-5-as-acceptance-dimension check: is direction 2's relative OE
ranking (the same relative_vs_baseline_delta the UI now shows as its
headline comparison, see app/ui/views/oe_capacity.py._render_row_comparison)
stable across a couple of different growth conditions, for a small,
already-attested gene shortlist?

Deliberately NOT a Streamlit page or a new service: EXECUTION_PLAN.md's
"未授权" section only authorizes direction 5 "作为方向 2-3 的局部验收维度"
(a local acceptance dimension of direction 2-3), not a full cross-condition
ranking product. This script calls the exact same, already-parameterized
submit_oe_capacity_screen() the OE capacity page calls, once per
(carbon_source_id, growth_rate) pair; nothing here is a new capability.
"""


DEFAULT_GENE_IDS = ("PAS_chr2-1_0308", "PAS_chr1-4_0458", "PAS_chr2-1_0047")
# glucose/mu=0.1 is the reviewed-formulation default every other OE capacity
# run in this repo uses; mu=0.15 changes only growth rate on the same
# reference-quality carbon source; glycerol/mu=0.1 changes only carbon
# source. media.list_carbon_source_formulations() marks glycerol
# "draft_carbon_source_boundary" (not glucose's "corrected_reference"), so
# that leg's result must be read as a boundary probe, not a calibrated
# alternate condition.
DEFAULT_CONDITIONS: tuple[tuple[str, float], ...] = (
    ("glucose", 0.1),
    ("glucose", 0.15),
    ("glycerol", 0.1),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Re-run direction 2's relative OE screen for a small gene set across a "
            "few (carbon_source, growth_rate) contexts and report ranking stability."
        )
    )
    parser.add_argument("--target-id", default="hLF")
    parser.add_argument("--gene-id", dest="gene_ids", action="append", default=None)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--run-name-prefix", default="condition-robustness")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    gene_ids = tuple(args.gene_ids) if args.gene_ids else DEFAULT_GENE_IDS

    cases = []
    per_context_results = {}
    for carbon_source_id, growth_rate in DEFAULT_CONDITIONS:
        context_id = f"{carbon_source_id}_mu_{growth_rate:g}"
        run_name = f"{args.run_name_prefix}-{args.target_id}-{carbon_source_id}-mu{growth_rate:g}"
        result = submit_oe_capacity_screen(
            target_id=args.target_id,
            gene_ids=gene_ids,
            dose_payload={
                "dose_id": "2x",
                "dose_mode": "explicit_multiplier",
                "expression_multiplier": 2.0,
            },
            parameter_scenarios=("nominal",),
            execution_mode="comparison",
            product_mode="relative_uncalibrated",
            feature_enabled=True,
            compare_proxy=False,
            growth_rate=growth_rate,
            carbon_source_id=carbon_source_id,
            run_name=run_name,
            output_root=args.output_root,
        )
        per_context_results[context_id] = result
        cases.append(
            ranking_case_from_screen_result(
                context_id=context_id,
                carbon_source_id=carbon_source_id,
                growth_rate=growth_rate,
                screen_result=result,
            )
        )

    robustness = evaluate_condition_robustness(target_id=args.target_id, gene_ids=gene_ids, cases=tuple(cases))

    report = {
        "target_id": robustness.target_id,
        "gene_ids": list(robustness.gene_ids),
        "full_order_is_stable": robustness.full_order_is_stable,
        "top1_is_stable": robustness.top1_is_stable,
        "baseline_context": robustness.baseline_case.context_id,
        "baseline_ranking": list(robustness.baseline_case.ranked_gene_ids),
        "cases": [
            {
                "context_id": case.context_id,
                "carbon_source_id": case.carbon_source_id,
                "growth_rate": case.growth_rate,
                "ranked_gene_ids": list(case.ranked_gene_ids),
                "gene_deltas": [[gene_id, delta] for gene_id, delta in case.gene_deltas],
                "missing_gene_ids": list(case.missing_gene_ids),
            }
            for case in robustness.cases
        ],
        "unstable_context_ids": [case.context_id for case in robustness.unstable_cases],
        "top1_unstable_context_ids": [case.context_id for case in robustness.top1_unstable_cases],
    }
    print(json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True))

    report_path = Path(args.output_root) / f"{args.run_name_prefix}-{args.target_id}-report.json"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
