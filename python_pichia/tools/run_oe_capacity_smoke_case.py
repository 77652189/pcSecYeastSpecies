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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run one controlled Phase 2 OE-capacity acceptance smoke case."
    )
    parser.add_argument("--target-id", required=True)
    parser.add_argument("--gene-id", required=True)
    parser.add_argument("--case-kind", choices=("executable", "boundary"), required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--run-name", required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    result = submit_oe_capacity_screen(
        target_id=args.target_id,
        gene_ids=(args.gene_id,),
        dose_payload={
            "dose_id": "2x",
            "dose_mode": "explicit_multiplier",
            "expression_multiplier": 2.0,
        },
        parameter_scenarios=("low", "nominal", "high"),
        execution_mode="comparison",
        product_mode="absolute_capacity",
        feature_enabled=True,
        compare_proxy=True,
        run_name=args.run_name,
        output_root=args.output_root,
        case_kind=args.case_kind,
    )
    print(json.dumps(result, ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
