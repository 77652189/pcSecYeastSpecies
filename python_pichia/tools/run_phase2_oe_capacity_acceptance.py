from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Verify Phase 2 OE-capacity artifacts and emit the formal gate."
    )
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=Path(__file__).resolve().parents[2],
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=Path("local_runs/oe_capacity/round6/acceptance"),
    )
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    repo_root = args.repo_root.resolve()
    source_root = repo_root / "python_pichia" / "src"
    sys.path.insert(0, str(source_root))
    from pcsec_pichia.oe_capacity import run_phase2_oe_capacity_acceptance

    output_root = args.output_root
    if not output_root.is_absolute():
        output_root = repo_root / output_root
    summary = run_phase2_oe_capacity_acceptance(repo_root, output_root)
    print(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True))
    return 0 if summary["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
