from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
PYTHON_SRC = REPO_ROOT / "python_pichia" / "src"
if str(PYTHON_SRC) not in sys.path:
    sys.path.insert(0, str(PYTHON_SRC))

from pcsec_pichia.external_refs import (
    default_external_model_inventory_records,
    write_external_model_inventory,
)


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    output_dir = _resolve_output_dir(args.output_dir)
    records = default_external_model_inventory_records()
    outputs = write_external_model_inventory(records, output_dir)
    print(
        json.dumps(
            {
                "output_dir": str(output_dir),
                "jsonl_path": str(outputs.jsonl_path),
                "tsv_path": str(outputs.tsv_path),
                "report_path": str(outputs.report_path),
                "record_count": outputs.record_count,
                "download_status_counts": _status_counts(records),
            },
            ensure_ascii=False,
            indent=2,
        )
    )
    return 0


def _parse_args(argv: list[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build external GEM/GPR source inventory for pcSecPichia.")
    parser.add_argument(
        "--output-dir",
        default="",
        help="Output directory. Defaults to local_runs/external_model_gpr_inventory/<timestamp>.",
    )
    return parser.parse_args(argv)


def _resolve_output_dir(value: str) -> Path:
    if value:
        path = Path(value)
        return path if path.is_absolute() else REPO_ROOT / path
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return REPO_ROOT / "local_runs" / "external_model_gpr_inventory" / stamp


def _status_counts(records: tuple[object, ...]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for record in records:
        status = str(getattr(record, "download_status"))
        counts[status] = counts.get(status, 0) + 1
    return dict(sorted(counts.items()))


if __name__ == "__main__":
    raise SystemExit(main())
