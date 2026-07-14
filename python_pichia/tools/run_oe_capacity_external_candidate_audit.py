from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
for path in (REPO_ROOT, REPO_ROOT / "python_pichia" / "src"):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from pcsec_pichia.oe_capacity.external_candidate_audit import (
    G6PDH2_GENE_ID,
    ExternalCapacityAuditRequest,
    ExternalCapacitySourceType,
    run_external_capacity_candidate_audit,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Audit Round 6A G6PDH2 external baseline-capacity candidates."
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--offline-replay", action="store_true")
    parser.add_argument(
        "--identity-cache-dir",
        type=Path,
        help="Existing UniProt identity cache to reuse for offline replay.",
    )
    parser.add_argument("--measurement-file", type=Path)
    parser.add_argument("--source-id", default="")
    parser.add_argument(
        "--source-type",
        choices=tuple(item.value for item in ExternalCapacitySourceType),
        default=ExternalCapacitySourceType.QUANTITATIVE_PROTEOMICS.value,
    )
    parser.add_argument("--source-version", default="")
    parser.add_argument("--source-url", default="")
    parser.add_argument("--license-id", default="")
    parser.add_argument("--license-url", default="")
    parser.add_argument("--query", default=f"{G6PDH2_GENE_ID} glucose mu=0.1")
    parser.add_argument("--expected-sha256", default="")
    parser.add_argument("--terms-reviewed", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    outputs = run_external_capacity_candidate_audit(
        ExternalCapacityAuditRequest(
            repo_root=REPO_ROOT,
            output_dir=args.output_dir,
            offline_replay=args.offline_replay,
            identity_cache_dir=args.identity_cache_dir,
            measurement_file=args.measurement_file,
            source_id=args.source_id,
            source_type=ExternalCapacitySourceType(args.source_type),
            source_version=args.source_version,
            source_url=args.source_url,
            license_id=args.license_id,
            license_url=args.license_url,
            query=args.query,
            expected_sha256=args.expected_sha256,
            terms_reviewed=args.terms_reviewed,
        )
    )
    print(json.dumps(outputs.summary()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
