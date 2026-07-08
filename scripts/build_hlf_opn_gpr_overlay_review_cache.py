from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from app.services.pichia_gene_catalog_service import (
    hlf_opn_gpr_overlay_review_cache_path,
    hlf_opn_gpr_overlay_review_summary,
    load_hlf_opn_gpr_overlay_review,
)
from pcsec_pichia.core.paths import ProjectPaths


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    paths = ProjectPaths(repo_root=Path(args.root).resolve())
    rows = load_hlf_opn_gpr_overlay_review(force_refresh=args.force_refresh, paths=paths)
    summary = hlf_opn_gpr_overlay_review_summary(paths=paths)
    print(
        json.dumps(
            {
                "cache_path": str(hlf_opn_gpr_overlay_review_cache_path(paths)),
                "row_count": len(rows),
                **summary,
            },
            ensure_ascii=False,
            indent=2,
        )
    )
    return 0


def _parse_args(argv: list[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build review-only hLF/OPN GPR overlay candidate cache.")
    parser.add_argument("--root", default=str(REPO_ROOT))
    parser.add_argument("--force-refresh", action="store_true")
    return parser.parse_args(argv)


if __name__ == "__main__":
    raise SystemExit(main())
