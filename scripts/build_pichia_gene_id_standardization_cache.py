from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from app.services.pichia_gene_catalog_service import (
    load_pichia_gene_id_standardization,
    pichia_gene_id_standardization_cache_path,
    pichia_gene_id_standardization_summary,
)
from pcsec_pichia.core.paths import ProjectPaths


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    paths = ProjectPaths(repo_root=Path(args.root).resolve())
    rows = load_pichia_gene_id_standardization(force_refresh=args.force_refresh, paths=paths)
    summary = pichia_gene_id_standardization_summary(paths=paths)
    print(
        json.dumps(
            {
                "cache_path": str(pichia_gene_id_standardization_cache_path(paths)),
                "row_count": len(rows),
                **summary,
            },
            ensure_ascii=False,
            indent=2,
        )
    )
    return 0


def _parse_args(argv: list[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build Pichia gene_id standard naming cache.")
    parser.add_argument("--root", default=str(REPO_ROOT))
    parser.add_argument("--force-refresh", action="store_true")
    return parser.parse_args(argv)


if __name__ == "__main__":
    raise SystemExit(main())
