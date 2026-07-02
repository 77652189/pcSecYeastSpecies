from __future__ import annotations

import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from app.services.pichia_gene_catalog_service import build_pichia_gene_evidence_cache


def _progress(message: str) -> None:
    print(message, flush=True)


def main() -> None:
    summary = build_pichia_gene_evidence_cache(paths=None, progress=_progress, refresh_full_catalog=False)
    output = Path("local_runs") / "gene_evidence_cache" / "gene_evidence_summary.json"
    print(json.dumps({"summary_path": str(output), **summary}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
