from __future__ import annotations

import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
PYTHON_PICHIA_SRC = REPO_ROOT / "python_pichia" / "src"
if str(PYTHON_PICHIA_SRC) not in sys.path:
    sys.path.insert(0, str(PYTHON_PICHIA_SRC))

from pcsec_pichia.loading import load_pcsec_pichia_inputs
from pcsec_pichia.services.gene_rule_overlay import (
    DEFAULT_GENE_RULE_EVIDENCE_CACHE,
    DEFAULT_GENE_RULE_EVIDENCE_REPORT,
    DEFAULT_GENE_RULE_EVIDENCE_SUMMARY,
    build_gene_rule_evidence_cache,
)


def _progress(message: str) -> None:
    print(message, flush=True)


def main() -> None:
    inputs = load_pcsec_pichia_inputs(REPO_ROOT)
    summary = build_gene_rule_evidence_cache(
        output_path=REPO_ROOT / DEFAULT_GENE_RULE_EVIDENCE_CACHE,
        summary_path=REPO_ROOT / DEFAULT_GENE_RULE_EVIDENCE_SUMMARY,
        report_path=REPO_ROOT / DEFAULT_GENE_RULE_EVIDENCE_REPORT,
        model=inputs.prepared_model,
        progress=_progress,
    )
    print(
        json.dumps(
            {
                "cache_path": str(REPO_ROOT / DEFAULT_GENE_RULE_EVIDENCE_CACHE),
                "summary_path": str(REPO_ROOT / DEFAULT_GENE_RULE_EVIDENCE_SUMMARY),
                "report_path": str(REPO_ROOT / DEFAULT_GENE_RULE_EVIDENCE_REPORT),
                "total_records": summary.get("total_records"),
                "high_confidence_count": summary.get("high_confidence_count"),
                "executable_overlay_entry_count": summary.get("executable_overlay_entry_count"),
            },
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
