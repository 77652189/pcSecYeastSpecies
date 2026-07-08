from __future__ import annotations

import csv
import json
from pathlib import Path

from pcsec_pichia.core.paths import ProjectPaths

from app.services.screen_report_fact_pack import build_screen_report_fact_pack


def _paths(tmp_path: Path) -> ProjectPaths:
    for directory in ("local_runs", "Results", "Data", "Model", "Enzymedata"):
        (tmp_path / directory).mkdir(parents=True, exist_ok=True)
    return ProjectPaths(tmp_path)


def _write_run(tmp_path: Path, run_name: str, rows: list[dict[str, object]]) -> Path:
    run_dir = tmp_path / "local_runs" / run_name
    run_dir.mkdir(parents=True)
    csv_path = run_dir / "gene_tradeoff_rows.csv"
    fieldnames = sorted({key for row in rows for key in row})
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    (run_dir / "status.json").write_text(
        json.dumps({"status": "done", "targets": ["hLF", "OPN_ALPHA_FULL_PROJECT"], "csv_path": str(csv_path)}),
        encoding="utf-8",
    )
    return csv_path


def _row(target_id: str, gene_id: str, intervention: str, ratio: float, *, kind: str = "gene") -> dict[str, object]:
    return {
        "target_id": target_id,
        "gene_id": gene_id,
        "candidate_kind": kind,
        "intervention_type": intervention,
        "support_status": "ko_runnable_gpr_gene_deletion" if intervention == "KO" else "reaction_level_proxy",
        "secretion_ratio_vs_wildtype": ratio,
        "growth_retention_ratio": 1.0,
        "max_feasible_mu": 0.1,
        "standard_symbol": gene_id.replace("PAS_", "SYM_"),
        "protein_name": f"Protein {gene_id}",
    }


def test_fact_pack_groups_hlf_and_opn_without_target_mixup(tmp_path: Path) -> None:
    paths = _paths(tmp_path)
    csv_path = _write_run(
        tmp_path,
        "fixture_run",
        [
            _row("hLF", "PAS_HLF_KO", "KO", 1.2),
            _row("OPN_ALPHA_FULL_PROJECT", "PAS_OPN_OE", "OE", 1.3),
        ],
    )

    fact_pack = build_screen_report_fact_pack(paths, csv_paths=(csv_path,))

    assert fact_pack["schema_version"] == 1
    assert fact_pack["targets"]["hLF"]["candidate_count"] == 1
    assert fact_pack["targets"]["OPN"]["candidate_count"] == 1
    evidence_ids = [item["evidence_id"] for item in fact_pack["evidence_items"]]
    assert evidence_ids == ["hLF-KO-0001", "OPN-OE-0001"]
    assert fact_pack["targets"]["hLF"]["useful_ko_candidates"][0]["gene_id"] == "PAS_HLF_KO"
    assert fact_pack["targets"]["OPN"]["useful_oe_candidates"][0]["gene_id"] == "PAS_OPN_OE"


def test_fact_pack_keeps_homology_auxiliary_out_of_executable_recommendations(tmp_path: Path) -> None:
    paths = _paths(tmp_path)
    csv_path = _write_run(
        tmp_path,
        "fixture_run",
        [
            _row("hLF", "PAS_MODEL_EXTERNAL", "KO", 1.5, kind="homology_auxiliary"),
        ],
    )

    fact_pack = build_screen_report_fact_pack(paths, csv_paths=(csv_path,))

    assert fact_pack["targets"]["hLF"]["useful_ko_candidates"] == []
    assert fact_pack["targets"]["hLF"]["manual_review_candidates"][0]["gene_id"] == "PAS_MODEL_EXTERNAL"
    assert fact_pack["evidence_items"][0]["numeric_fields"]["secretion_ratio_vs_wildtype"] == 1.5
    assert any("Results/" in warning for warning in fact_pack["warnings"])
