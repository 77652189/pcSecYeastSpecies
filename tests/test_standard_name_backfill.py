from __future__ import annotations

import csv
import json
from pathlib import Path

import scripts.backfill_pichia_standard_gene_names as backfill


def _fake_standard_rows(*, paths=None):
    return [
        {
            "gene_id": "PAS_chr2-1_0140",
            "display_name": "KAR2",
            "standard_symbol": "KAR2",
            "protein_name": "BiP molecular chaperone",
            "external_ids": {"uniprot": "C4R8K4"},
            "annotation_sources": ["UniProt"],
            "annotation_confidence": "high_exact_locus_tag",
        }
    ]


def test_standard_name_backfill_dry_run_apply_idempotent_and_protects_results(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(backfill, "REPO_ROOT", tmp_path)
    monkeypatch.setattr(backfill, "REPORT_PATH", tmp_path / "local_runs" / "standard_name_backfill_report.json")
    monkeypatch.setattr(backfill, "load_pichia_gene_id_standardization", _fake_standard_rows)
    local_runs = tmp_path / "local_runs"
    local_runs.mkdir()
    results = tmp_path / "Results"
    results.mkdir()
    (results / "historical.csv").write_text("gene_id,objective_value\nPAS_chr2-1_0140,9\n", encoding="utf-8")

    csv_path = local_runs / "gene_tradeoff_rows.csv"
    csv_path.write_text(
        "\n".join(
            [
                "target_id,gene_id,candidate_kind,intervention_type,objective_value,secretion_ratio_vs_wildtype",
                "hLF,PAS_chr2-1_0140,gene,KO,1.5,1.2",
                "hLF,sec_Kar2p_complex_formation,catalog_reaction,OE,1.0,1.0",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    json_path = local_runs / "hLF_summary.json"
    json_payload = {
        "screen_results": [
            {
                "rows": [
                    {
                        "target_id": "hLF",
                        "gene_id": "PAS_chr2-1_0140",
                        "candidate_kind": "gene",
                        "intervention_type": "KO",
                        "objective_value": 1.5,
                    }
                ]
            }
        ]
    }
    json_path.write_text(json.dumps(json_payload), encoding="utf-8")
    before_csv = csv_path.read_text(encoding="utf-8")
    before_json = json_path.read_text(encoding="utf-8")

    dry_report = backfill.run_backfill(dry_run=True)

    assert csv_path.read_text(encoding="utf-8") == before_csv
    assert json_path.read_text(encoding="utf-8") == before_json
    assert not backfill.REPORT_PATH.exists()
    assert dry_report["protected_results_detected"] is True
    assert dry_report["updated_file_count"] == 2

    apply_report = backfill.run_backfill(dry_run=False)
    second_apply_report = backfill.run_backfill(dry_run=False)

    with csv_path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    updated_json = json.loads(json_path.read_text(encoding="utf-8"))
    json_row = updated_json["screen_results"][0]["rows"][0]

    assert rows[0]["gene_display_name"] == "KAR2"
    assert rows[0]["standard_symbol"] == "KAR2"
    assert rows[0]["external_ids"] == '{"uniprot": "C4R8K4"}'
    assert rows[0]["objective_value"] == "1.5"
    assert rows[0]["secretion_ratio_vs_wildtype"] == "1.2"
    assert rows[1]["standard_name_status"] == "not_gene_candidate"
    assert json_row["standard_symbol"] == "KAR2"
    assert json_row["objective_value"] == 1.5
    assert apply_report["numeric_invariance_status"] == "passed"
    assert second_apply_report["updated_file_count"] == 0
    assert backfill.REPORT_PATH.exists()
