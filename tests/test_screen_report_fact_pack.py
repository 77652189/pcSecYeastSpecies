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
    assert [item["rank"] for item in fact_pack["evidence_items"]] == [1, 1]
    assert fact_pack["targets"]["hLF"]["useful_ko_candidates"][0]["gene_id"] == "PAS_HLF_KO"
    assert fact_pack["targets"]["OPN"]["useful_oe_candidates"][0]["gene_id"] == "PAS_OPN_OE"


def test_fact_pack_assigns_deterministic_target_local_prediction_ranks(tmp_path: Path) -> None:
    paths = _paths(tmp_path)
    csv_path = _write_run(
        tmp_path,
        "fixture_run",
        [
            _row("hLF", "PAS_HLF_LOW", "KO", 1.1),
            _row("OPN_ALPHA_FULL_PROJECT", "PAS_OPN", "OE", 1.2),
            _row("hLF", "PAS_HLF_HIGH", "KO", 1.4),
        ],
    )

    fact_pack = build_screen_report_fact_pack(paths, csv_paths=(csv_path,))

    ranks = {
        (item["target_id"], item["gene_id"]): item["rank"]
        for item in fact_pack["evidence_items"]
    }
    assert ranks[("hLF", "PAS_HLF_HIGH")] == 1
    assert ranks[("hLF", "PAS_HLF_LOW")] == 2
    assert ranks[("OPN_ALPHA_FULL_PROJECT", "PAS_OPN")] == 1
    assert fact_pack["targets"]["hLF"]["top_candidates"][0]["rank"] == 1


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


def test_fact_pack_preserves_external_model_gpr_fields_without_tier_upgrade(tmp_path: Path) -> None:
    paths = _paths(tmp_path)
    row = _row("hLF", "PAS_EXTERNAL_GPR", "KO", 1.4)
    row.update(
        {
            "recommendation_tier": "model_executable",
            "database_annotation_sources": json.dumps(["UniProt"]),
            "external_model_sources": json.dumps(["Kp.1.0"]),
            "gpr_source_priority": json.dumps({"best_priority_tier": "pichia_literature_model_gpr"}),
            "external_gpr_candidate_count": 2,
            "best_external_gpr_source": "biomodels:Kp.1.0",
            "external_gpr_mapping_status": json.dumps({"gene_mapping_required": 2}),
            "external_gpr_conflict_warnings": json.dumps(["conflicting external GPR rules"]),
            "manual_review_reasons": json.dumps(["external GPR candidate requires mapped current model gene"]),
        }
    )
    csv_path = _write_run(tmp_path, "fixture_run", [row])

    fact_pack = build_screen_report_fact_pack(paths, csv_paths=(csv_path,))

    item = fact_pack["evidence_items"][0]
    brief = fact_pack["targets"]["hLF"]["manual_review_candidates"][0]
    assert item["recommendation_tier"] == "model_executable"
    assert item["recommendation_tier"] != "experiment_calibrated"
    assert item["external_model_sources"] == ["Kp.1.0"]
    assert item["gpr_source_priority"]["best_priority_tier"] == "pichia_literature_model_gpr"
    assert item["external_gpr_candidate_count"] == 2
    assert item["external_gpr_mapping_status"] == {"gene_mapping_required": 2}
    assert item["external_gpr_conflict_warnings"] == ["conflicting external GPR rules"]
    assert brief["best_external_gpr_source"] == "biomodels:Kp.1.0"
    assert brief["manual_review_reasons"] == ["external GPR candidate requires mapped current model gene"]
