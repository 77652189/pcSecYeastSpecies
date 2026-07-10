from __future__ import annotations

import json

from pcsec_pichia.external_refs import (
    ExternalModelInventoryRecord,
    default_external_model_inventory_records,
    load_external_model_inventory,
    write_external_model_inventory,
)


def test_default_external_model_inventory_covers_required_round_a_resources() -> None:
    records = default_external_model_inventory_records()
    by_id = {record.model_id: record for record in records}

    assert {
        "iPichia",
        "ecPichia",
        "Kp.1.0",
        "iAUKM",
        "Yeast8_Yeast9",
        "BioModels_Kp.1.0_MODEL1703150000",
        "GPRuler",
    }.issubset(by_id)
    assert by_id["Kp.1.0"].has_sbml is True
    assert by_id["Kp.1.0"].has_gpr is True
    assert "SBML" in by_id["Kp.1.0"].available_artifact_types
    assert by_id["iPichia"].download_status == "needs_manual_access"
    assert by_id["ecPichia"].download_status == "needs_manual_access"
    assert by_id["GPRuler"].download_status == "tool_only_not_primary_gem"
    assert "supplemental_gpr_tool_only" in by_id["GPRuler"].warnings


def test_inventory_records_require_source_or_warning_and_no_fake_local_artifacts() -> None:
    for record in default_external_model_inventory_records():
        record.validate()
        assert record.source_url or record.publication_url
        if record.download_status in {"needs_manual_access", "publication_only", "tool_only_not_primary_gem"}:
            assert record.local_path == ""
            assert record.checksum_sha256 == ""


def test_write_external_model_inventory_outputs_jsonl_tsv_and_report(tmp_path) -> None:
    records = (
        ExternalModelInventoryRecord(
            model_id="toy_model",
            model_name="Toy GEM",
            organism="Komagataella phaffii",
            source_database_or_repository="test repository",
            source_url="https://example.test/toy",
            publication_url="https://example.test/paper",
            license="unknown",
            available_artifact_types=("SBML",),
            download_status="downloadable",
            local_path="",
            checksum_sha256="",
            has_gpr=True,
            has_gene_ids=True,
            has_reaction_ids=True,
            has_sbml=True,
            notes="toy inventory record",
            warnings=("not_downloaded_in_round_a",),
        ),
    )

    outputs = write_external_model_inventory(records, tmp_path)

    assert outputs.jsonl_path.name == "external_model_inventory.jsonl"
    assert outputs.tsv_path.name == "external_model_inventory.tsv"
    assert outputs.report_path.name == "external_model_inventory_report.md"
    payload = json.loads(outputs.jsonl_path.read_text(encoding="utf-8").splitlines()[0])
    assert payload["model_id"] == "toy_model"
    assert payload["has_gpr"] is True
    assert load_external_model_inventory(tmp_path) == records
    assert load_external_model_inventory(outputs.jsonl_path) == records
    assert "toy_model" in outputs.tsv_path.read_text(encoding="utf-8")
    report = outputs.report_path.read_text(encoding="utf-8")
    assert "External GEM / GPR Resource Inventory" in report
    assert "does not import external GPR rules into the current Pichia GEM" in report
