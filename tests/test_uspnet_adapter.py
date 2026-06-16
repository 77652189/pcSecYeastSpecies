from __future__ import annotations

from pathlib import Path

from app.adapters.uspnet import USPNetAdapter, parse_uspnet_results


def test_uspnet_parser_maps_results_to_fasta_order(tmp_path: Path) -> None:
    results = tmp_path / "results.csv"
    results.write_text(
        "sequence,predicted_type,predicted_cleavage\n"
        "MKALLLALLALAAASAGAQREST,SP,MKALLLALLALAAASAGA\n"
        "MNNNNNNNNNNNNNNNNNNNN,NO_SP,\n",
        encoding="utf-8",
    )

    predictions = parse_uspnet_results(results, ["OPN_UNIPROT_X12345", "OPN_UNIPROT_Y12345"])

    assert predictions[0].candidate_id == "OPN_UNIPROT_X12345"
    assert predictions[0].predicted_type == "SP"
    assert predictions[0].passed is True
    assert predictions[1].candidate_id == "OPN_UNIPROT_Y12345"
    assert predictions[1].passed is False


def test_uspnet_adapter_reports_missing_repo(tmp_path: Path) -> None:
    adapter = USPNetAdapter(repo_dir=tmp_path / "missing-USPNet")

    status = adapter.status()

    assert status.available is False
    assert "未检测到 USPNet 本地仓库" in status.message
    assert "MIT License" in status.message
