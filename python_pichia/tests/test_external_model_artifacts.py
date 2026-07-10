from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

from pcsec_pichia.external_refs import (
    ExternalModelArtifactRequest,
    ExternalModelInventoryRecord,
    build_artifact_requests_from_inventory,
    cache_external_model_artifacts,
    write_external_model_inventory,
)


def test_cache_external_model_artifacts_downloads_with_checksum(tmp_path) -> None:
    payload = b"<sbml>toy</sbml>\n"
    expected = hashlib.sha256(payload).hexdigest()
    request = ExternalModelArtifactRequest(
        model_id="toy_model",
        artifact_url="https://example.test/toy.xml",
        artifact_type="SBML",
        filename="toy.xml",
        expected_sha256=expected,
    )

    outputs = cache_external_model_artifacts(
        (request,),
        tmp_path,
        fetcher=lambda url, timeout: payload,
    )

    assert outputs.manifest.downloaded_count == 1
    assert outputs.manifest.failed_count == 0
    result = outputs.manifest.results[0]
    assert result.download_status == "downloaded"
    assert result.checksum_sha256 == expected
    assert result.local_path
    assert tuple(result.local_path.replace("\\", "/").split("/")[-3:]) == ("artifacts", "toy_model", "toy.xml")
    assert (tmp_path / "artifacts" / "toy_model" / "toy.xml").read_bytes() == payload


def test_cache_external_model_artifacts_records_manual_required_without_fetching(tmp_path) -> None:
    called = False

    def _fetcher(url: str, timeout: float) -> bytes:
        nonlocal called
        called = True
        return b"unexpected"

    outputs = cache_external_model_artifacts(
        (
            ExternalModelArtifactRequest(
                model_id="iPichia",
                artifact_url="",
                artifact_type="SBML",
                filename="iPichia.xml",
                requires_manual_access=True,
                source_page_url="https://doi.org/10.1016/j.bej.2025.109940",
            ),
        ),
        tmp_path,
        fetcher=_fetcher,
    )

    assert called is False
    assert outputs.manifest.manual_required_count == 1
    assert outputs.manifest.results[0].download_status == "manual_download_required"
    assert outputs.manifest.results[0].local_path == ""
    assert outputs.manifest.results[0].checksum_sha256 == ""


def test_cache_external_model_artifacts_flags_checksum_mismatch(tmp_path) -> None:
    request = ExternalModelArtifactRequest(
        model_id="toy_model",
        artifact_url="https://example.test/toy.xml",
        artifact_type="SBML",
        filename="toy.xml",
        expected_sha256="0" * 64,
    )

    outputs = cache_external_model_artifacts(
        (request,),
        tmp_path,
        fetcher=lambda url, timeout: b"different",
    )

    assert outputs.manifest.downloaded_count == 0
    assert outputs.manifest.failed_count == 1
    assert outputs.manifest.checksum_mismatch_count == 1
    assert outputs.manifest.results[0].download_status == "checksum_mismatch"
    assert not (tmp_path / "artifacts" / "toy_model" / "toy.xml").exists()


def test_build_artifact_requests_from_inventory_needs_direct_artifact_url(tmp_path) -> None:
    records = (
        ExternalModelInventoryRecord(
            model_id="Kp.1.0",
            model_name="Kp.1.0",
            organism="Komagataella phaffii",
            source_database_or_repository="repository",
            source_url="https://www.repository.cam.ac.uk/items/02da7483-3966-4d96-b90d-eda1e890e104",
            publication_url="https://doi.org/10.1002/bit.26380",
            license="CC BY 4.0",
            available_artifact_types=("SBML",),
            download_status="downloadable",
            local_path="",
            checksum_sha256="",
            has_gpr=True,
            has_gene_ids=True,
            has_reaction_ids=True,
            has_sbml=True,
            notes="repository page, not a direct artifact file",
            warnings=("not_downloaded_in_round_a",),
        ),
    )

    requests = build_artifact_requests_from_inventory(records)
    assert len(requests) == 1
    assert requests[0].requires_manual_access is True
    assert requests[0].artifact_url == ""
    assert requests[0].source_page_url == records[0].source_url

    outputs = cache_external_model_artifacts(requests, tmp_path)
    manifest_payload = json.loads(outputs.manifest_path.read_text(encoding="utf-8"))
    assert manifest_payload["manual_required_count"] == 1


def test_cache_external_model_artifacts_cli_accepts_inventory_dir(tmp_path) -> None:
    inventory_dir = tmp_path / "inventory"
    output_dir = tmp_path / "artifacts"
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
    write_external_model_inventory(records, inventory_dir)

    repo_root = Path(__file__).resolve().parents[2]
    result = subprocess.run(
        [
            sys.executable,
            str(repo_root / "scripts" / "cache_external_model_artifacts.py"),
            "--inventory-dir",
            str(inventory_dir),
            "--output-dir",
            str(output_dir),
        ],
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr
    payload = json.loads(result.stdout)
    assert payload["request_count"] == 1
    assert payload["manual_required_count"] == 1
    assert (output_dir / "external_model_artifact_manifest.json").exists()
