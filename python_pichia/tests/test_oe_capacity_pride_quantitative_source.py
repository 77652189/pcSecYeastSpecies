from __future__ import annotations

import hashlib
import json
from dataclasses import replace
from pathlib import Path

import pytest

from pcsec_pichia.external_refs.capacity_sources import (
    PXD055501_G6PDH2_PROFILE,
    ExternalCapacitySourceType,
    RetrievalMode,
    cache_pride_maxquant_source,
    cache_uniprot_identity_source,
)
from pcsec_pichia.external_refs.clients import ExternalHttpResponse
from pcsec_pichia.oe_capacity.external_candidate_audit import (
    ExternalCapacityAuditRequest,
    run_external_capacity_candidate_audit,
)
from pcsec_pichia.oe_capacity.external_candidate_io import (
    load_external_capacity_candidate_bundle,
    parse_pride_maxquant_g6pdh2_evidence,
)
from pcsec_pichia.oe_capacity.external_candidate_schema import (
    CapacityCandidateStatus,
    CapacityParameterKind,
)
from pcsec_pichia.oe_capacity.schema import OECapacityValidationError


REPO_ROOT = Path(__file__).resolve().parents[2]
SEQUENCE = "MSTNPKPQRKTKRNTNRRPQDVKFPGGGQIVGGV"


def _pride_responses(*, license_text: str = "Creative Commons Public Domain (CC0)", raw_value: str = "12868000"):
    profile = PXD055501_G6PDH2_PROFILE
    project_url = f"https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{profile.project_accession}"
    files_url = project_url + "/files"
    protein_groups = (
        "Protein IDs\tMajority protein IDs"
        + "".join(f"\tiBAQ {sample_id}" for sample_id in profile.sample_ids)
        + "\n"
        + "F2QTE5\tF2QTE5\t"
        + "\t".join((raw_value, "10476000", "21552000"))
        + "\n"
    )
    sample_map = "Raw file\tprotein groups.txt\n" + "".join(
        f"raw-{sample_id}\t{sample_id}\n" for sample_id in profile.sample_ids
    )
    database_fasta = f">sp|F2QTE5|ZWF1_KOMPG Glucose-6-phosphate dehydrogenase\n{SEQUENCE}\n"
    query_fasta = f">tr|C4R099|C4R099_KOMPG G6PDH2\n{SEQUENCE}\n"
    file_text = {
        profile.protein_groups_filename: protein_groups,
        profile.sample_map_filename: sample_map,
        profile.database_fasta_filename: database_fasta,
    }
    entries = []
    responses: dict[str, str] = {
        project_url: json.dumps(
            {
                "accession": profile.project_accession,
                "title": "Quantitative proteomics of oxygen limitation in K. phaffii",
                "publicationDate": "2025-03-01T00:00:00Z",
                "license": license_text,
            }
        ),
        profile.query_fasta_url: query_fasta,
    }
    for index, (filename, text) in enumerate(file_text.items(), start=1):
        ftp_url = f"ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2025/03/PXD055501/{filename}"
        https_url = ftp_url.replace("ftp://ftp.pride.ebi.ac.uk/", "https://ftp.pride.ebi.ac.uk/")
        responses[https_url] = text
        entries.append(
            {
                "accession": f"PRIDE_FILE_{index}",
                "fileName": filename,
                "fileSizeBytes": len(text.encode("utf-8")),
                "publicFileLocations": [{"value": ftp_url}],
            }
        )
    responses[files_url] = json.dumps(entries)
    return responses


def _fake_get(responses: dict[str, str]):
    def get(url, config):
        if url not in responses:
            raise AssertionError(f"unexpected external URL: {url}")
        return ExternalHttpResponse(status_code=200, text=responses[url], url=url)

    return get


def _cache_identity(output_dir: Path) -> None:
    raw = json.dumps(
        {
            "results": [
                {
                    "primaryAccession": "C4R099",
                    "uniProtkbId": "A0A1E4RTV1_PICPA",
                    "entryType": "UniProtKB unreviewed (TrEMBL)",
                    "genes": [
                        {
                            "geneName": {"value": "G6PDH2"},
                            "orderedLocusNames": [{"value": "PAS_chr2-1_0308"}],
                        }
                    ],
                    "organism": {
                        "taxonId": 644223,
                        "scientificName": "Komagataella phaffii",
                    },
                    "proteinDescription": {
                        "recommendedName": {
                            "fullName": {
                                "value": "Glucose-6-phosphate 1-dehydrogenase"
                            }
                        }
                    },
                }
            ]
        }
    )
    cache_uniprot_identity_source(
        "PAS_chr2-1_0308",
        output_dir,
        http_get=lambda url, config: ExternalHttpResponse(
            status_code=200,
            text=raw,
            url=url,
            headers={"X-UniProt-Release": "2026_02"},
        ),
    )


def test_pride_adapter_parses_real_shaped_ibaq_and_replays_offline(tmp_path: Path) -> None:
    responses = _pride_responses()
    source = cache_pride_maxquant_source(
        PXD055501_G6PDH2_PROFILE,
        tmp_path,
        http_get=_fake_get(responses),
    )
    evidence = parse_pride_maxquant_g6pdh2_evidence(source)

    assert source.source_id == "pride:PXD055501"
    assert source.license_id == "CC0-1.0"
    assert source.terms_reviewed is True
    assert evidence.measurement.parameter_kind is CapacityParameterKind.ABUNDANCE
    assert evidence.measurement.unit == "ibaq_intensity"
    assert evidence.raw_values == (12868000.0, 10476000.0, 21552000.0)
    assert evidence.measurement.nominal_value == 12868000.0
    assert evidence.measurement.condition.growth_rate_per_h == pytest.approx(0.075)
    assert f"exact_sequence_identity:{len(SEQUENCE)}/{len(SEQUENCE)}" in evidence.mapping_evidence
    assert evidence.quantitative_boundary["absolute_abundance_available"] is False
    assert evidence.artifact_sha256s["protein_groups"]

    replayed = cache_pride_maxquant_source(
        PXD055501_G6PDH2_PROFILE,
        tmp_path,
        http_get=lambda *_: pytest.fail("offline replay must not call the network"),
        retrieval_mode=RetrievalMode.OFFLINE_REPLAY,
    )
    assert replayed.retrieval_mode is RetrievalMode.OFFLINE_REPLAY
    assert parse_pride_maxquant_g6pdh2_evidence(replayed).raw_values == evidence.raw_values


def test_pride_bundle_tamper_is_rejected_before_parsing(tmp_path: Path) -> None:
    source = cache_pride_maxquant_source(
        PXD055501_G6PDH2_PROFILE,
        tmp_path,
        http_get=_fake_get(_pride_responses()),
    )
    bundle = json.loads(Path(source.cache_path).read_text(encoding="utf-8"))
    protein_path = Path(source.cache_path).parent / bundle["artifacts"]["protein_groups"]["filename"]
    protein_path.write_text(protein_path.read_text(encoding="utf-8") + "tamper\n", encoding="utf-8")

    with pytest.raises(OECapacityValidationError, match="artifact sha256 mismatch"):
        parse_pride_maxquant_g6pdh2_evidence(source)


def test_pride_offline_replay_rejects_forged_license_metadata(tmp_path: Path) -> None:
    cache_pride_maxquant_source(
        PXD055501_G6PDH2_PROFILE,
        tmp_path,
        http_get=_fake_get(_pride_responses()),
    )
    metadata_path = tmp_path / "pride-PXD055501.source.json"
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    metadata["license_id"] = "CC-BY-4.0"
    metadata_path.write_text(json.dumps(metadata), encoding="utf-8")

    with pytest.raises(OECapacityValidationError, match="license_id"):
        cache_pride_maxquant_source(
            PXD055501_G6PDH2_PROFILE,
            tmp_path,
            retrieval_mode=RetrievalMode.OFFLINE_REPLAY,
        )


@pytest.mark.parametrize(
    ("profile", "message"),
    (
        (replace(PXD055501_G6PDH2_PROFILE, metric_name=""), "metric_name"),
        (replace(PXD055501_G6PDH2_PROFILE, growth_rate_per_h=0), "condition values"),
    ),
)
def test_pride_profile_requires_unit_and_condition(profile, message: str, tmp_path: Path) -> None:
    with pytest.raises(OECapacityValidationError, match=message):
        cache_pride_maxquant_source(
            profile,
            tmp_path,
            http_get=_fake_get(_pride_responses()),
        )


def test_pride_selected_raw_values_must_all_be_positive(tmp_path: Path) -> None:
    source = cache_pride_maxquant_source(
        PXD055501_G6PDH2_PROFILE,
        tmp_path,
        http_get=_fake_get(_pride_responses(raw_value="0")),
    )
    with pytest.raises(OECapacityValidationError, match="must be positive"):
        parse_pride_maxquant_g6pdh2_evidence(source)


def test_uniprot_identity_cannot_be_parsed_as_quantitative_capacity(tmp_path: Path) -> None:
    _cache_identity(tmp_path)
    source = cache_uniprot_identity_source(
        "PAS_chr2-1_0308",
        tmp_path,
        retrieval_mode=RetrievalMode.OFFLINE_REPLAY,
    )
    assert source.source_type is ExternalCapacitySourceType.IDENTITY_REFERENCE
    with pytest.raises(OECapacityValidationError, match="quantitative_proteomics"):
        parse_pride_maxquant_g6pdh2_evidence(source)


def test_pride_relative_evidence_audit_never_becomes_promotion_ready(tmp_path: Path) -> None:
    identity_cache = tmp_path / "identity"
    quantitative_cache = tmp_path / "quantitative"
    output_dir = tmp_path / "local_runs" / "oe_capacity" / "round6a" / "a0b"
    _cache_identity(identity_cache)
    cache_pride_maxquant_source(
        PXD055501_G6PDH2_PROFILE,
        quantitative_cache,
        http_get=_fake_get(_pride_responses()),
    )

    outputs = run_external_capacity_candidate_audit(
        ExternalCapacityAuditRequest(
            repo_root=REPO_ROOT,
            output_dir=output_dir,
            offline_replay=True,
            identity_cache_dir=identity_cache,
            quantitative_cache_dir=quantitative_cache,
            pride_pxd055501=True,
        )
    )
    bundle = load_external_capacity_candidate_bundle(outputs.candidate_manifest_path)
    candidate = bundle.candidates[0]
    audit = json.loads(outputs.audit_json_path.read_text(encoding="utf-8"))
    report = outputs.audit_markdown_path.read_text(encoding="utf-8")

    assert outputs.candidate_count == 1
    assert outputs.promotion_ready_count == 0
    assert candidate.status is CapacityCandidateStatus.REVIEW_REQUIRED
    assert candidate.nominal_capacity is None
    assert candidate.conversion_steps[0].step_id == "retain-raw-abundance-evidence"
    assert "absolute_abundance_calibration" in candidate.conversion_steps[0].missing_metadata
    assert "baseline_capacity_or_abundance_and_kcat" in candidate.missing_information
    assert "formal_glucose_mu_0.1_condition_match" in candidate.missing_information
    assert audit["source_assessments"][0]["raw_values"] == [12868000.0, 10476000.0, 21552000.0]
    assert audit["source_assessments"][0]["formal_context_match"] is False
    assert "10.1111/1751-7915.70106" in audit["source_assessments"][0]["condition_source_ref"]
    assert audit["source_assessments"][0]["artifact_sha256s"]["protein_groups"]
    quantitative_row = next(
        row for row in audit["sources_checked"]
        if row.get("source") == "same-host quantitative Pichia proteomics"
    )
    assert quantitative_row["quantitative_value_available"] is True
    assert quantitative_row["capacity_value_available"] is False
    assert "relative quantitative evidence only" in report
    assert "Unit chain: iBAQ intensity" in report
    assert output_dir.is_relative_to(tmp_path / "local_runs")


def test_unreviewed_license_does_not_upgrade_relative_evidence(tmp_path: Path) -> None:
    source = cache_pride_maxquant_source(
        PXD055501_G6PDH2_PROFILE,
        tmp_path,
        http_get=_fake_get(_pride_responses(license_text="terms not declared")),
    )
    evidence = parse_pride_maxquant_g6pdh2_evidence(source)

    assert source.terms_reviewed is False
    assert source.license_id == "unreviewed"
    assert "license_review_required" in source.warnings
    assert evidence.measurement.nominal_value > 0
    assert evidence.quantitative_boundary["model_flux_conversion_available"] is False
    assert hashlib.sha256(Path(source.cache_path).read_bytes()).hexdigest() == source.raw_sha256
