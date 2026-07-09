from __future__ import annotations

import json

import pytest

from pcsec_pichia.external_refs import (
    CACHE_SCHEMA_VERSION,
    ExternalGeneFunctionEvidence,
    ExternalGprCandidateEvidence,
    ExternalReactionAssociation,
    ExternalReferenceCacheManifest,
    ExternalReferenceProvenance,
    ExternalReferenceRecord,
    ExternalReferenceSchemaError,
    build_external_reference_manifest,
    load_external_reference_cache,
    load_external_reference_manifest,
    record_from_dict,
    record_to_dict,
    validate_external_reference_cache,
    validate_no_duplicate_cache_keys,
    write_external_reference_cache,
    write_external_reference_cache_bundle,
)


def test_external_reference_schema_roundtrips_all_round1_record_types(tmp_path) -> None:
    records = (
        _reference_record("P12345"),
        _gene_function("PAS_chr1-1_0001"),
        _reaction_association("YBR160W"),
        _gpr_candidate("YBR160W and YDR123C"),
    )
    path = tmp_path / "external_reference_records.jsonl"

    write_external_reference_cache(records, path)
    loaded = load_external_reference_cache(path)

    assert loaded == records
    assert [record.record_type for record in loaded] == [
        "external_reference",
        "gene_function",
        "reaction_association",
        "gpr_candidate",
    ]
    assert all(record.cache_key for record in loaded)


def test_external_reference_manifest_and_bundle_validate_cache_counts(tmp_path) -> None:
    records = (_reference_record("P12345"), _gene_function("PAS_chr1-1_0001"))

    manifest = write_external_reference_cache_bundle(
        records,
        tmp_path,
        query_count=3,
        failed_query_count=1,
        input_cache_fingerprint="homology-cache-fingerprint",
        warnings=("one query failed",),
    )
    loaded_manifest = load_external_reference_manifest(tmp_path / "external_reference_manifest.json")
    validated_manifest = validate_external_reference_cache(
        tmp_path / "external_reference_records.jsonl",
        manifest_path=tmp_path / "external_reference_manifest.json",
    )

    assert manifest == loaded_manifest
    assert manifest.cache_schema_version == CACHE_SCHEMA_VERSION
    assert manifest.query_count == 3
    assert manifest.record_count == 2
    assert manifest.failed_query_count == 1
    assert manifest.source_counts == {"uniprot": 2}
    assert manifest.record_type_counts == {"external_reference": 1, "gene_function": 1}
    assert validated_manifest.record_count == 2
    assert validated_manifest.source_counts == {"uniprot": 2}


def test_external_reference_cache_rejects_duplicate_cache_keys(tmp_path) -> None:
    records = (_reference_record("P12345"), _reference_record("P12345"))

    with pytest.raises(ExternalReferenceSchemaError, match="Duplicate external reference cache key"):
        write_external_reference_cache(records, tmp_path / "records.jsonl")

    with pytest.raises(ExternalReferenceSchemaError, match="Duplicate external reference cache key"):
        validate_no_duplicate_cache_keys(records)


def test_external_reference_cache_rejects_mismatched_manifest_before_writing(tmp_path) -> None:
    records = (_reference_record("P12345"),)
    manifest = ExternalReferenceCacheManifest(
        generated_at="2026-07-09T00:00:00Z",
        cache_schema_version=CACHE_SCHEMA_VERSION,
        query_count=2,
        record_count=2,
        failed_query_count=0,
        source_counts={"uniprot": 2},
        record_type_counts={"external_reference": 2},
    )
    path = tmp_path / "records.jsonl"

    with pytest.raises(ExternalReferenceSchemaError, match="does not match JSONL record count"):
        write_external_reference_cache(records, path, manifest=manifest)

    assert not path.exists()


def test_external_reference_schema_requires_full_provenance() -> None:
    payload = record_to_dict(_reference_record("P12345"))
    payload["provenance"]["source_url"] = ""

    with pytest.raises(ExternalReferenceSchemaError, match="Missing provenance"):
        record_from_dict(payload)


def test_external_reference_cache_rejects_invalid_jsonl(tmp_path) -> None:
    path = tmp_path / "records.jsonl"
    path.write_text("{not-json}\n", encoding="utf-8")

    with pytest.raises(ExternalReferenceSchemaError, match="Invalid JSONL"):
        load_external_reference_cache(path)


def test_external_gpr_candidate_cannot_claim_model_executable_without_current_model_mapping() -> None:
    candidate = _gpr_candidate(
        "YBR160W",
        candidate_status="model_gpr_executable",
        mapped_pichia_reaction_id=None,
        mapped_pichia_gene_ids=(),
    )

    with pytest.raises(ExternalReferenceSchemaError, match="model_gpr_executable requires"):
        record_to_dict(candidate)


def test_external_reference_manifest_rejects_schema_version_drift(tmp_path) -> None:
    manifest = build_external_reference_manifest((_reference_record("P12345"),))
    payload = {
        **json.loads(json.dumps(manifest.__dict__, default=list)),
        "cache_schema_version": "external_refs.v999",
    }

    with pytest.raises(ExternalReferenceSchemaError, match="Unsupported cache schema version"):
        load_external_reference_manifest(_write_json(tmp_path / "manifest.json", payload))


def test_external_reference_manifest_rejects_count_drift(tmp_path) -> None:
    payload = {
        "generated_at": "2026-07-09T00:00:00Z",
        "cache_schema_version": CACHE_SCHEMA_VERSION,
        "query_count": 1,
        "record_count": 2,
        "failed_query_count": 0,
        "source_counts": {"uniprot": 1},
        "record_type_counts": {"external_reference": 2},
        "duplicate_key_count": 0,
        "warnings": [],
    }

    with pytest.raises(ExternalReferenceSchemaError, match="source_counts total"):
        load_external_reference_manifest(_write_json(tmp_path / "manifest.json", payload))


def _provenance(query: str = "P12345") -> ExternalReferenceProvenance:
    return ExternalReferenceProvenance(
        source_database="uniprot",
        source_version="2026_03",
        source_url=f"https://rest.uniprot.org/uniprotkb/{query}.json",
        source_query=query,
        retrieved_at="2026-07-09T00:00:00Z",
        raw_record_sha256="a" * 64,
    )


def _reference_record(accession: str) -> ExternalReferenceRecord:
    return ExternalReferenceRecord(
        provenance=_provenance(accession),
        taxon_id="4922",
        organism="Komagataella phaffii",
        primary_accession=accession,
        gene_id="PAS_chr1-1_0001",
        gene_name="SEC1",
        locus_tag="PAS_chr1-1_0001",
        aliases=("SEC1", "secretory protein 1"),
        protein_name="Secretory pathway protein",
        reviewed=True,
    )


def _gene_function(gene_id: str) -> ExternalGeneFunctionEvidence:
    return ExternalGeneFunctionEvidence(
        provenance=_provenance(gene_id),
        gene_id=gene_id,
        protein_name="Secretory pathway protein",
        function_description="Annotation-only external function evidence.",
        ec_numbers=("1.1.1.1",),
        go_terms=("GO:0006886",),
        pathways=("secretory pathway",),
        reviewed=True,
    )


def _reaction_association(gene_id: str) -> ExternalReactionAssociation:
    return ExternalReactionAssociation(
        provenance=_provenance(gene_id),
        external_model_id="yeast-GEM",
        external_reaction_id="r_1234",
        external_reaction_name="Example secretion reaction",
        external_gene_ids=(gene_id,),
        gene_rule=gene_id,
        ec_numbers=("2.7.1.1",),
    )


def _gpr_candidate(
    gene_rule: str,
    *,
    candidate_status: str = "external_gpr_candidate",
    mapped_pichia_reaction_id: str | None = None,
    mapped_pichia_gene_ids: tuple[str, ...] = (),
) -> ExternalGprCandidateEvidence:
    return ExternalGprCandidateEvidence(
        provenance=_provenance(gene_rule),
        external_model_id="yeast-GEM",
        external_reaction_id="r_5678",
        external_gene_rule=gene_rule,
        candidate_status=candidate_status,
        mapped_pichia_reaction_id=mapped_pichia_reaction_id,
        mapped_pichia_gene_ids=mapped_pichia_gene_ids,
    )


def _write_json(path, payload):
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path
