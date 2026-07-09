from __future__ import annotations

from pcsec_pichia.external_refs import (
    EXTERNAL_ALIAS_CONFIRMED,
    EXTERNAL_CONFLICT,
    EXTERNAL_LOCUS_CONFIRMED,
    EXTERNAL_MATCH_CONFIRMED,
    EXTERNAL_REFERENCE_MISSING,
    ExternalFetchResult,
    ExternalReferenceQuery,
    attach_external_references_to_name_audit,
    merge_external_fetch_results,
)
from pcsec_pichia.external_refs.schema import ExternalReferenceProvenance, ExternalReferenceRecord
from pcsec_pichia.homology.cache_schema import NameAuditRow


def test_merge_external_fetch_results_dedupes_successful_records_only() -> None:
    query = _query("P12345")
    record = _record("uniprot", "P12345", gene_name="SEC1", locus_tag="PAS_chr1-1_0001")
    records = merge_external_fetch_results(
        (
            ExternalFetchResult(source_database="uniprot", query=query, success=True, records=(record,)),
            ExternalFetchResult(source_database="uniprot", query=query, success=True, records=(record,)),
            ExternalFetchResult(source_database="sgd", query=query, success=False),
        )
    )

    assert records == (record,)


def test_attach_external_references_to_name_audit_preserves_internal_model_facts() -> None:
    rows = (
        _name_row("PAS_chr1-1_0001", "SEC1", "P12345", "SEC1", "PAS_chr1-1_0001", is_rbh=True),
        _name_row("PAS_chr1-1_0002", "DOA10", "P22222", "SSM4", "PAS_chr1-1_0002", is_rbh=False),
        _name_row("PAS_chr1-1_0003", "VPS1", "P33333", "VPSX", "PAS_chr1-1_0003", in_model=False),
        _name_row("PAS_chr1-1_0004", "HRD1", "P44444", "HRD1", "PAS_chr1-1_0004"),
        _name_row("PAS_chr1-1_0005", "KAR2", "P55555", "KAR2", "PAS_chr1-1_0005"),
    )
    records = (
        _record("uniprot", "P12345", gene_name="SEC1", locus_tag="PAS_chr1-1_0001"),
        _record("uniprot", "P22222", gene_name="SSM4", aliases=("DOA10",)),
        _record("ncbi", "333", gene_name="uncharacterized", locus_tag="PAS_chr1-1_0003"),
        _record("uniprot", "P44444", gene_name="PEP4", locus_tag="PAS_chr9_9999", warnings=("name conflict",)),
    )

    merged = attach_external_references_to_name_audit(rows, records)
    by_gene = {row["internal_gene_id"]: row for row in merged}

    assert by_gene["PAS_chr1-1_0001"]["external_name_status"] == EXTERNAL_MATCH_CONFIRMED
    assert by_gene["PAS_chr1-1_0002"]["external_name_status"] == EXTERNAL_ALIAS_CONFIRMED
    assert by_gene["PAS_chr1-1_0003"]["external_name_status"] == EXTERNAL_LOCUS_CONFIRMED
    assert by_gene["PAS_chr1-1_0004"]["external_name_status"] == EXTERNAL_CONFLICT
    assert by_gene["PAS_chr1-1_0005"]["external_name_status"] == EXTERNAL_REFERENCE_MISSING

    assert by_gene["PAS_chr1-1_0002"]["is_rbh"] is False
    assert by_gene["PAS_chr1-1_0003"]["in_model_gene_index"] is False
    assert by_gene["PAS_chr1-1_0004"]["review_status"] == "model_ready_rbh_high_confidence"
    assert "recommendation_tier" not in by_gene["PAS_chr1-1_0001"]
    assert by_gene["PAS_chr1-1_0004"]["external_reference_warnings"] == ("name conflict",)
    assert by_gene["PAS_chr1-1_0005"]["external_manual_review_reasons"]


def test_attach_external_references_to_mapping_rows_does_not_overwrite_existing_recommendation_tier() -> None:
    rows = (
        {
            "internal_gene_id": "PAS_chr1-1_0001",
            "internal_common_name": "SEC1",
            "external_accession": "P12345",
            "external_gene_name": "SEC1",
            "external_locus_tag": "PAS_chr1-1_0001",
            "recommendation_tier": "model_executable",
        },
    )
    merged = attach_external_references_to_name_audit(
        rows,
        (_record("uniprot", "P12345", gene_name="SEC1", locus_tag="PAS_chr1-1_0001"),),
    )

    assert merged[0]["external_name_status"] == EXTERNAL_MATCH_CONFIRMED
    assert merged[0]["recommendation_tier"] == "model_executable"


def _query(value: str) -> ExternalReferenceQuery:
    return ExternalReferenceQuery(
        query_type="external_accession",
        query_value=value,
        source_context="test",
        source_id="test",
    )


def _name_row(
    gene_id: str,
    common_name: str,
    accession: str,
    external_gene_name: str,
    locus_tag: str,
    *,
    is_rbh: bool = True,
    in_model: bool = True,
) -> NameAuditRow:
    return NameAuditRow(
        internal_gene_id=gene_id,
        internal_common_name=common_name,
        internal_sequence_id="YDR164C",
        external_accession=accession,
        external_gene_name=external_gene_name,
        external_locus_tag=locus_tag,
        external_aliases=(),
        identity_pct=80.0,
        query_coverage=95.0,
        subject_coverage=94.0,
        evalue=1e-40,
        is_rbh=is_rbh,
        in_model_gene_index=in_model,
        name_consistency_status="name_confirmed_by_rbh",
        review_status="model_ready_rbh_high_confidence",
    )


def _record(
    source: str,
    accession: str,
    *,
    gene_name: str,
    locus_tag: str | None = None,
    aliases: tuple[str, ...] = (),
    warnings: tuple[str, ...] = (),
) -> ExternalReferenceRecord:
    return ExternalReferenceRecord(
        provenance=ExternalReferenceProvenance(
            source_database=source,
            source_version="test",
            source_url=f"https://example.test/{accession}",
            source_query=accession,
            retrieved_at="2026-07-09T00:00:00Z",
            raw_record_sha256="b" * 64,
            warnings=warnings,
        ),
        taxon_id="4922",
        organism="Komagataella phaffii",
        primary_accession=accession,
        gene_id=locus_tag,
        gene_name=gene_name,
        locus_tag=locus_tag,
        aliases=aliases,
    )
