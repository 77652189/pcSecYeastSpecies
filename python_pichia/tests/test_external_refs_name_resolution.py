from __future__ import annotations

from pcsec_pichia.external_refs import (
    EXTERNAL_ALIAS_CONFIRMED,
    EXTERNAL_CONFLICT,
    EXTERNAL_LOCUS_CONFIRMED,
    EXTERNAL_MATCH_CONFIRMED,
    EXTERNAL_REFERENCE_MISSING,
    classify_external_name_consistency,
    select_external_records_for_name_audit_row,
)
from pcsec_pichia.external_refs.schema import ExternalReferenceProvenance, ExternalReferenceRecord


def test_external_name_resolution_classifies_exact_alias_locus_conflict_and_missing() -> None:
    exact = classify_external_name_consistency(
        internal_gene_id="PAS_chr1-1_0001",
        internal_common_name="SEC1",
        internal_aliases=(),
        external_records=(_record("uniprot", "P12345", gene_name="SEC1", locus_tag="PAS_chr1-1_0001"),),
    )
    alias = classify_external_name_consistency(
        internal_gene_id="PAS_chr1-1_0002",
        internal_common_name="DOA10",
        internal_aliases=(),
        external_records=(_record("uniprot", "P22222", gene_name="SSM4", aliases=("DOA10",)),),
    )
    locus = classify_external_name_consistency(
        internal_gene_id="PAS_chr1-1_0003",
        internal_common_name="VPS1",
        internal_aliases=(),
        external_records=(_record("ncbi", "333", gene_name="uncharacterized", locus_tag="PAS_chr1-1_0003"),),
    )
    conflict = classify_external_name_consistency(
        internal_gene_id="PAS_chr1-1_0004",
        internal_common_name="HRD1",
        internal_aliases=(),
        external_records=(_record("uniprot", "P44444", gene_name="PEP4", locus_tag="PAS_chr9_9999"),),
    )
    missing = classify_external_name_consistency(
        internal_gene_id="PAS_chr1-1_0005",
        internal_common_name="KAR2",
        internal_aliases=(),
        external_records=(),
    )

    assert exact.external_name_status == EXTERNAL_MATCH_CONFIRMED
    assert alias.external_name_status == EXTERNAL_ALIAS_CONFIRMED
    assert locus.external_name_status == EXTERNAL_LOCUS_CONFIRMED
    assert conflict.external_name_status == EXTERNAL_CONFLICT
    assert conflict.manual_review_reasons
    assert missing.external_name_status == EXTERNAL_REFERENCE_MISSING
    assert missing.manual_review_reasons == ("no external reference matched the current name-audit row",)


def test_external_name_resolution_selects_records_by_accession_name_or_locus() -> None:
    row = {
        "internal_gene_id": "PAS_chr1-1_0001",
        "internal_common_name": "SEC1",
        "internal_sequence_id": "YDR164C",
        "external_accession": "P12345",
        "external_gene_name": "SEC1",
        "external_locus_tag": "PAS_chr1-1_0001",
        "external_aliases": ("secretory protein 1",),
    }
    records = (
        _record("uniprot", "P12345", gene_name="SEC1", locus_tag="PAS_chr1-1_0001"),
        _record("sgd", "S000", gene_name="SCE1", locus_tag="YDR164C"),
        _record("uniprot", "P99999", gene_name="PEP4", locus_tag="PAS_chr9_9999"),
    )

    selected = select_external_records_for_name_audit_row(row, records)

    assert {record.primary_accession for record in selected} == {"P12345", "S000"}


def _record(
    source: str,
    accession: str,
    *,
    gene_name: str,
    locus_tag: str | None = None,
    aliases: tuple[str, ...] = (),
) -> ExternalReferenceRecord:
    return ExternalReferenceRecord(
        provenance=ExternalReferenceProvenance(
            source_database=source,
            source_version="test",
            source_url=f"https://example.test/{accession}",
            source_query=accession,
            retrieved_at="2026-07-09T00:00:00Z",
            raw_record_sha256="a" * 64,
        ),
        taxon_id="4922",
        organism="Komagataella phaffii",
        primary_accession=accession,
        gene_id=locus_tag,
        gene_name=gene_name,
        locus_tag=locus_tag,
        aliases=aliases,
    )
