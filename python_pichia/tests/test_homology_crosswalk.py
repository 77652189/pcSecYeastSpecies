from __future__ import annotations

from pathlib import Path

from pcsec_pichia.homology.cache_schema import (
    BlastConfig,
    BlastHit,
    CatalogHomologyQuery,
    ProteinRecord,
    ReciprocalBestHit,
)
from pcsec_pichia.homology.crosswalk import (
    build_homology_crosswalk,
    build_name_audit_rows,
    load_homology_cache,
    write_homology_cache,
)
from pcsec_pichia.homology.review_rules import MODEL_READY_RBH_HIGH_CONFIDENCE, RBH_NOT_IN_MODEL


def _hit(query: str, subject: str) -> BlastHit:
    return BlastHit(
        query_id=query,
        subject_id=subject,
        identity_pct=70,
        alignment_length=100,
        query_length=100,
        subject_length=100,
        evalue=1e-80,
        bitscore=250,
        query_coverage=100,
        subject_coverage=100,
    )


def test_build_homology_crosswalk_keeps_model_operability_separate() -> None:
    queries = (
        CatalogHomologyQuery(internal_common_name="KAR2 / BiP", query_symbol="KAR2"),
        CatalogHomologyQuery(internal_common_name="PDI1", query_symbol="PDI1"),
    )
    sce_records = (
        ProteinRecord(organism="sce", gene_id="YJL034W", symbol="KAR2", sequence="M"),
        ProteinRecord(organism="sce", gene_id="YCL043C", symbol="PDI1", sequence="M"),
    )
    pichia_records = (
        ProteinRecord(organism="pichia", gene_id="PAS_chr2-1_0140", symbol="KAR2", sequence="M"),
        ProteinRecord(organism="pichia", gene_id="PAS_chr4_0844", symbol="PDI1", sequence="M"),
    )
    rbh_calls = (
        ReciprocalBestHit("YJL034W", "PAS_chr2-1_0140", True, _hit("YJL034W", "PAS_chr2-1_0140"), _hit("PAS_chr2-1_0140", "YJL034W")),
        ReciprocalBestHit("YCL043C", "PAS_chr4_0844", True, _hit("YCL043C", "PAS_chr4_0844"), _hit("PAS_chr4_0844", "YCL043C")),
    )

    rows = build_homology_crosswalk(
        queries,
        sce_records,
        pichia_records,
        {"PAS_chr2-1_0140"},
        rbh_calls,
        BlastConfig(),
    )

    by_symbol = {row.query_symbol: row for row in rows}
    assert by_symbol["KAR2"].review_status == MODEL_READY_RBH_HIGH_CONFIDENCE
    assert by_symbol["KAR2"].pichia_model_gene_id == "PAS_chr2-1_0140"
    assert by_symbol["PDI1"].review_status == RBH_NOT_IN_MODEL
    assert by_symbol["PDI1"].pichia_model_gene_id is None


def test_cache_output_jsonl_tsv_fields_are_stable(tmp_path: Path) -> None:
    row = build_homology_crosswalk(
        (CatalogHomologyQuery(internal_common_name="KAR2", query_symbol="KAR2"),),
        (ProteinRecord(organism="sce", gene_id="YJL034W", symbol="KAR2", sequence="M"),),
        (ProteinRecord(organism="pichia", gene_id="PAS_chr2-1_0140", symbol="KAR2", sequence="M"),),
        {"PAS_chr2-1_0140"},
        (ReciprocalBestHit("YJL034W", "PAS_chr2-1_0140", True, _hit("YJL034W", "PAS_chr2-1_0140"), _hit("PAS_chr2-1_0140", "YJL034W")),),
        BlastConfig(),
    )
    jsonl = tmp_path / "cache.jsonl"
    tsv = tmp_path / "cache.tsv"

    result = write_homology_cache(row, jsonl, tsv)

    assert result.row_count == 1
    assert tsv.read_text(encoding="utf-8").splitlines()[0].split("\t") == [
        "internal_common_name",
        "query_symbol",
        "sce_orf",
        "pichia_gene_id",
        "pichia_model_gene_id",
        "is_rbh",
        "identity_pct",
        "evalue",
        "query_coverage",
        "subject_coverage",
        "in_model_gene_index",
        "review_status",
        "warnings",
        "external_accession",
        "external_gene_name",
        "external_locus_tag",
        "external_aliases",
    ]
    loaded = load_homology_cache(jsonl)
    assert loaded[0].query_symbol == "KAR2"


def test_build_name_audit_rows_flags_name_conflict() -> None:
    rows = build_homology_crosswalk(
        (CatalogHomologyQuery(internal_common_name="SSA1", query_symbol="SSA1"),),
        (ProteinRecord(organism="sce", gene_id="YAL005C", symbol="SSA1", sequence="M"),),
        (ProteinRecord(organism="pichia", gene_id="PAS_chr4_0552", symbol="SSA2", sequence="M"),),
        {"PAS_chr4_0552"},
        (ReciprocalBestHit("YAL005C", "PAS_chr4_0552", True, _hit("YAL005C", "PAS_chr4_0552"), _hit("PAS_chr4_0552", "YAL005C")),),
        BlastConfig(),
    )

    audit = build_name_audit_rows(rows)

    assert audit[0].name_consistency_status == "sequence_name_conflict"
