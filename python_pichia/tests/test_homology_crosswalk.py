from __future__ import annotations

from pathlib import Path

from pcsec_pichia.homology.cache_schema import (
    BlastConfig,
    BlastHit,
    CatalogHomologyQuery,
    ExternalNameReference,
    NameAuditRow,
    ProteinRecord,
    ReciprocalBestHit,
)
from pcsec_pichia.homology.crosswalk import (
    build_external_database_crosschecks,
    build_homology_crosswalk,
    build_name_audit_rows,
    build_rule_transfer_audit_rows,
    load_external_name_reference_cache,
    load_homology_cache,
    load_name_audit_cache,
    load_rule_transfer_audit_cache,
    merge_external_crosschecks_into_name_audit,
    summarize_homology_audits,
    write_homology_cache,
    write_name_audit_cache,
    write_rule_transfer_audit_cache,
)
from pcsec_pichia.homology.review_rules import (
    EXTERNAL_ALIAS_CONFIRMED,
    EXTERNAL_CONFLICT,
    EXTERNAL_LOCUS_CONFIRMED,
    EXTERNAL_MATCH_CONFIRMED,
    LOW_IDENTITY_REVIEW_REQUIRED,
    MODEL_READY_RBH_HIGH_CONFIDENCE,
    NO_RECIPROCAL_HIT,
    PICHIA_LOCUS_CONFIRMED_BY_RBH,
    RBH_NOT_IN_MODEL,
    RULE_TRANSFER_LOW_CONFIDENCE,
    RULE_TRANSFER_NOT_SUPPORTED,
    RULE_TRANSFER_READY,
    RULE_TRANSFER_SUPPORTED_NOT_MODEL_OPERABLE,
    RULE_TRANSFER_UNRESOLVED,
)
from pcsec_pichia.services.homology_evidence import build_homology_evidence_map


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


def _low_identity_hit(query: str, subject: str) -> BlastHit:
    hit = _hit(query, subject)
    return BlastHit(
        query_id=hit.query_id,
        subject_id=hit.subject_id,
        identity_pct=20,
        alignment_length=hit.alignment_length,
        query_length=hit.query_length,
        subject_length=hit.subject_length,
        evalue=hit.evalue,
        bitscore=hit.bitscore,
        query_coverage=hit.query_coverage,
        subject_coverage=hit.subject_coverage,
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


def test_build_name_audit_rows_treats_pichia_locus_as_stable_name() -> None:
    rows = build_homology_crosswalk(
        (CatalogHomologyQuery(internal_common_name="CDC48", query_symbol="CDC48"),),
        (ProteinRecord(organism="sce", gene_id="YDL126C", symbol="CDC48", sequence="M"),),
        (
            ProteinRecord(
                organism="pichia",
                gene_id="PAS_FragD_0026",
                symbol="PAS_FragD_0026",
                aliases=("C4R9A6",),
                sequence="M",
            ),
        ),
        {"PAS_FragD_0026"},
        (
            ReciprocalBestHit(
                "YDL126C",
                "PAS_FragD_0026",
                True,
                _hit("YDL126C", "PAS_FragD_0026"),
                _hit("PAS_FragD_0026", "YDL126C"),
            ),
        ),
        BlastConfig(),
    )

    audit = build_name_audit_rows(rows)

    assert audit[0].name_consistency_status == PICHIA_LOCUS_CONFIRMED_BY_RBH


def test_external_reference_cache_loads_jsonl_and_tsv(tmp_path: Path) -> None:
    jsonl = tmp_path / "external_refs.jsonl"
    tsv = tmp_path / "external_refs.tsv"
    jsonl.write_text(
        (
            '{"source_database":"UniProt","source_version":"2026_01","taxon":"Komagataella phaffii",'
            '"accession":"C4R","gene_name":"KAR2","locus_tag":"PAS_chr2-1_0140",'
            '"aliases":["BiP"],"retrieved_at":"2026-07-07","warnings":["offline snapshot"]}\n'
        ),
        encoding="utf-8",
    )
    tsv.write_text(
        "\t".join(
            [
                "source_database",
                "source_version",
                "taxon",
                "accession",
                "gene_name",
                "locus_tag",
                "aliases",
                "retrieved_at",
                "warnings",
            ]
        )
        + "\n"
        + "\t".join(["SGD", "R64", "Saccharomyces cerevisiae", "S000001", "KAR2", "YJL034W", "BiP;GRP78", "2026-07-07", ""])
        + "\n",
        encoding="utf-8",
    )

    jsonl_refs = load_external_name_reference_cache(jsonl)
    tsv_refs = load_external_name_reference_cache(tsv)

    assert jsonl_refs[0].source_database == "UniProt"
    assert jsonl_refs[0].aliases == ("BiP",)
    assert jsonl_refs[0].warnings == ("offline snapshot",)
    assert tsv_refs[0].source_database == "SGD"
    assert tsv_refs[0].aliases == ("BiP", "GRP78")


def test_external_crosscheck_merges_status_without_overriding_rbh_facts() -> None:
    name_rows = (
        _name_row("KAR2 / BiP", "YJL034W", "C4R", "KAR2", "PAS_chr2-1_0140"),
        _name_row("DOA10", "YIL030C", "C4D", "DOA10", "PAS_chr3_0123"),
        _name_row("VPS1", "YOR069W", "", "", "PAS_chr1_0440"),
        _name_row("HRD1", "YOL013C", "C4H", "HRD1", "PAS_chr4_0156"),
        _name_row("CDC48", "YDL126C", "C4R9A6", "PAS_FragD_0026", "PAS_FragD_0026"),
    )
    references = (
        ExternalNameReference("UniProt", "2026_01", "Komagataella phaffii", "C4R", "KAR2", "PAS_chr2-1_0140"),
        ExternalNameReference("UniProt", "2026_01", "Komagataella phaffii", "C4D", "SSM4", "PAS_chr3_0123", ("DOA10",)),
        ExternalNameReference("NCBI", "2026_01", "Komagataella phaffii", "", "", "PAS_chr1_0440"),
        ExternalNameReference("UniProt", "2026_01", "Komagataella phaffii", "C4H", "PEP4", "PAS_chr4_9999"),
        ExternalNameReference("NCBI", "gene", "Komagataella phaffii", "8200528", "PAS_FragD_0026", ""),
        ExternalNameReference("UniProt", "2026_02", "Komagataella phaffii", "C4R9A6", "C4R9A6_KOMPG", "PAS_FragD_0026"),
    )

    crosschecks = build_external_database_crosschecks(name_rows, references)
    merged = merge_external_crosschecks_into_name_audit(name_rows, crosschecks)
    by_name = {row.internal_common_name: row for row in merged}

    assert by_name["KAR2 / BiP"].external_crosscheck_status == EXTERNAL_MATCH_CONFIRMED
    assert by_name["DOA10"].external_crosscheck_status == EXTERNAL_ALIAS_CONFIRMED
    assert by_name["VPS1"].external_crosscheck_status == EXTERNAL_LOCUS_CONFIRMED
    assert by_name["HRD1"].external_crosscheck_status == EXTERNAL_CONFLICT
    assert by_name["CDC48"].external_crosscheck_status == EXTERNAL_MATCH_CONFIRMED
    assert by_name["HRD1"].review_status == MODEL_READY_RBH_HIGH_CONFIDENCE
    assert by_name["HRD1"].is_rbh is True
    assert by_name["HRD1"].in_model_gene_index is True
    assert by_name["HRD1"].external_crosscheck_warnings


def test_rule_transfer_audit_rows_cover_ready_not_model_low_confidence_and_no_rbh() -> None:
    queries = (
        CatalogHomologyQuery(internal_common_name="KAR2", query_symbol="KAR2"),
        CatalogHomologyQuery(internal_common_name="PDI1", query_symbol="PDI1"),
        CatalogHomologyQuery(internal_common_name="HRD1", query_symbol="HRD1"),
        CatalogHomologyQuery(internal_common_name="SEC12", query_symbol="SEC12"),
        CatalogHomologyQuery(internal_common_name="MISSING", query_symbol="MISSING"),
    )
    sce_records = (
        ProteinRecord(organism="sce", gene_id="YJL034W", symbol="KAR2", sequence="M"),
        ProteinRecord(organism="sce", gene_id="YCL043C", symbol="PDI1", sequence="M"),
        ProteinRecord(organism="sce", gene_id="YOL013C", symbol="HRD1", sequence="M"),
        ProteinRecord(organism="sce", gene_id="YNR026C", symbol="SEC12", sequence="M"),
    )
    pichia_records = (
        ProteinRecord(organism="pichia", gene_id="PAS_chr2-1_0140", symbol="KAR2", sequence="M"),
        ProteinRecord(organism="pichia", gene_id="PAS_chr4_0844", symbol="PDI1", sequence="M"),
        ProteinRecord(organism="pichia", gene_id="PAS_chr4_0156", symbol="HRD1", sequence="M"),
        ProteinRecord(organism="pichia", gene_id="PAS_chr4_0606", symbol="SEC12", sequence="M"),
    )
    rbh_calls = (
        ReciprocalBestHit("YJL034W", "PAS_chr2-1_0140", True, _hit("YJL034W", "PAS_chr2-1_0140"), _hit("PAS_chr2-1_0140", "YJL034W")),
        ReciprocalBestHit("YCL043C", "PAS_chr4_0844", True, _hit("YCL043C", "PAS_chr4_0844"), _hit("PAS_chr4_0844", "YCL043C")),
        ReciprocalBestHit("YOL013C", "PAS_chr4_0156", True, _low_identity_hit("YOL013C", "PAS_chr4_0156"), _low_identity_hit("PAS_chr4_0156", "YOL013C")),
        ReciprocalBestHit("YNR026C", "PAS_chr4_0606", False, _hit("YNR026C", "PAS_chr4_0606"), _hit("PAS_chr4_0606", "YCR067C"), "reverse_best_is_YCR067C"),
    )

    crosswalk = build_homology_crosswalk(
        queries,
        sce_records,
        pichia_records,
        {"PAS_chr2-1_0140"},
        rbh_calls,
        BlastConfig(),
    )
    rule_rows = build_rule_transfer_audit_rows(crosswalk)

    by_symbol = {row.query_symbol: row for row in rule_rows}
    assert by_symbol["KAR2"].rule_transfer_status == RULE_TRANSFER_READY
    assert by_symbol["PDI1"].rule_transfer_status == RULE_TRANSFER_SUPPORTED_NOT_MODEL_OPERABLE
    assert by_symbol["HRD1"].homology_review_status == LOW_IDENTITY_REVIEW_REQUIRED
    assert by_symbol["HRD1"].rule_transfer_status == RULE_TRANSFER_LOW_CONFIDENCE
    assert by_symbol["SEC12"].homology_review_status == NO_RECIPROCAL_HIT
    assert by_symbol["SEC12"].rule_transfer_status == RULE_TRANSFER_NOT_SUPPORTED
    assert by_symbol["MISSING"].rule_transfer_status == RULE_TRANSFER_UNRESOLVED


def test_name_and_rule_transfer_cache_outputs_are_stable(tmp_path: Path) -> None:
    crosswalk = build_homology_crosswalk(
        (CatalogHomologyQuery(internal_common_name="KAR2", query_symbol="KAR2"),),
        (ProteinRecord(organism="sce", gene_id="YJL034W", symbol="KAR2", sequence="M"),),
        (ProteinRecord(organism="pichia", gene_id="PAS_chr2-1_0140", symbol="KAR2", accession="C4Q", sequence="M"),),
        {"PAS_chr2-1_0140"},
        (ReciprocalBestHit("YJL034W", "PAS_chr2-1_0140", True, _hit("YJL034W", "PAS_chr2-1_0140"), _hit("PAS_chr2-1_0140", "YJL034W")),),
        BlastConfig(),
    )
    name_rows = build_name_audit_rows(crosswalk)
    rule_rows = build_rule_transfer_audit_rows(crosswalk)

    write_name_audit_cache(name_rows, tmp_path / "name.jsonl", tmp_path / "name.tsv")
    write_rule_transfer_audit_cache(rule_rows, tmp_path / "rule.jsonl", tmp_path / "rule.tsv")

    name_header = (tmp_path / "name.tsv").read_text(encoding="utf-8").splitlines()[0].split("\t")
    rule_header = (tmp_path / "rule.tsv").read_text(encoding="utf-8").splitlines()[0].split("\t")
    assert "name_consistency_status" in name_header
    assert "rule_transfer_status" in rule_header
    assert load_name_audit_cache(tmp_path / "name.jsonl")[0].name_consistency_status == "name_confirmed_by_rbh"
    assert load_rule_transfer_audit_cache(tmp_path / "rule.jsonl")[0].rule_transfer_status == RULE_TRANSFER_READY

    evidence = build_homology_evidence_map(name_rows=name_rows, rule_rows=rule_rows)
    kar2 = evidence["pas_chr2-1_0140"]
    assert kar2.name_consistency_status == "name_confirmed_by_rbh"
    assert kar2.rule_transfer_status == RULE_TRANSFER_READY
    assert kar2.homology_review_status == MODEL_READY_RBH_HIGH_CONFIDENCE


def test_homology_audit_summary_counts_all_three_outputs() -> None:
    crosswalk = build_homology_crosswalk(
        (CatalogHomologyQuery(internal_common_name="KAR2", query_symbol="KAR2"),),
        (ProteinRecord(organism="sce", gene_id="YJL034W", symbol="KAR2", sequence="M"),),
        (ProteinRecord(organism="pichia", gene_id="PAS_chr2-1_0140", symbol="KAR2", sequence="M"),),
        {"PAS_chr2-1_0140"},
        (ReciprocalBestHit("YJL034W", "PAS_chr2-1_0140", True, _hit("YJL034W", "PAS_chr2-1_0140"), _hit("PAS_chr2-1_0140", "YJL034W")),),
        BlastConfig(),
    )
    name_rows = build_name_audit_rows(crosswalk)
    rule_rows = build_rule_transfer_audit_rows(crosswalk)

    summary = summarize_homology_audits(
        blast_status="completed",
        homology_rows=crosswalk,
        name_audit_rows=name_rows,
        rule_transfer_rows=rule_rows,
    )

    assert summary.homology_row_count == 1
    assert summary.name_audit_row_count == 1
    assert summary.rule_transfer_status_counts == {RULE_TRANSFER_READY: 1}


def _name_row(
    common_name: str,
    sce_orf: str,
    accession: str,
    gene_name: str,
    locus_tag: str,
) -> NameAuditRow:
    return NameAuditRow(
        internal_gene_id=locus_tag,
        internal_common_name=common_name,
        internal_sequence_id=sce_orf,
        external_accession=accession,
        external_gene_name=gene_name,
        external_locus_tag=locus_tag,
        external_aliases=(),
        identity_pct=75.0,
        query_coverage=95.0,
        subject_coverage=95.0,
        evalue=1e-100,
        is_rbh=True,
        in_model_gene_index=True,
        name_consistency_status="name_confirmed_by_rbh",
        review_status=MODEL_READY_RBH_HIGH_CONFIDENCE,
        warnings=(),
    )
