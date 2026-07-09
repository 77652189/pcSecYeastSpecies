from __future__ import annotations

from pcsec_pichia.external_refs import (
    ExternalReferenceQuery,
    build_external_reference_queries,
    build_external_reference_queries_from_gene_catalog,
    build_external_reference_queries_from_homology_cache,
    build_external_reference_queries_from_ko_oe_candidate_rows,
    build_external_reference_queries_from_name_audit,
    dedupe_external_reference_queries,
    external_reference_query_fingerprint,
    normalize_external_query_name,
)
from pcsec_pichia.homology.cache_schema import HomologyCrosswalkRow, NameAuditRow
from pcsec_pichia.homology.crosswalk import write_homology_cache, write_name_audit_cache


def test_homology_cache_queries_cover_four_source_types_and_stable_fingerprint() -> None:
    rows = (
        _homology_row(
            internal_common_name="Sec1",
            query_symbol="SEC1",
            sce_orf="YDR164C",
            pichia_gene_id="PAS_chr1-1_0001",
            pichia_model_gene_id="PAS_chr1-1_0001",
            external_accession="P12345",
        ),
    )

    queries = build_external_reference_queries_from_homology_cache(rows)
    by_type = {query.query_type: query for query in queries}

    assert set(by_type) == {"sce_homolog", "pichia_gene", "model_gene", "external_accession"}
    assert by_type["sce_homolog"].query_value == "YDR164C"
    assert by_type["pichia_gene"].query_value == "PAS_chr1-1_0001"
    assert by_type["model_gene"].source_context == "homology_cache"
    assert by_type["external_accession"].preferred_sources == ("uniprot", "ncbi", "sgd")
    assert by_type["sce_homolog"].metadata["review_status"] == "model_ready_rbh_high_confidence"
    assert by_type["sce_homolog"].query_fingerprint == external_reference_query_fingerprint(by_type["sce_homolog"])
    repeated_by_type = {query.query_type: query for query in build_external_reference_queries_from_homology_cache(rows)}
    assert by_type["sce_homolog"].query_fingerprint == repeated_by_type["sce_homolog"].query_fingerprint


def test_homology_cache_path_loading_and_sce_fallback_warning(tmp_path) -> None:
    row = _homology_row(
        internal_common_name="Unresolved",
        query_symbol="ERO1",
        sce_orf=None,
        pichia_gene_id=None,
        pichia_model_gene_id=None,
        external_accession="",
        warnings=("manual review required",),
    )
    jsonl_path = tmp_path / "homology.jsonl"
    write_homology_cache((row,), jsonl_path, tmp_path / "homology.tsv")

    queries = build_external_reference_queries_from_homology_cache(jsonl_path)

    assert len(queries) == 1
    assert queries[0].query_type == "sce_homolog"
    assert queries[0].query_value == "ERO1"
    assert queries[0].warnings == (
        "manual review required",
        "homology cache row has no SCE ORF; falling back to query symbol.",
    )


def test_name_audit_queries_preserve_name_status_metadata_and_path_loading(tmp_path) -> None:
    row = _name_audit_row()
    jsonl_path = tmp_path / "name_audit.jsonl"
    write_name_audit_cache((row,), jsonl_path, tmp_path / "name_audit.tsv")

    queries = build_external_reference_queries_from_name_audit(jsonl_path)
    values_by_type = {query.query_type: query.query_value for query in queries}

    assert values_by_type == {
        "external_accession": "P12345",
        "model_gene": "PAS_chr1-1_0001",
        "pichia_gene": "PAS_chr1-1_0001",
        "sce_homolog": "YDR164C",
    }
    assert all(query.source_context == "name_audit" for query in queries)
    assert queries[0].metadata["name_consistency_status"] == "name_confirmed_by_rbh"


def test_gene_catalog_queries_use_model_gene_symbol_and_external_ids() -> None:
    rows = (
        {
            "canonical_gene_id": "PAS_chr2-1_0002",
            "standard_gene_symbol": "KAR2",
            "common_name": "Kar2",
            "external_ids": {"uniprot": ["Q11111", "Q22222"], "sgd": "YJL034W"},
            "recommended_use": "manual_review",
        },
    )

    queries = build_external_reference_queries_from_gene_catalog(rows)
    values = {(query.query_type, query.query_value) for query in queries}

    assert ("model_gene", "PAS_chr2-1_0002") in values
    assert ("pichia_gene", "KAR2") in values
    assert ("external_accession", "Q11111") in values
    assert ("external_accession", "Q22222") in values
    assert ("external_accession", "YJL034W") in values
    assert all(query.source_context == "gene_catalog" for query in queries)


def test_ko_oe_candidate_rows_keep_context_separate_from_gene_catalog() -> None:
    candidate_rows = (
        {
            "gene_id": "PAS_chr3-1_0003",
            "common_name": "PDI1",
            "intervention_type": "OE",
            "target_id": "hLF",
            "recommendation_tier": "model_executable",
            "external_ids": {"uniprot": "Q33333;Q44444"},
        },
    )

    queries = build_external_reference_queries_from_ko_oe_candidate_rows(candidate_rows)

    assert {query.source_context for query in queries} == {"ko_oe_candidate_rows"}
    assert {query.query_type for query in queries} == {"external_accession", "model_gene", "pichia_gene"}
    assert any(query.metadata.get("intervention_type") == "OE" for query in queries)


def test_dedupe_is_stable_and_merges_source_rows_and_warnings() -> None:
    left = ExternalReferenceQuery(
        query_type="pichia_gene",
        query_value="sec1",
        source_context="gene_catalog",
        source_id="gene_catalog",
        source_row_id="row-a",
        warnings=("left warning",),
    )
    right = ExternalReferenceQuery(
        query_type="pichia_gene",
        query_value="SEC1",
        source_context="gene_catalog",
        source_id="gene_catalog",
        source_row_id="row-b",
        warnings=("right warning",),
    )

    deduped = dedupe_external_reference_queries((right, left))

    assert len(deduped) == 1
    assert deduped[0].source_row_id == "row-b;row-a"
    assert deduped[0].metadata["merged_source_row_ids"] == "row-b;row-a"
    assert deduped[0].warnings == ("right warning", "left warning")
    assert normalize_external_query_name(" sec1 ") == "SEC1"


def test_combined_builder_does_not_perform_network_or_change_query_tiers() -> None:
    queries = build_external_reference_queries(
        homology_cache=(
            _homology_row(
                internal_common_name="Sec1",
                query_symbol="SEC1",
                sce_orf="YDR164C",
                pichia_gene_id="PAS_chr1-1_0001",
                pichia_model_gene_id="PAS_chr1-1_0001",
                external_accession="P12345",
            ),
        ),
        gene_catalog_rows=(
            {
                "canonical_gene_id": "PAS_chr1-1_0001",
                "standard_gene_symbol": "SEC1",
                "external_ids": {"uniprot": "P12345"},
            },
        ),
        ko_oe_candidate_rows=(
            {
                "gene_id": "PAS_chr1-1_0001",
                "common_name": "SEC1",
                "intervention_type": "KO",
                "recommendation_tier": "not_recommended_growth_risk",
            },
        ),
    )

    assert queries
    assert all(not hasattr(query, "recommendation_tier") for query in queries)
    assert {query.source_context for query in queries} == {
        "gene_catalog",
        "homology_cache",
        "ko_oe_candidate_rows",
    }


def _homology_row(
    *,
    internal_common_name: str,
    query_symbol: str,
    sce_orf: str | None,
    pichia_gene_id: str | None,
    pichia_model_gene_id: str | None,
    external_accession: str,
    warnings: tuple[str, ...] = (),
) -> HomologyCrosswalkRow:
    return HomologyCrosswalkRow(
        internal_common_name=internal_common_name,
        query_symbol=query_symbol,
        sce_orf=sce_orf,
        pichia_gene_id=pichia_gene_id,
        pichia_model_gene_id=pichia_model_gene_id,
        is_rbh=bool(sce_orf and pichia_gene_id),
        identity_pct=80.0 if pichia_gene_id else None,
        evalue=1e-40 if pichia_gene_id else None,
        query_coverage=95.0 if pichia_gene_id else None,
        subject_coverage=94.0 if pichia_gene_id else None,
        in_model_gene_index=bool(pichia_model_gene_id),
        review_status="model_ready_rbh_high_confidence" if pichia_model_gene_id else "unresolved_query_symbol",
        warnings=warnings,
        external_accession=external_accession,
        external_gene_name=query_symbol,
        external_locus_tag=pichia_gene_id or "",
        external_aliases=(query_symbol,),
    )


def _name_audit_row() -> NameAuditRow:
    return NameAuditRow(
        internal_gene_id="PAS_chr1-1_0001",
        internal_common_name="Sec1",
        internal_sequence_id="YDR164C",
        external_accession="P12345",
        external_gene_name="SEC1",
        external_locus_tag="PAS_chr1-1_0001",
        external_aliases=("SEC1",),
        identity_pct=80.0,
        query_coverage=95.0,
        subject_coverage=94.0,
        evalue=1e-40,
        is_rbh=True,
        in_model_gene_index=True,
        name_consistency_status="name_confirmed_by_rbh",
        review_status="model_ready_rbh_high_confidence",
    )
