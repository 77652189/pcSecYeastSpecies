from __future__ import annotations

from pcsec_pichia.services.gene_id_standardization import (
    build_standard_name_lookup,
    build_pichia_gene_id_standardization_rows,
    enrich_gene_standard_name_fields,
    load_pichia_gene_id_standardization_cache,
    summarize_pichia_gene_id_standardization_rows,
    standard_name_fields_for_csv,
    write_pichia_gene_id_standardization_cache,
)


def test_gene_id_standardization_uses_gene_id_primary_key_and_annotation_tiers(tmp_path) -> None:
    source_rows = [
        {
            "gene_id": "PAS_chr2-1_0140",
            "display_name": "KAR2",
            "standard_gene_symbol": "KAR2",
            "protein_name": "BiP molecular chaperone",
            "external_ids": {"uniprot": "C4R8K4", "ncbi_gene": "8199000"},
            "evidence_sources": ["UniProt", "NCBI Gene", "KEGG", "RefSeq"],
            "evidence_confidence": "high_exact_locus_tag",
            "ko_support_status": "ko_runnable_gpr_gene_deletion",
            "oe_support_status": "oe_runnable_reaction_proxy",
            "affected_reactions": ["sec_Kar2p_complex_formation"],
            "gpr_role": "single_gene",
        },
        {
            "gene_id": "AT250_GQ_6803479",
            "display_name": "AT250_GQ_6803479",
            "standard_gene_symbol": "",
            "protein_name": "",
            "external_ids": {},
            "evidence_sources": [],
            "evidence_confidence": "low_model_only",
            "ko_support_status": "ko_no_gpr_effect",
            "oe_support_status": "oe_no_gpr_effect",
            "affected_reactions": [],
            "gpr_role": "unresolved",
        },
        {
            "gene_id": "PAS_chr1-4_0601",
            "display_name": "KEGG-only gene",
            "external_ids": {},
            "evidence_sources": ["KEGG"],
            "evidence_confidence": "medium_exact_kegg_locus_tag",
            "ko_support_status": "ko_runnable_gpr_gene_deletion",
            "oe_support_status": "oe_no_gpr_effect",
        },
    ]

    rows = build_pichia_gene_id_standardization_rows(source_rows)
    output_path = tmp_path / "pichia_gene_id_standardization.json"
    write_pichia_gene_id_standardization_cache(rows, output_path)
    loaded = load_pichia_gene_id_standardization_cache(output_path)
    summary = summarize_pichia_gene_id_standardization_rows(loaded)

    assert [row.gene_id for row in loaded] == [
        "PAS_chr2-1_0140",
        "AT250_GQ_6803479",
        "PAS_chr1-4_0601",
    ]
    assert loaded[0].standard_symbol == "KAR2"
    assert loaded[0].annotation_sources == ("UniProt", "NCBI Gene", "KEGG", "RefSeq")
    assert loaded[0].annotation_confidence == "high_exact_locus_tag"
    assert loaded[0].model_operable is True
    assert loaded[0].gpr_status == "ko_and_oe_model_executable"
    assert loaded[1].display_name == "AT250_GQ_6803479"
    assert loaded[1].annotation_sources == ("model_only",)
    assert loaded[1].annotation_confidence == "low_model_only"
    assert loaded[1].gpr_status == "model_gene_no_gpr_effect"
    assert summary["total_genes"] == 3
    assert summary["annotated_gene_count"] == 2
    assert summary["model_only_count"] == 1
    assert summary["model_only_gene_ids"] == ["AT250_GQ_6803479"]


def test_gene_id_standardization_rejects_missing_or_duplicate_gene_ids() -> None:
    missing_gene_rows = [{"gene_id": "", "display_name": "missing"}]
    duplicate_rows = [{"gene_id": "PAS1"}, {"gene_id": "PAS1"}]

    try:
        build_pichia_gene_id_standardization_rows(missing_gene_rows)
    except ValueError as exc:
        assert "gene_id" in str(exc)
    else:  # pragma: no cover - documents the expected failure path
        raise AssertionError("missing gene_id should fail")

    try:
        build_pichia_gene_id_standardization_rows(duplicate_rows)
    except ValueError as exc:
        assert "duplicate" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("duplicate gene_id should fail")


def test_standard_name_enrichment_statuses_and_csv_serialization() -> None:
    rows = build_pichia_gene_id_standardization_rows(
        [
            {
                "gene_id": "PAS_chr2-1_0140",
                "display_name": "KAR2",
                "standard_gene_symbol": "KAR2",
                "protein_name": "BiP molecular chaperone",
                "external_ids": {"uniprot": "C4R8K4"},
                "evidence_sources": ["UniProt"],
                "evidence_confidence": "high_exact_locus_tag",
            },
            {
                "gene_id": "AT250_GQ_6803479",
                "display_name": "AT250_GQ_6803479",
                "external_ids": {},
                "evidence_sources": [],
            },
        ]
    )
    lookup = build_standard_name_lookup(rows)

    annotated = enrich_gene_standard_name_fields({"gene_id": "PAS_chr2-1_0140", "candidate_kind": "gene"}, lookup)
    model_only = enrich_gene_standard_name_fields({"gene_id": "AT250_GQ_6803479", "candidate_kind": "gene"}, lookup)
    missing = enrich_gene_standard_name_fields({"gene_id": "PAS_missing", "candidate_kind": "gene"}, lookup)
    reaction = enrich_gene_standard_name_fields(
        {"gene_id": "sec_Kar2p_complex_formation", "candidate_kind": "catalog_reaction"},
        lookup,
    )
    csv_fields = standard_name_fields_for_csv(annotated)

    assert annotated["gene_display_name"] == "KAR2"
    assert annotated["standard_symbol"] == "KAR2"
    assert annotated["protein_name"] == "BiP molecular chaperone"
    assert annotated["annotation_confidence"] == "high_exact_locus_tag"
    assert annotated["standard_name_status"] == "annotated"
    assert model_only["standard_name_status"] == "model_only"
    assert missing["standard_name_status"] == "missing_standard_name"
    assert reaction["standard_name_status"] == "not_gene_candidate"
    assert csv_fields["external_ids"] == '{"uniprot": "C4R8K4"}'
    assert csv_fields["annotation_sources"] == '["UniProt"]'
    assert "{" in csv_fields["external_ids"] and "'" not in csv_fields["external_ids"]
