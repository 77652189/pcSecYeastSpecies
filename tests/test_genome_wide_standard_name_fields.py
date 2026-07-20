from __future__ import annotations

from pcsec_pichia.services.gene_id_standardization import PichiaGeneIdStandardName, build_standard_name_lookup

from python_pichia.tools.run_genome_wide_ko_oe_screen_parallel import _row_to_csv_record
from app.ui.views.simulation import _prefill_field_values


def test_genome_wide_csv_record_includes_standard_name_fields() -> None:
    lookup = build_standard_name_lookup(
        (
            PichiaGeneIdStandardName(
                gene_id="PAS_chr2-1_0140",
                display_name="KAR2",
                standard_symbol="KAR2",
                protein_name="BiP molecular chaperone",
                external_ids={"uniprot": "C4R8K4"},
                annotation_sources=("UniProt",),
                annotation_confidence="high_exact_locus_tag",
            ),
        )
    )
    record = _row_to_csv_record(
        {
            "target_id": "hLF",
            "gene_id": "PAS_chr2-1_0140",
            "candidate_kind": "gene",
            "intervention_type": "KO",
            "support_status": "ko_runnable_gpr_gene_deletion",
            "secretory_process": "ER folding",
            "gpr_role": "single_gene",
            "mapping_confidence": "high",
            "max_feasible_mu": 0.1,
            "secretion_at_max_feasible_mu": 1.2,
            "wildtype_max_feasible_mu": 0.1,
            "wildtype_secretion_at_max_feasible_mu": 1.0,
            "growth_retention_ratio": 1.0,
            "secretion_ratio_vs_wildtype": 1.2,
            "skipped_reason": None,
        },
        lookup,
    )

    assert record["gene_display_name"] == "KAR2"
    assert record["standard_symbol"] == "KAR2"
    assert record["protein_name"] == "BiP molecular chaperone"
    assert record["external_ids"] == '{"uniprot": "C4R8K4"}'
    assert record["annotation_sources"] == '["UniProt"]'
    assert record["standard_name_status"] == "annotated"
    assert record["gene_id"] == "PAS_chr2-1_0140"
    assert record["secretion_ratio_vs_wildtype"] == 1.2


def test_genome_wide_non_gene_candidate_is_not_standard_named() -> None:
    record = _row_to_csv_record(
        {
            "target_id": "hLF",
            "gene_id": "sec_Kar2p_complex_formation",
            "candidate_kind": "catalog_reaction",
            "intervention_type": "OE",
            "support_status": "reaction_level_diagnostic",
            "secretory_process": "ER folding",
            "gpr_role": "reaction_level",
            "mapping_confidence": "medium",
            "max_feasible_mu": None,
            "secretion_at_max_feasible_mu": None,
            "wildtype_max_feasible_mu": 0.1,
            "wildtype_secretion_at_max_feasible_mu": 1.0,
            "growth_retention_ratio": None,
            "secretion_ratio_vs_wildtype": None,
            "skipped_reason": None,
        },
        {},
    )

    assert record["standard_name_status"] == "not_gene_candidate"
    assert record["standard_symbol"] == ""


def test_genome_wide_verify_prefill_keeps_gene_id_not_standard_symbol() -> None:
    values = _prefill_field_values("PAS_chr2-1_0140", "KO", "gene")

    assert values["pichia_draft_ko_genes"] == "PAS_chr2-1_0140"
    assert "KAR2" not in values.values()
