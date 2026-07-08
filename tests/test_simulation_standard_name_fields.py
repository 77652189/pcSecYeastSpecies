from __future__ import annotations

from pcsec_pichia.reports import CANDIDATE_COLUMNS, normalize_candidate_explanation_row
from pcsec_pichia.services.gene_id_standardization import PichiaGeneIdStandardName, build_standard_name_lookup

from app.services.pichia_screen_preview_service import _gene_row, _reaction_row
from app.ui.views.simulation_display import CANDIDATE_DISPLAY_COLUMNS, candidate_row_label
from app.ui.views.simulation_gene_inputs import gene_mapping_rows_for_display


def test_preview_gene_rows_pass_through_standard_name_fields() -> None:
    lookup = build_standard_name_lookup(
        (
            PichiaGeneIdStandardName(
                gene_id="PAS_chr2-1_0140",
                display_name="KAR2",
                standard_symbol="KAR2",
                protein_name="BiP molecular chaperone",
                annotation_sources=("UniProt",),
                annotation_confidence="high_exact_locus_tag",
            ),
        )
    )

    row = _gene_row(
        "KAR2",
        "KO",
        resolved=True,
        canonical_gene_id="PAS_chr2-1_0140",
        standard_name_lookup=lookup,
    )
    reaction = _reaction_row("sec_Kar2p_complex_formation", "KO_reaction", resolved=True)

    assert row["input_id"] == "KAR2"
    assert row["canonical_gene_id"] == "PAS_chr2-1_0140"
    assert row["standard_symbol"] == "KAR2"
    assert row["gene_display_name"] == "KAR2"
    assert row["standard_name_status"] == "annotated"
    assert reaction["standard_name_status"] == "not_gene_candidate"
    assert reaction["standard_symbol"] == ""


def test_candidate_report_and_ui_display_standard_name_fields() -> None:
    for column in (
        "gene_display_name",
        "standard_symbol",
        "protein_name",
        "annotation_confidence",
        "standard_name_status",
    ):
        assert column in CANDIDATE_COLUMNS
        assert column in CANDIDATE_DISPLAY_COLUMNS

    normalized = normalize_candidate_explanation_row(
        {
            "candidate_id": "PAS_chr2-1_0140",
            "gene_id": "PAS_chr2-1_0140",
            "canonical_gene_id": "PAS_chr2-1_0140",
            "standard_symbol": "KAR2",
            "gene_display_name": "KAR2",
            "intervention_type": "KO",
            "effect_label": "提升分泌",
            "secretory_process": "ER folding",
            "mapping_confidence": "high",
            "gpr_role": "single_gene",
            "simulation_basis": "gpr_gene_deletion",
            "delta_objective": 0.2,
            "status": "0",
            "success": True,
        }
    )

    assert normalized["candidate_id"] == "PAS_chr2-1_0140"
    assert normalized["gene_id"] == "PAS_chr2-1_0140"
    assert normalized["standard_symbol"] == "KAR2"
    assert "`KAR2 (PAS_chr2-1_0140)`" in normalized["summary"]


def test_gene_mapping_display_includes_standard_name_fields() -> None:
    frame = gene_mapping_rows_for_display(
        [
            {
                "input_gene_id": "KAR2",
                "canonical_gene_id": "PAS_chr2-1_0140",
                "gene_display_name": "KAR2",
                "standard_symbol": "KAR2",
                "protein_name": "BiP molecular chaperone",
                "annotation_confidence": "high_exact_locus_tag",
                "standard_name_status": "annotated",
                "reaction_id": "sec_Kar2p_complex_formation",
            }
        ]
    )

    assert frame.loc[0, "标准符号"] == "KAR2"
    assert frame.loc[0, "蛋白名称"] == "BiP molecular chaperone"


def test_candidate_row_label_prefers_standard_symbol() -> None:
    import pandas as pd

    label = candidate_row_label(
        0,
        pd.Series(
            {
                "standard_symbol": "KAR2",
                "gene_display_name": "BiP",
                "gene_id": "PAS_chr2-1_0140",
                "effect_label": "提升分泌",
                "delta_objective": 0.2,
            }
        ),
    )

    assert label.startswith("1. KAR2 |")
