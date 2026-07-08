from __future__ import annotations

from pcsec_pichia.services.gene_rule_overlay import GeneRuleEvidence, HIGH_CONFIDENCE, LOW_CONFIDENCE
from pcsec_pichia.services.target_gpr_overlay_review import (
    CANDIDATE_GPR_OVERLAY_REVIEW,
    MANUAL_REVIEW_REQUIRED,
    MODEL_EXPANSION_REQUIRED,
    build_hlf_opn_gpr_overlay_review_rows,
    filter_hlf_opn_gpr_overlay_review_rows,
    load_hlf_opn_gpr_overlay_review_cache,
    summarize_hlf_opn_gpr_overlay_review_rows,
    write_hlf_opn_gpr_overlay_review_cache,
)


def test_high_confidence_model_external_candidate_enters_overlay_review(tmp_path) -> None:
    rows = build_hlf_opn_gpr_overlay_review_rows(
        candidate_rows=_candidate_rows(),
        gene_rule_evidence_by_name=_evidence(),
        model_reaction_ids={"sec_PDI1_ERV2_Ero1p_complex_formation"},
    )
    output = tmp_path / "overlay_review.json"
    write_hlf_opn_gpr_overlay_review_cache(rows, output)
    loaded = load_hlf_opn_gpr_overlay_review_cache(output)
    summary = summarize_hlf_opn_gpr_overlay_review_rows(loaded)
    pdi = next(row for row in loaded if row.source_common_name == "PDI1")

    assert loaded == rows
    assert pdi.review_status == CANDIDATE_GPR_OVERLAY_REVIEW
    assert pdi.gene_id == "PAS_PDI1_EXTERNAL"
    assert pdi.source_candidate_gene_id == "PAS_PDI1_HOMOLOGY"
    assert pdi.existing_model_reaction_ids == ("sec_PDI1_ERV2_Ero1p_complex_formation",)
    assert pdi.external_ids["uniprot"] == "U-PDI1"
    assert "formal GPR" in pdi.warnings[0]
    assert "not phenotype evidence" in " ".join(pdi.warnings)
    assert summary["candidate_gpr_overlay_review_count"] == 2


def test_low_confidence_or_incomplete_evidence_stays_manual_review() -> None:
    rows = build_hlf_opn_gpr_overlay_review_rows(
        candidate_rows=_candidate_rows(),
        gene_rule_evidence_by_name=_evidence(),
        model_reaction_ids={"sec_PDI1_ERV2_Ero1p_complex_formation"},
    )
    erv2 = next(row for row in rows if row.source_common_name == "ERV2")

    assert erv2.review_status == MANUAL_REVIEW_REQUIRED
    assert erv2.evidence_confidence == LOW_CONFIDENCE
    assert erv2.existing_model_reaction_ids == ("sec_PDI1_ERV2_Ero1p_complex_formation",)
    assert any("manual_review_required" in warning for warning in erv2.warnings)


def test_high_confidence_without_existing_reaction_requires_model_expansion() -> None:
    rows = build_hlf_opn_gpr_overlay_review_rows(
        candidate_rows=_candidate_rows(),
        gene_rule_evidence_by_name=_evidence(),
        model_reaction_ids={"sec_PDI1_ERV2_Ero1p_complex_formation"},
    )
    gap = next(row for row in rows if row.source_common_name == "GAP1")

    assert gap.review_status == MODEL_EXPANSION_REQUIRED
    assert gap.missing_model_reaction_ids == ("sec_missing_gap_reaction",)
    assert any("No current model reaction" in warning for warning in gap.warnings)


def test_overlay_review_filters_hlf_opn_and_skips_model_operable_candidates() -> None:
    rows = build_hlf_opn_gpr_overlay_review_rows(
        candidate_rows=_candidate_rows(),
        gene_rule_evidence_by_name=_evidence(),
        model_reaction_ids={"sec_PDI1_ERV2_Ero1p_complex_formation"},
    )
    hlf_rows = filter_hlf_opn_gpr_overlay_review_rows(rows, target_context="hLF")
    opn_rows = filter_hlf_opn_gpr_overlay_review_rows(rows, target_context="OPN")

    assert {row.source_common_name for row in hlf_rows} >= {"PDI1", "ERV2", "GAP1"}
    assert {row.source_common_name for row in opn_rows} == {"PMT9"}
    assert "KAR2 / BiP" not in {row.source_common_name for row in rows}


def _candidate_rows() -> list[dict[str, object]]:
    return [
        {
            "target_context": "hLF",
            "source_common_name": "PDI1",
            "gene_id": "PAS_PDI1_HOMOLOGY",
            "candidate_role": "disulfide_bond_folding",
            "recommended_intervention": "OE",
            "model_operable": False,
            "homology_review_status": "rbh_not_in_model",
            "rule_transfer_status": "rule_transfer_supported_not_model_operable",
        },
        {
            "target_context": "hLF",
            "source_common_name": "ERV2",
            "gene_id": "PAS_ERV2_HOMOLOGY",
            "candidate_role": "disulfide_bond_folding",
            "recommended_intervention": "OE",
            "model_operable": False,
        },
        {
            "target_context": "hLF",
            "source_common_name": "GAP1",
            "gene_id": "PAS_GAP1",
            "candidate_role": "manual_review",
            "recommended_intervention": "OE",
            "model_operable": False,
        },
        {
            "target_context": "OPN",
            "source_common_name": "PMT9",
            "gene_id": "PAS_PMT9",
            "candidate_role": "o_glycosylation_processing",
            "recommended_intervention": "OE",
            "model_operable": False,
        },
        {
            "target_context": "shared",
            "source_common_name": "KAR2 / BiP",
            "gene_id": "PAS_KAR2",
            "candidate_role": "er_folding_chaperone",
            "recommended_intervention": "OE",
            "model_operable": True,
        },
    ]


def _evidence() -> dict[str, GeneRuleEvidence]:
    return {
        "PDI1": GeneRuleEvidence(
            common_name="PDI1",
            candidate_locus_tag="PAS_PDI1_EXTERNAL",
            external_ids={"uniprot": "U-PDI1", "kegg": "ppa:PAS_PDI1_EXTERNAL"},
            protein_name="protein disulfide-isomerase",
            evidence_sources=("UniProt GS115 proteome exact locus", "KEGG ppa exact locus"),
            confidence=HIGH_CONFIDENCE,
            target_reaction_ids=("sec_PDI1_ERV2_Ero1p_complex_formation",),
            rule_status="high_confidence_locus_candidate",
            recommended_action="eligible_for_overlay_if_all_complex_subunits_are_confirmed",
        ),
        "ERV2": GeneRuleEvidence(
            common_name="ERV2",
            candidate_locus_tag="PAS_ERV2_EXTERNAL",
            evidence_sources=("KEGG ppa name search",),
            confidence=LOW_CONFIDENCE,
            target_reaction_ids=("sec_PDI1_ERV2_Ero1p_complex_formation",),
            rule_status="not_executable",
            recommended_action="manual_locus_review_required",
        ),
        "GAP1": GeneRuleEvidence(
            common_name="GAP1",
            candidate_locus_tag="PAS_GAP1",
            evidence_sources=("UniProt GS115 proteome exact locus", "KEGG ppa exact locus"),
            confidence=HIGH_CONFIDENCE,
            target_reaction_ids=("sec_missing_gap_reaction",),
            rule_status="high_confidence_locus_candidate",
            recommended_action="eligible_for_overlay_if_model_reaction_is_added",
        ),
        "PMT9": GeneRuleEvidence(
            common_name="PMT9",
            candidate_locus_tag="PAS_PMT9",
            evidence_sources=("UniProt GS115 proteome exact locus", "KEGG ppa exact locus"),
            confidence=HIGH_CONFIDENCE,
            target_reaction_ids=("sec_PDI1_ERV2_Ero1p_complex_formation",),
            rule_status="high_confidence_locus_candidate",
            recommended_action="eligible_for_overlay_if_model_reaction_is_added",
        ),
        "KAR2": GeneRuleEvidence(
            common_name="KAR2",
            candidate_locus_tag="PAS_KAR2_EXTERNAL",
            confidence=HIGH_CONFIDENCE,
            target_reaction_ids=("sec_PDI1_ERV2_Ero1p_complex_formation",),
        ),
    }
