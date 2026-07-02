from __future__ import annotations

from pcsec_pichia.analysis import (
    analyze_yield_improvement_candidates,
    summarize_yield_improvement_recommendations,
)


CURATED_ROWS = (
    {
        "common_name": "PEP4",
        "mapped_model_gene_id": "PAS_chr2-2_0107",
        "declared_model_gene_id": "PAS_chr2-2_0107",
        "mapping_status": "model_gpr_gene_available",
        "recommended_use": "gene_level_gpr_perturbation",
        "gene_level_ready": True,
    },
    {
        "common_name": "PRB1",
        "mapped_model_gene_id": "PAS_chr2-1_0785",
        "declared_model_gene_id": "PAS_chr2-1_0785",
        "mapping_status": "model_gpr_gene_available",
        "recommended_use": "gene_level_gpr_perturbation",
        "gene_level_ready": True,
    },
    {
        "common_name": "PDI1",
        "mapped_model_gene_id": "",
        "declared_model_gene_id": "",
        "mapping_status": "reaction_proxy_only",
        "recommended_use": "reaction_level_proxy_requires_locus_review",
        "oe_reaction_id": "sec_PDI1_ERV2_Ero1p_complex_formation",
        "reaction_proxy_ready": True,
    },
    {
        "common_name": "ERO1",
        "mapped_model_gene_id": "",
        "declared_model_gene_id": "",
        "mapping_status": "reaction_proxy_only",
        "recommended_use": "reaction_level_proxy_requires_locus_review",
        "oe_reaction_id": "sec_PDI1_ERV2_Ero1p_complex_formation",
        "reaction_proxy_ready": True,
    },
)


def test_yield_recommendations_prioritize_gene_level_ko_over_proxy() -> None:
    result = analyze_yield_improvement_candidates(
        [
            {
                "target_id": "OPN_ALPHA_FULL_PROJECT",
                "candidate_id": "PAS_chr2-2_0107",
                "gene_id": "PAS_chr2-2_0107",
                "canonical_gene_id": "PAS_chr2-2_0107",
                "input_gene_id": "PEP4",
                "intervention_type": "KO",
                "ko_support_status": "ko_runnable_gpr_gene_deletion",
                "success": True,
                "status": "optimal",
                "delta_objective": 0.00002,
                "effect_label": "提升分泌",
                "secretory_process": "proteasome",
            },
            {
                "target_id": "OPN_ALPHA_FULL_PROJECT",
                "candidate_id": "sec_PDI1_ERV2_Ero1p_complex_formation",
                "reaction_id": "sec_PDI1_ERV2_Ero1p_complex_formation",
                "input_gene_id": "PDI1",
                "intervention_type": "OE_gene_proxy",
                "success": True,
                "status": "optimal",
                "delta_objective": 0.00005,
                "effect_label": "提升分泌",
                "secretory_process": "folding_dsb",
            },
        ],
        target_id="OPN_ALPHA_FULL_PROJECT",
        curated_gene_evidence=CURATED_ROWS,
    )
    summary = summarize_yield_improvement_recommendations(result)
    recommended = summary["recommended_candidates"]

    assert recommended[0]["display_name"] == "PEP4"
    assert recommended[0]["recommendation_label"] == "strong_model_candidate"
    assert recommended[0]["execution_mode"] == "gene_level_ko"
    assert recommended[1]["display_name"] == "PDI1"
    assert recommended[1]["recommendation_label"] == "promising_but_proxy_only"
    assert recommended[1]["execution_mode"] == "reaction_level_oe_proxy"
    assert "reaction-level proxy" in recommended[1]["warnings"][0]


def test_yield_recommendations_do_not_mark_pdi_ero_as_gene_level_ko() -> None:
    result = analyze_yield_improvement_candidates(
        [
            {
                "target_id": "hLF",
                "candidate_id": "sec_PDI1_ERV2_Ero1p_complex_formation",
                "reaction_id": "sec_PDI1_ERV2_Ero1p_complex_formation",
                "input_gene_id": "ERO1",
                "intervention_type": "OE_gene_proxy",
                "success": True,
                "status": "optimal",
                "delta_objective": 0.00001,
                "effect_label": "提升分泌",
                "secretory_process": "folding_dsb",
            }
        ],
        target_id="hLF",
        curated_gene_evidence=CURATED_ROWS,
    )

    row = summarize_yield_improvement_recommendations(result)["recommended_candidates"][0]

    assert row["display_name"] == "ERO1"
    assert row["execution_mode"] == "reaction_level_oe_proxy"
    assert row["recommendation_label"] == "promising_but_proxy_only"
    assert row["model_gene_id"] == ""


def test_yield_recommendations_exclude_failed_and_unresolved_from_top_list() -> None:
    result = analyze_yield_improvement_candidates(
        [
            {
                "target_id": "OPN_ALPHA_FULL_PROJECT",
                "candidate_id": "PRB1",
                "gene_id": "PAS_chr2-1_0785",
                "input_gene_id": "PRB1",
                "intervention_type": "KO",
                "ko_support_status": "ko_runnable_gpr_gene_deletion",
                "success": False,
                "status": "2",
                "effect_label": "约束不可行",
                "delta_objective": "",
            },
            {
                "target_id": "OPN_ALPHA_FULL_PROJECT",
                "candidate_id": "NO_SUCH",
                "input_gene_id": "NO_SUCH",
                "intervention_type": "KO",
                "success": False,
                "status": "unresolved_gene",
                "mapping_level": "unresolved",
                "effect_label": "未解析",
            },
        ],
        target_id="OPN_ALPHA_FULL_PROJECT",
        curated_gene_evidence=CURATED_ROWS,
    )
    summary = summarize_yield_improvement_recommendations(result)

    assert summary["recommended_candidates"] == ()
    assert summary["not_recommended_candidates"][0]["recommendation_label"] == "not_recommended_solver_failed"
    assert summary["unresolved_candidates"][0]["recommendation_label"] == "unresolved_not_actionable"
    assert summary["summary_counts"] == {"recommended": 0, "not_recommended": 1, "unresolved": 1}


def test_yield_recommendations_mark_essential_ko_as_growth_risk_without_oe_inheritance() -> None:
    result = analyze_yield_improvement_candidates(
        [
            {
                "target_id": "OPN_ALPHA_FULL_PROJECT",
                "candidate_id": "ESSENTIAL_TEST",
                "gene_id": "ESSENTIAL_TEST",
                "input_gene_id": "ESSENTIAL_TEST",
                "intervention_type": "KO",
                "ko_support_status": "ko_runnable_gpr_gene_deletion",
                "success": True,
                "status": "optimal",
                "delta_objective": 0.0001,
            },
            {
                "target_id": "OPN_ALPHA_FULL_PROJECT",
                "candidate_id": "ESSENTIAL_TEST",
                "gene_id": "ESSENTIAL_TEST",
                "input_gene_id": "ESSENTIAL_TEST",
                "intervention_type": "OE_gene_proxy",
                "success": True,
                "status": "optimal",
                "delta_objective": 0.0001,
            },
        ],
        target_id="OPN_ALPHA_FULL_PROJECT",
    )
    summary = summarize_yield_improvement_recommendations(result)

    assert summary["not_recommended_candidates"][0]["recommendation_tier"] == "not_recommended_growth_risk"
    assert summary["not_recommended_candidates"][0]["recommendation_label"] == "not_recommended_growth_risk"
    assert summary["not_recommended_candidates"][0]["growth_risk_label"] == "essential_gene_ko"
    assert summary["recommended_candidates"][0]["intervention_type"] == "OE_gene_proxy"
    assert summary["recommended_candidates"][0]["recommendation_tier"] == "model_executable"


def test_yield_recommendation_tiers_require_source_intervention_and_context_match() -> None:
    result = analyze_yield_improvement_candidates(
        [
            {
                "target_id": "OPN_ALPHA_FULL_PROJECT",
                "candidate_id": "POSITIVE_SCREEN_TEST",
                "gene_id": "POSITIVE_SCREEN_TEST",
                "input_gene_id": "POSITIVE_SCREEN_TEST",
                "intervention_type": "OE_gene_proxy",
                "success": True,
                "status": "optimal",
                "delta_objective": 0.0001,
            },
            {
                "target_id": "OPN_ALPHA_FULL_PROJECT",
                "candidate_id": "CALIBRATED_TEST",
                "gene_id": "CALIBRATED_TEST",
                "input_gene_id": "CALIBRATED_TEST",
                "intervention_type": "OE_gene_proxy",
                "success": True,
                "status": "optimal",
                "delta_objective": 0.00009,
            },
        ],
        target_id="OPN_ALPHA_FULL_PROJECT",
    )
    rows = summarize_yield_improvement_recommendations(result)["recommended_candidates"]
    by_gene = {row["display_name"]: row for row in rows}

    assert by_gene["POSITIVE_SCREEN_TEST"]["recommendation_tier"] == "evidence_supported"
    assert by_gene["CALIBRATED_TEST"]["recommendation_tier"] == "experiment_calibrated"
    assert by_gene["CALIBRATED_TEST"]["phenotype_evidence"]["evidence_source"] == (
        "same_target_same_host_same_intervention"
    )


def test_yield_recommendation_context_mismatch_caps_experiment_evidence() -> None:
    result = analyze_yield_improvement_candidates(
        [
            {
                "target_id": "hLF",
                "candidate_id": "CALIBRATED_TEST",
                "gene_id": "CALIBRATED_TEST",
                "input_gene_id": "CALIBRATED_TEST",
                "intervention_type": "OE_gene_proxy",
                "success": True,
                "status": "optimal",
                "delta_objective": 0.0001,
            }
        ],
        target_id="hLF",
    )
    row = summarize_yield_improvement_recommendations(result)["recommended_candidates"][0]

    assert row["recommendation_tier"] == "evidence_supported"


def test_yield_recommendation_unresolved_candidate_stays_manual_review_despite_phenotype_evidence() -> None:
    result = analyze_yield_improvement_candidates(
        [
            {
                "target_id": "OPN_ALPHA_FULL_PROJECT",
                "candidate_id": "POSITIVE_SCREEN_TEST",
                "gene_id": "POSITIVE_SCREEN_TEST",
                "input_gene_id": "POSITIVE_SCREEN_TEST",
                "intervention_type": "OE_gene_proxy",
                "success": False,
                "status": "unresolved_gene",
                "delta_objective": "",
            }
        ],
        target_id="OPN_ALPHA_FULL_PROJECT",
    )
    row = summarize_yield_improvement_recommendations(result)["unresolved_candidates"][0]

    assert row["recommendation_tier"] == "manual_review_required"
    assert row["phenotype_evidence"]["secretion_screen_effect"] == "screen_positive"
    assert row["recommendation_tier"] != "evidence_supported"


def test_yield_recommendation_annotation_only_is_not_experiment_calibrated() -> None:
    result = analyze_yield_improvement_candidates(
        [
            {
                "target_id": "OPN_ALPHA_FULL_PROJECT",
                "candidate_id": "ANNOTATION_ONLY",
                "gene_id": "ANNOTATION_ONLY",
                "input_gene_id": "ANNOTATION_ONLY",
                "intervention_type": "KO",
                "ko_support_status": "ko_runnable_gpr_gene_deletion",
                "success": True,
                "status": "optimal",
                "delta_objective": 0.0001,
                "evidence_sources": ["UniProt", "KEGG"],
                "evidence_confidence": "high_exact_locus_tag",
            }
        ],
        target_id="OPN_ALPHA_FULL_PROJECT",
    )
    row = summarize_yield_improvement_recommendations(result)["recommended_candidates"][0]

    assert row["recommendation_tier"] == "model_executable"
    assert row["recommendation_tier"] != "experiment_calibrated"
    assert row["database_annotation_sources"] == ("UniProt", "KEGG")


def test_yield_recommendation_not_run_oe_proxy_is_manual_review() -> None:
    result = analyze_yield_improvement_candidates(
        [
            {
                "target_id": "OPN_ALPHA_FULL_PROJECT",
                "candidate_id": "NO_GPR_OE",
                "gene_id": "NO_GPR_OE",
                "input_gene_id": "NO_GPR_OE",
                "intervention_type": "OE_gene_proxy",
                "success": False,
                "status": "not_run_no_gpr_effect",
                "delta_objective": "",
            }
        ],
        target_id="OPN_ALPHA_FULL_PROJECT",
    )
    row = summarize_yield_improvement_recommendations(result)["not_recommended_candidates"][0]

    assert row["recommendation_tier"] == "manual_review_required"
    assert row["oe_reaction_proxy"] is False
