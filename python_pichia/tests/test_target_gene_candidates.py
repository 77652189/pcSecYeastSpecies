from __future__ import annotations

from pcsec_pichia.services.gene_catalog import SecretionGeneEntry
from pcsec_pichia.services.homology_evidence import GeneHomologyEvidence
from pcsec_pichia.services.target_gene_candidates import (
    MODEL_KO_EXECUTABLE,
    MODEL_OE_PROXY_EXECUTABLE,
    NOT_IN_MODEL,
    UNRESOLVED_NAME,
    build_hlf_opn_candidate_gene_rows,
    executable_inputs_for_hlf_opn_candidates,
    filter_hlf_opn_candidate_gene_rows,
    load_hlf_opn_candidate_gene_cache,
    summarize_hlf_opn_candidate_gene_rows,
    write_hlf_opn_candidate_gene_cache,
)


def test_hlf_opn_candidates_filter_contexts_and_keep_required_fields(tmp_path) -> None:
    rows = build_hlf_opn_candidate_gene_rows(
        full_model_rows=_full_model_rows(),
        standard_name_rows=_standard_rows(),
        secretion_gene_catalog=_catalog_entries(),
        homology_evidence_by_gene=_homology_evidence(),
    )

    hlf_rows = filter_hlf_opn_candidate_gene_rows(rows, target_context="hLF")
    opn_rows = filter_hlf_opn_candidate_gene_rows(rows, target_context="OPN")
    summary = summarize_hlf_opn_candidate_gene_rows(rows)
    cache_path = tmp_path / "hlf_opn_candidate_genes.json"
    write_hlf_opn_candidate_gene_cache(rows, cache_path)
    loaded = load_hlf_opn_candidate_gene_cache(cache_path)

    assert loaded == rows
    assert summary["target_candidate_counts"]["hLF"] >= 2
    assert summary["target_candidate_counts"]["OPN"] >= 2
    assert {row.target_context for row in rows} == {"hLF", "OPN", "shared"}
    assert {row.gene_id for row in hlf_rows} >= {"PAS_KAR2", "PAS_PDI1"}
    assert {row.gene_id for row in opn_rows} >= {"PAS_PMT1", "PAS_KAR2"}
    required = {
        "target_context",
        "gene_id",
        "candidate_role",
        "evidence_type",
        "evidence_confidence",
        "model_operable",
        "recommended_intervention",
        "reason",
        "warnings",
    }
    assert required <= set(rows[0].to_dict())


def test_operability_audit_separates_executable_and_review_only_candidates() -> None:
    rows = build_hlf_opn_candidate_gene_rows(
        full_model_rows=_full_model_rows(),
        standard_name_rows=_standard_rows(),
        secretion_gene_catalog=_catalog_entries(),
        homology_evidence_by_gene=_homology_evidence(),
    )
    by_name = {row.source_common_name: row for row in rows}

    assert by_name["KAR2 / BiP"].operability_status == MODEL_OE_PROXY_EXECUTABLE
    assert by_name["KAR2 / BiP"].model_operable is True
    assert any("reaction-level proxy" in warning for warning in by_name["KAR2 / BiP"].warnings)
    assert by_name["PEP4"].operability_status == MODEL_KO_EXECUTABLE
    assert by_name["PDI1"].operability_status == NOT_IN_MODEL
    assert by_name["PDI1"].model_operability_label == "not_model_operable"
    assert by_name["UNRESOLVED"].operability_status == UNRESOLVED_NAME
    assert by_name["UNRESOLVED"].model_operability_label == "not_model_operable"

    executable = executable_inputs_for_hlf_opn_candidates(rows, target_context="hLF")

    assert executable["ko_gene_ids"] == ["PAS_PEP4"]
    assert executable["oe_gene_ids"] == ["PAS_KAR2"]
    assert "PAS_PDI1" not in executable["oe_gene_ids"]
    assert executable["excluded_count"] >= 2


def test_homology_is_auxiliary_not_phenotype_evidence() -> None:
    rows = build_hlf_opn_candidate_gene_rows(
        full_model_rows=_full_model_rows(),
        standard_name_rows=_standard_rows(),
        secretion_gene_catalog=_catalog_entries(),
        homology_evidence_by_gene=_homology_evidence(),
    )
    pdi = next(row for row in rows if row.source_common_name == "PDI1")

    assert pdi.evidence_type == "homology_auxiliary"
    assert pdi.homology_review_status == "rbh_not_in_model"
    assert all("phenotype" not in warning.lower() or "not phenotype evidence" in warning for warning in pdi.warnings)


def _catalog_entries() -> tuple[SecretionGeneEntry, ...]:
    return (
        SecretionGeneEntry(
            category="ER 折叠与分子伴侣",
            common_name="KAR2 / BiP",
            description="ER chaperone",
            intervention="OE",
            oe_reaction_id="sec_Kar2p_complex_formation",
            evidence="已报道 Kar2 过表达可提升毕赤酵母外源蛋白分泌",
        ),
        SecretionGeneEntry(
            category="二硫键 (DSB)",
            common_name="PDI1",
            description="disulfide folding",
            intervention="OE",
            oe_reaction_id="sec_Pdi1p_complex_formation",
            evidence="模型 sec_Pdi1p 复合体",
        ),
        SecretionGeneEntry(
            category="O-糖基化",
            common_name="PMT1/PMT2/PMT4-6",
            description="O-glycosylation",
            intervention="OE",
            oe_reaction_id="sec_Pmt_complex_formation",
            evidence="模型 sec_Pmt 复合体",
        ),
        SecretionGeneEntry(
            category="蛋白酶体与降解",
            common_name="PEP4",
            description="protease",
            gene_id="PAS_PEP4",
            intervention="KO",
            evidence="毕赤酵母蛋白表达常用 KO 靶点",
        ),
        SecretionGeneEntry(
            category="ER 转运",
            common_name="UNRESOLVED",
            description="review",
            intervention="OE",
            oe_reaction_id="sec_missing_complex",
            evidence="模型 proxy",
        ),
    )


def _full_model_rows() -> list[dict[str, object]]:
    return [
        {
            "gene_id": "PAS_KAR2",
            "display_name": "KAR2",
            "ko_support_status": "ko_runnable_gpr_gene_deletion",
            "oe_support_status": "oe_runnable_reaction_proxy",
            "affected_reactions": ["NTP1er_no_1_fwd"],
            "oe_executable_reactions": ["NTP1er_no_1_fwd"],
            "inactive_reactions_if_ko": ["NTP1er_no_1_fwd"],
            "evidence_sources": ["UniProt"],
            "evidence_confidence": "high_exact_locus_tag",
        },
        {
            "gene_id": "PAS_PMT1",
            "display_name": "PMT1",
            "ko_support_status": "ko_runnable_gpr_gene_deletion",
            "oe_support_status": "oe_runnable_reaction_proxy",
            "affected_reactions": ["OGLYCOS_no_1_fwd"],
            "oe_executable_reactions": ["OGLYCOS_no_1_fwd"],
            "evidence_sources": ["UniProt"],
            "evidence_confidence": "high_exact_locus_tag",
        },
        {
            "gene_id": "PAS_PEP4",
            "display_name": "PEP4 curated model row",
            "ko_support_status": "ko_runnable_gpr_gene_deletion",
            "oe_support_status": "oe_runnable_reaction_proxy",
            "affected_reactions": ["VAC_PROTEASE_no_1_fwd"],
            "inactive_reactions_if_ko": ["VAC_PROTEASE_no_1_fwd"],
            "evidence_sources": ["UniProt"],
            "evidence_confidence": "high_exact_locus_tag",
        },
    ]


def _standard_rows() -> list[dict[str, object]]:
    return [
        {
            "gene_id": "PAS_KAR2",
            "display_name": "KAR2",
            "standard_symbol": "KAR2",
            "protein_name": "ER chaperone BiP",
            "external_ids": {"uniprot": "U-KAR2"},
            "annotation_sources": ["UniProt", "KEGG"],
            "annotation_confidence": "high_exact_locus_tag",
            "model_operable": True,
            "gpr_status": "ko_and_oe_model_executable",
        },
        {
            "gene_id": "PAS_PMT1",
            "display_name": "PMT1",
            "standard_symbol": "PMT1",
            "protein_name": "O-mannosyltransferase",
            "external_ids": {"uniprot": "U-PMT1"},
            "annotation_sources": ["UniProt", "KEGG"],
            "annotation_confidence": "high_exact_locus_tag",
            "model_operable": True,
            "gpr_status": "ko_and_oe_model_executable",
        },
        {
            "gene_id": "PAS_PEP4",
            "display_name": "PEP4 curated model row",
            "external_ids": {"uniprot": "U-PEP4"},
            "annotation_sources": ["UniProt"],
            "annotation_confidence": "high_exact_locus_tag",
            "model_operable": True,
            "gpr_status": "ko_and_oe_model_executable",
        },
    ]


def _homology_evidence() -> dict[str, GeneHomologyEvidence]:
    return {
        "kar2": GeneHomologyEvidence(
            gene_id="PAS_KAR2",
            internal_common_name="KAR2 / BiP",
            query_symbol="KAR2",
            pichia_gene_id="PAS_KAR2",
            pichia_model_gene_id="PAS_KAR2",
            is_rbh=True,
            in_model_gene_index=True,
            homology_review_status="model_ready_rbh_high_confidence",
            rule_transfer_status="rule_transfer_ready",
        ),
        "pdi1": GeneHomologyEvidence(
            gene_id="PAS_PDI1",
            internal_common_name="PDI1",
            query_symbol="PDI1",
            pichia_gene_id="PAS_PDI1",
            pichia_model_gene_id="",
            is_rbh=True,
            in_model_gene_index=False,
            homology_review_status="rbh_not_in_model",
            rule_transfer_status="rule_transfer_supported_not_model_operable",
        ),
        "pmt1": GeneHomologyEvidence(
            gene_id="PAS_PMT1",
            internal_common_name="PMT1/PMT2/PMT4-6",
            query_symbol="PMT1",
            pichia_gene_id="PAS_PMT1",
            pichia_model_gene_id="PAS_PMT1",
            is_rbh=True,
            in_model_gene_index=True,
            homology_review_status="model_ready_rbh_high_confidence",
            rule_transfer_status="rule_transfer_ready",
        ),
        "pmt1/pmt2/pmt4-6": GeneHomologyEvidence(
            gene_id="PAS_PMT1",
            internal_common_name="PMT1/PMT2/PMT4-6",
            query_symbol="PMT1",
            pichia_gene_id="PAS_PMT1",
            pichia_model_gene_id="PAS_PMT1",
            is_rbh=True,
            in_model_gene_index=True,
            homology_review_status="model_ready_rbh_high_confidence",
            rule_transfer_status="rule_transfer_ready",
        ),
    }
