from __future__ import annotations

import json
from pathlib import Path

from pcsec_pichia.screens.candidate_resolution import reactions_for_gene
from pcsec_pichia.screens.planning import build_screen_plan
from pcsec_pichia.services.gene_rule_overlay import (
    DEFAULT_GENE_RULE_EVIDENCE_REPORT,
    HIGH_CONFIDENCE,
    LOW_CONFIDENCE,
    GeneRuleEvidence,
    apply_gpr_overlay_for_analysis,
    build_gene_rule_evidence_cache,
    build_gpr_overlay,
    load_gene_rule_evidence_cache,
    overlay_aliases_for_executable_rules,
)
from pcsec_pichia.services import gene_rule_overlay


class TinyModel:
    rxns = ["sec_PDI1_ERV2_Ero1p_complex_formation", "R_OTHER"]
    rules = ["", ""]
    gr_rules = ["", ""]
    genes = ["G_ORIGINAL"]
    gene_index = {"G_ORIGINAL": 0}
    reaction_index = {"sec_PDI1_ERV2_Ero1p_complex_formation": 0, "R_OTHER": 1}


def test_high_confidence_pdi_complex_overlay_adds_gene_level_gpr_rules() -> None:
    evidence = {
        "PDI1": _high("PDI1", "PAS_PDI1"),
        "ERO1": _high("ERO1", "PAS_ERO1"),
        "ERV2": _high("ERV2", "PAS_ERV2"),
    }

    overlay = build_gpr_overlay(TinyModel(), evidence)
    overlaid = apply_gpr_overlay_for_analysis(TinyModel(), overlay)

    assert overlay.entries[0].reaction_id == "sec_PDI1_ERV2_Ero1p_complex_formation"
    assert reactions_for_gene(overlaid, "PAS_PDI1") == ["sec_PDI1_ERV2_Ero1p_complex_formation"]
    assert reactions_for_gene(overlaid, "PAS_ERO1") == ["sec_PDI1_ERV2_Ero1p_complex_formation"]
    assert overlaid.rules[0] == "x(2) & x(3) & x(4)"
    assert overlaid.gr_rules[0] == "PAS_PDI1 and PAS_ERO1 and PAS_ERV2"


def test_overlay_common_name_inputs_resolve_to_locus_tags_for_screen_plan() -> None:
    evidence = {
        "PDI1": _high("PDI1", "PAS_PDI1"),
        "ERO1": _high("ERO1", "PAS_ERO1"),
        "ERV2": _high("ERV2", "PAS_ERV2"),
    }
    overlay = build_gpr_overlay(TinyModel(), evidence)
    overlaid = apply_gpr_overlay_for_analysis(TinyModel(), overlay)

    class Request:
        screen_candidate_limit = 20
        ko_gene_ids = ("PDI1",)
        ko_candidates = ()
        ko_reaction_ids = ()
        oe_gene_ids = ("ERO1",)
        oe_reaction_ids = ()
        oe_candidates = ()
        enable_gene_rule_overlay = True

    screen_plan = build_screen_plan(
        overlaid,
        Request(),
        gene_rule_evidence_by_name=evidence,
        gene_rule_overlay=overlay,
    )

    assert screen_plan.ko_gene_ids == ["PAS_PDI1"]
    assert screen_plan.ko_input_by_gene["PAS_PDI1"] == "PDI1"
    assert screen_plan.oe_gene_plans_by_gene["PAS_ERO1"].affected_reactions == (
        "sec_PDI1_ERV2_Ero1p_complex_formation",
    )
    assert screen_plan.oe_input_by_gene["PAS_ERO1"] == "ERO1"


def test_low_confidence_or_incomplete_overlay_is_display_only() -> None:
    evidence = {
        "PDI1": _high("PDI1", "PAS_PDI1"),
        "ERO1": GeneRuleEvidence(
            common_name="ERO1",
            candidate_locus_tag="PAS_ERO1",
            confidence=LOW_CONFIDENCE,
            target_reaction_ids=("sec_PDI1_ERV2_Ero1p_complex_formation",),
        ),
    }

    overlay = build_gpr_overlay(TinyModel(), evidence)
    overlaid = apply_gpr_overlay_for_analysis(TinyModel(), overlay)

    assert overlay.entries == ()
    assert overlay.skipped
    assert reactions_for_gene(overlaid, "PAS_PDI1") == []


def test_pdi_complex_overlay_requires_each_subunit_to_target_the_same_reaction() -> None:
    evidence = {
        "PDI1": _high("PDI1", "PAS_PDI1"),
        "ERO1": GeneRuleEvidence(
            common_name="ERO1",
            candidate_locus_tag="PAS_ERO1",
            confidence=HIGH_CONFIDENCE,
            target_reaction_ids=("R_OTHER",),
        ),
        "ERV2": _high("ERV2", "PAS_ERV2"),
    }

    overlay = build_gpr_overlay(TinyModel(), evidence)

    assert overlay.entries == ()
    assert overlay.skipped[0]["reason"] == "shared_complex_subunit_targets_different_reaction"


def test_pdi_complex_overlay_normalizes_common_name_case() -> None:
    incomplete = {
        "pdi1": GeneRuleEvidence(
            common_name="pdi1",
            candidate_locus_tag="PAS_PDI1",
            confidence=HIGH_CONFIDENCE,
            target_reaction_ids=("sec_PDI1_ERV2_Ero1p_complex_formation",),
        )
    }

    incomplete_overlay = build_gpr_overlay(TinyModel(), incomplete)

    assert incomplete_overlay.entries == ()
    assert incomplete_overlay.skipped[0]["reason"] == "shared_complex_requires_locus_for_all_subunits"

    complete = {
        "pdi1": _high("pdi1", "PAS_PDI1"),
        "ero1": _high("ero1", "PAS_ERO1"),
        "erv2": _high("erv2", "PAS_ERV2"),
    }

    complete_overlay = build_gpr_overlay(TinyModel(), complete)

    assert complete_overlay.entries[0].common_names == ("PDI1", "ERO1", "ERV2")
    assert complete_overlay.entries[0].gene_locus_tags == ("PAS_PDI1", "PAS_ERO1", "PAS_ERV2")


def test_generic_overlay_cannot_write_protected_pdi_complex_reaction() -> None:
    evidence = {
        "PDI_ALT": GeneRuleEvidence(
            common_name="PDI_ALT",
            candidate_locus_tag="PAS_ALT",
            confidence=HIGH_CONFIDENCE,
            target_reaction_ids=("sec_PDI1_ERV2_Ero1p_complex_formation",),
        )
    }

    overlay = build_gpr_overlay(TinyModel(), evidence)

    assert overlay.entries == ()
    assert overlay.skipped[0]["reason"] == "protected_complex_reaction_requires_complex_rule"


def test_overlay_aliases_require_actual_overlay_entries() -> None:
    evidence = {
        "PDI_ALT": GeneRuleEvidence(
            common_name="PDI_ALT",
            candidate_locus_tag="PAS_ALT",
            confidence=HIGH_CONFIDENCE,
            target_reaction_ids=("sec_PDI1_ERV2_Ero1p_complex_formation",),
        )
    }
    overlay = build_gpr_overlay(TinyModel(), evidence)

    assert overlay.entries == ()
    assert overlay_aliases_for_executable_rules(evidence, overlay) == {}
    assert overlay_aliases_for_executable_rules(evidence, None) == {}


def test_overlay_does_not_replace_existing_model_gpr_rule() -> None:
    class ExistingRuleModel:
        rxns = ["sec_PDI1_ERV2_Ero1p_complex_formation"]
        rules = ["x(1)"]
        gr_rules = ["G_ORIGINAL"]
        genes = ["G_ORIGINAL"]
        gene_index = {"G_ORIGINAL": 0}
        reaction_index = {"sec_PDI1_ERV2_Ero1p_complex_formation": 0}

    evidence = {
        "PDI1": _high("PDI1", "PAS_PDI1"),
        "ERO1": _high("ERO1", "PAS_ERO1"),
        "ERV2": _high("ERV2", "PAS_ERV2"),
    }

    overlay = build_gpr_overlay(ExistingRuleModel(), evidence)
    overlaid = apply_gpr_overlay_for_analysis(ExistingRuleModel(), overlay)

    assert overlay.entries == ()
    assert overlay.skipped[0]["reason"] == "target_reaction_already_has_model_gpr_rule"
    assert overlaid.rules == ["x(1)"]
    assert overlaid.gr_rules == ["G_ORIGINAL"]


def test_overlay_does_not_mutate_original_model() -> None:
    model = TinyModel()
    evidence = {
        "PDI1": _high("PDI1", "PAS_PDI1"),
        "ERO1": _high("ERO1", "PAS_ERO1"),
        "ERV2": _high("ERV2", "PAS_ERV2"),
    }

    overlay = build_gpr_overlay(model, evidence)
    overlaid = apply_gpr_overlay_for_analysis(model, overlay)

    assert model.genes == ["G_ORIGINAL"]
    assert model.rules == ["", ""]
    assert overlaid.genes != model.genes


def test_gene_rule_evidence_cache_writes_structured_display_only_records(tmp_path: Path) -> None:
    output_path = tmp_path / "gene_rule_evidence.json"
    summary_path = tmp_path / "summary.json"

    summary = build_gene_rule_evidence_cache(
        ("PDI1", "ERO1", "OCH1"),
        output_path=output_path,
        summary_path=summary_path,
        enable_online=False,
    )
    loaded = load_gene_rule_evidence_cache(output_path)
    payload = json.loads(output_path.read_text(encoding="utf-8"))

    assert summary["total_records"] == 3
    assert payload["schema_version"] == 1
    assert set(loaded) == {"PDI1", "ERO1", "OCH1"}
    assert loaded["PDI1"].candidate_locus_tag
    assert loaded["PDI1"].rule_status == "display_only_requires_multi_source_confirmation"
    assert loaded["OCH1"].confidence == "unresolved"


def test_gene_rule_evidence_cache_enriches_exact_locus_without_executable_complex(
    tmp_path: Path,
    monkeypatch,
) -> None:
    output_path = tmp_path / "gene_rule_evidence.json"
    summary_path = tmp_path / "summary.json"
    report_path = tmp_path / DEFAULT_GENE_RULE_EVIDENCE_REPORT.name

    monkeypatch.setattr(
        gene_rule_overlay,
        "_fetch_source_rows_for_records",
        lambda *args, **kwargs: {
            "PAS_chr1-1_0160": {
                "uniprot": {
                    "Entry": "C4QWA2",
                    "Gene Names": "PAS_chr1-1_0160",
                    "Protein names": "protein disulfide-isomerase (EC 5.3.4.1)",
                    "KEGG": "ppa:PAS_chr1-1_0160;",
                    "GeneID": "8196664;",
                },
                "kegg_description": "Protein disulfide isomerase",
            },
            "PAS_chr1-1_0011": {
                "uniprot": {
                    "Entry": "C4QVU1",
                    "Gene Names": "PAS_chr1-1_0011",
                    "Protein names": "thiol oxidase (EC 1.8.3.2)",
                    "KEGG": "ppa:PAS_chr1-1_0011;",
                    "GeneID": "8197528;",
                },
                "kegg_description": "Thiol oxidase required for oxidative protein folding",
            },
        },
    )

    summary = build_gene_rule_evidence_cache(
        ("PDI1", "ERO1", "ERV2"),
        output_path=output_path,
        summary_path=summary_path,
        report_path=report_path,
        model=TinyModel(),
    )
    loaded = load_gene_rule_evidence_cache(output_path)

    assert summary["high_confidence_count"] == 2
    assert summary["executable_overlay_entry_count"] == 0
    assert loaded["PDI1"].confidence == HIGH_CONFIDENCE
    assert loaded["ERO1"].confidence == HIGH_CONFIDENCE
    assert loaded["ERV2"].confidence == "unresolved"
    assert loaded["PDI1"].external_ids["uniprot"] == "C4QWA2"
    assert "UniProt GS115 proteome exact locus" in loaded["PDI1"].evidence_sources
    assert "PDI1/ERO1/ERV2 共享复合体反应" in report_path.read_text(encoding="utf-8")


def _high(common_name: str, locus_tag: str) -> GeneRuleEvidence:
    return GeneRuleEvidence(
        common_name=common_name,
        candidate_locus_tag=locus_tag,
        evidence_sources=("unit_test_source_a", "unit_test_source_b"),
        confidence=HIGH_CONFIDENCE,
        target_reaction_ids=("sec_PDI1_ERV2_Ero1p_complex_formation",),
    )
