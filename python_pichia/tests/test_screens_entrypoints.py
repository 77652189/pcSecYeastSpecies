from __future__ import annotations

import ast
import os
from functools import lru_cache
from pathlib import Path
from types import SimpleNamespace

import pytest

from pcsec_pichia.loading import PcSecPichiaInputs, load_pcsec_pichia_inputs
from pcsec_pichia.screens import (
    ScreenResult,
    build_all_gene_capability_catalog,
    build_gene_capability_profile,
    build_gene_perturbation_map,
    explain_only_gene_overexpression_rows,
    plan_gene_knockout,
    plan_gene_overexpression,
    run_knockout_screen,
    run_overexpression_screen,
    run_reaction_knockout_screen,
    summarize_screen_result,
)
from pcsec_pichia.screens.candidate_resolution import reactions_for_gene
from pcsec_pichia.services.gene_evidence import (
    GeneExternalEvidence,
    GenePhenotypeEvidence,
    load_gene_evidence_cache,
    phenotype_evidence_for_candidate,
    recommendation_tier_for_candidate,
)
from pcsec_pichia.targets import TargetSpec, load_builtin_targets


REPO_ROOT = Path(__file__).resolve().parents[2]
KO_GENES = ["AT250_GQ_6803479", "AT250_GQ_6803809"]
OE_REACTIONS = ["sec_BIP_NEFS_complex_formation", "sec_Kar2p_complex_formation"]
KO_REACTIONS = ["sec_Och1p_complex_formation"]
slow_screen = pytest.mark.skipif(
    os.environ.get("PCSEC_RUN_SLOW_SCREEN_TESTS") != "1",
    reason="slow pcSec screen solve; set PCSEC_RUN_SLOW_SCREEN_TESTS=1 to run",
)
REQUIRED_ROW_FIELDS = {
    "candidate_id",
    "intervention_type",
    "objective_value",
    "baseline_objective_value",
    "delta_objective",
    "canonical_gene_id",
    "input_gene_id",
    "resolved_reaction_id",
    "effect_label",
    "solver_status_label",
    "failure_reason",
    "secretory_process",
    "mapping_level",
    "mapping_confidence",
    "mapping_interpretation",
    "complex_id",
    "affected_reactions",
    "inactive_reactions",
    "inactive_reaction_count",
    "gpr_rules",
    "gpr_role",
    "capacity_effect",
    "simulation_basis",
    "ko_support_status",
    "oe_support_status",
    "support_reason",
    "missing_information",
    "complex_subunit_ids",
    "complex_subunit_stoichiometry",
}


def test_gene_perturbation_map_explains_gene_reaction_process_and_confidence() -> None:
    class TinyModel:
        rxns = [
            "sec_Kar2p_complex_formation",
            "sec_PDI1_ERV2_Ero1p_complex_formation",
            "sec_Kar2p_activity",
            "BIOMASS",
        ]
        rules = ["x(1)", "x(1)", "x(2)", "x(3)"]
        gr_rules = ["G1", "G1", "G2", "G3"]
        gene_index = {"G1": 0, "G2": 1, "G3": 2}
        reaction_index = {
            "sec_Kar2p_complex_formation": 0,
            "sec_PDI1_ERV2_Ero1p_complex_formation": 1,
            "sec_Kar2p_activity": 2,
            "BIOMASS": 3,
        }

    result = build_gene_perturbation_map(
        TinyModel(),
        ("G1", "G2", "G3", "NO_SUCH_GENE"),
        complex_subunits={
            "sec_Kar2p_complex": [
                {"subunit_id": "Kar2p", "stoichiometry": 1.0},
            ],
        },
    )
    rows = result.rows

    g1_rows = [row for row in rows if row["gene_id"] == "G1"]
    assert len(g1_rows) == 2
    assert all(row["reaction_count"] == 2 for row in g1_rows)
    assert any(row["secretory_process"] == "ER 折叠 / 分子伴侣" for row in g1_rows)

    complex_row = next(row for row in rows if row["reaction_id"] == "sec_Kar2p_complex_formation")
    assert complex_row["mapping_level"] == "complex_subunit"
    assert complex_row["mapping_confidence"] == "medium"
    assert complex_row["complex_subunit_ids"] == ["Kar2p"]
    assert complex_row["complex_subunit_stoichiometry"] == [1.0]

    direct_row = next(row for row in rows if row["reaction_id"] == "sec_Kar2p_activity")
    assert direct_row["mapping_level"] == "direct_gpr"
    assert direct_row["mapping_confidence"] == "high"

    metabolic_row = next(row for row in rows if row["reaction_id"] == "BIOMASS")
    assert metabolic_row["mapping_level"] == "metabolic_or_other"
    assert metabolic_row["mapping_confidence"] == "low"

    unresolved_row = next(row for row in rows if row["gene_id"] == "NO_SUCH_GENE")
    assert unresolved_row["mapping_level"] == "unresolved"
    assert unresolved_row["mapping_confidence"] == "unresolved"
    assert unresolved_row["resolved"] is False


def test_formal_screen_row_normalization_adds_mapping_explanation_fields() -> None:
    import pcsec_pichia.screens as screens_module

    row = screens_module._normalize_screen_row(
        {
            "reaction": "sec_Kar2p_complex_formation",
            "success": True,
            "status": "0",
            "objective_value": 0.002,
            "delta_vs_baseline": 0.0001,
        },
        target_id="OPN_ALPHA_FULL_PROJECT",
        screen_type="overexpression",
        intervention_type="OE_reaction",
        baseline_objective_value=0.0019,
        complex_subunits={
            "sec_Kar2p_complex": [
                {"subunit_id": "Kar2p", "stoichiometry": 1.0},
                {"subunit_id": "Sil1p", "stoichiometry": 1.0},
            ],
        },
    )

    assert row["mapping_level"] == "complex_subunit"
    assert row["mapping_confidence"] == "medium"
    assert row["mapping_interpretation"]
    assert row["complex_id"] == "sec_Kar2p_complex"
    assert row["complex_subunit_ids"] == ["Kar2p", "Sil1p"]
    assert row["secretory_process"] == "ER 折叠 / 分子伴侣"


def test_gene_intervention_plans_respect_gpr_and_or_rules() -> None:
    class TinyModel:
        rxns = ["R_AND", "R_OR", "R_MIXED", "R_SINGLE"]
        rules = ["x(1) & x(2)", "x(1) | x(2)", "(x(1) & x(2)) | x(3)", "x(4)"]
        gr_rules = ["G1 and G2", "G1 or G2", "(G1 and G2) or G3", "G4"]
        gene_index = {"G1": 0, "G2": 1, "G3": 2, "G4": 3}
        reaction_index = {"R_AND": 0, "R_OR": 1, "R_MIXED": 2, "R_SINGLE": 3}

    ko = plan_gene_knockout(TinyModel(), "G1")
    assert ko.resolved is True
    assert "R_AND" in ko.inactive_reactions
    assert "R_OR" not in ko.inactive_reactions
    assert "R_MIXED" not in ko.inactive_reactions
    assert ko.simulation_basis == "gpr_gene_deletion"
    assert ko.capacity_effect == "disables_reactions"

    oe_direct = plan_gene_overexpression(TinyModel(), "G4")
    assert oe_direct.executable_reactions == ("R_SINGLE",)
    assert oe_direct.explain_only_reactions == ()
    assert oe_direct.gpr_role == "single_gene"
    assert oe_direct.capacity_effect == "reaction_capacity_proxy"

    unresolved = plan_gene_knockout(TinyModel(), "NO_SUCH_GENE")
    assert unresolved.resolved is False
    assert unresolved.gpr_role == "unresolved"

    TinyModel.genes = ["G1", "G2", "G3", "G4", "G5"]
    TinyModel.gene_index = {"G1": 0, "G2": 1, "G3": 2, "G4": 3, "G5": 4}
    ko_no_gpr = plan_gene_knockout(TinyModel(), "G5")
    oe_no_gpr = plan_gene_overexpression(TinyModel(), "G5")
    assert ko_no_gpr.resolved is True
    assert ko_no_gpr.ko_support_status == "ko_no_gpr_effect"
    assert oe_no_gpr.resolved is True
    assert oe_no_gpr.oe_support_status == "oe_no_gpr_effect"


def test_gene_knockout_solver_uses_planned_inactive_reactions(monkeypatch) -> None:
    import pcsec_pichia.screens as screens_module

    captured: dict[str, object] = {}

    class TinyModel:
        def with_bounds(self, changes):
            captured["changes"] = changes
            return self

    def fake_solve(model, objective, **kwargs):
        captured["objective"] = objective
        captured["kwargs"] = kwargs
        return SimpleNamespace(status="0", success=True, objective_value=0.8), {"eq_total": 1}

    monkeypatch.setattr(screens_module, "solve_pcsec_maximize", fake_solve)
    plan = screens_module.GeneInterventionPlan(
        gene_id="G1",
        intervention_type="KO",
        resolved=True,
        affected_reactions=("R_AND", "R_OR"),
        inactive_reactions=("R_AND",),
        executable_reactions=("R_AND",),
        explain_only_reactions=(),
        gpr_rules=(),
        gpr_role="mixed",
        capacity_effect="disables_reactions",
        simulation_basis="gpr_gene_deletion",
        mapping_confidence="high",
        warnings=(),
        ko_support_status="ko_runnable_gpr_gene_deletion",
    )
    row = screens_module._solve_gene_knockout_plan(
        {
            "fixed_model": TinyModel(),
            "baseline": SimpleNamespace(objective_value=1.0),
            "secretory": object(),
            "combined": object(),
            "exchange_reaction_id": "r_target_exchange",
        },
        plan,
        metabolic=object(),
        growth_rate=0.10,
        write_ribosome_translation_constraint=False,
        write_misfolding_constraints=False,
    )

    assert captured["changes"] == {"R_AND": (0.0, 0.0)}
    assert captured["objective"] == "r_target_exchange"
    assert row["gene"] == "G1"
    assert row["inactive_reactions"] == ["R_AND"]
    assert row["inactive_reaction_count"] == 1
    assert row["delta_vs_baseline"] == pytest.approx(-0.2)


def test_gene_capability_profile_covers_gpr_statuses_for_all_model_genes() -> None:
    class TinyModel:
        rxns = ["R_AND", "R_OR", "R_SINGLE"]
        rules = ["x(1) & x(2)", "x(1) | x(2)", "x(4)"]
        gr_rules = ["G1 and G2", "G1 or G2", "G4"]
        genes = ["G1", "G2", "G3", "G4"]
        gene_index = {"G1": 0, "G2": 1, "G3": 2, "G4": 3}
        reaction_index = {"R_AND": 0, "R_OR": 1, "R_SINGLE": 2}

    g1 = build_gene_capability_profile(TinyModel(), "G1").to_dict()
    g3 = build_gene_capability_profile(TinyModel(), "G3").to_dict()
    missing = build_gene_capability_profile(TinyModel(), "NO_SUCH_GENE").to_dict()
    catalog = [profile.to_dict() for profile in build_all_gene_capability_catalog(TinyModel())]

    assert g1["ko_support_status"] == "ko_runnable_gpr_gene_deletion"
    assert g1["oe_support_status"] == "oe_runnable_reaction_proxy"
    assert g1["inactive_reactions_if_ko"] == ["R_AND"]
    assert g3["ko_support_status"] == "ko_no_gpr_effect"
    assert g3["oe_support_status"] == "oe_no_gpr_effect"
    assert "model_gpr_rule" in g3["missing_information"]
    assert missing["ko_support_status"] == "unresolved_gene"
    assert missing["oe_support_status"] == "unresolved_gene"
    assert {row["gene_id"] for row in catalog} == {"G1", "G2", "G3", "G4"}
    assert all(row["ko_support_status"] and row["oe_support_status"] for row in catalog)


def test_gene_capability_profile_keeps_ko_oe_phenotype_tiers_separate() -> None:
    class TinyModel:
        rxns = ["R1"]
        rules = ["x(1)"]
        gr_rules = ["ESSENTIAL_TEST"]
        genes = ["ESSENTIAL_TEST"]
        gene_index = {"ESSENTIAL_TEST": 0}
        reaction_index = {"R1": 0}

    profile = build_gene_capability_profile(
        TinyModel(),
        "ESSENTIAL_TEST",
        target_protein_context="OPN_ALPHA_FULL_PROJECT",
    ).to_dict()

    assert profile["recommendation_tier"]["KO"] == "not_recommended_growth_risk"
    assert profile["recommendation_tier"]["OE"] == "model_executable"
    assert profile["phenotype_evidence"]["KO"]["essentiality_status"] == "essential"
    assert "OE" not in profile["phenotype_evidence"]


def test_ko_essentiality_remains_highest_priority_when_multiple_phenotype_records() -> None:
    evidence_by_gene = {
        "RISKY": (
            GenePhenotypeEvidence(
                intervention_type="KO",
                essentiality_status="nonessential",
                secretion_screen_effect="screen_positive",
                evidence_source="same_target_same_host_same_intervention",
                evidence_confidence="high",
                target_protein_context=("OPN_ALPHA_FULL_PROJECT",),
            ),
            GenePhenotypeEvidence(
                intervention_type="KO",
                essentiality_status="essential",
                secretion_screen_effect="growth_risk",
                evidence_source="public_screen_positive",
                evidence_confidence="medium",
            ),
        )
    }
    phenotype = phenotype_evidence_for_candidate(
        "RISKY",
        "KO",
        target_protein_context="OPN_ALPHA_FULL_PROJECT",
        evidence_by_gene=evidence_by_gene,
    )
    tier, _reason, matched = recommendation_tier_for_candidate(
        gene_id="RISKY",
        intervention_type="KO",
        target_protein_context="OPN_ALPHA_FULL_PROJECT",
        model_gpr_executable=True,
        phenotype_evidence=phenotype,
    )

    assert matched is phenotype
    assert phenotype is not None
    assert phenotype.essentiality_status == "essential"
    assert tier == "not_recommended_growth_risk"


def test_phenotype_evidence_can_be_explicitly_disabled_for_isolated_callers() -> None:
    phenotype = phenotype_evidence_for_candidate(
        "ESSENTIAL_TEST",
        "KO",
        target_protein_context="OPN_ALPHA_FULL_PROJECT",
        evidence_by_gene={},
    )

    assert phenotype is None


def test_recommendation_tier_ignores_explicit_phenotype_with_wrong_intervention() -> None:
    wrong_intervention_evidence = GenePhenotypeEvidence(
        intervention_type="KO",
        essentiality_status="essential",
        secretion_screen_effect="growth_risk",
        evidence_source="internal_curated",
        evidence_confidence="high",
        target_protein_context=("OPN_ALPHA_FULL_PROJECT",),
    )

    tier, _reason, matched = recommendation_tier_for_candidate(
        gene_id="ESSENTIAL_TEST",
        intervention_type="OE_gene_proxy",
        target_protein_context="OPN_ALPHA_FULL_PROJECT",
        oe_reaction_proxy=True,
        phenotype_evidence=wrong_intervention_evidence,
    )

    assert matched is None
    assert tier == "model_executable"


def test_gene_capability_profile_resolves_alias_from_offline_evidence() -> None:
    class TinyModel:
        rxns = ["R1"]
        rules = ["x(1)"]
        gr_rules = ["G1"]
        genes = ["G1"]
        gene_index = {"G1": 0}
        reaction_index = {"R1": 0}

    evidence = {
        "G1": GeneExternalEvidence(
            gene_id="G1",
            canonical_gene_id="G1",
            aliases=("ALIAS1",),
            evidence_sources=("offline_cache",),
        )
    }

    profile = build_gene_capability_profile(TinyModel(), "ALIAS1", evidence_by_gene=evidence).to_dict()

    assert profile["gene_id"] == "ALIAS1"
    assert profile["canonical_gene_id"] == "G1"
    assert profile["aliases"] == ["ALIAS1"]
    assert profile["resolved"] is True
    assert profile["ko_support_status"] == "ko_runnable_gpr_gene_deletion"
    assert profile["oe_support_status"] == "oe_runnable_reaction_proxy"


def test_gene_evidence_cache_ignores_corrupt_json(tmp_path: Path) -> None:
    cache_path = tmp_path / "gene_evidence.json"
    cache_path.write_text("{not valid json", encoding="utf-8")

    assert load_gene_evidence_cache(cache_path) == {}


def test_gene_evidence_cache_ignores_corrupt_utf8_bytes(tmp_path: Path) -> None:
    cache_path = tmp_path / "gene_evidence.json"
    cache_path.write_bytes(b"\xff\xfe\x00bad")

    assert load_gene_evidence_cache(cache_path) == {}


def test_gene_evidence_mapping_preserves_wet_lab_annotation_fields() -> None:
    from pcsec_pichia.services.gene_evidence import gene_external_evidence_from_mapping

    record = gene_external_evidence_from_mapping(
        {
            "gene_id": "PAS_chr1-4_0141",
            "canonical_gene_id": "PAS_chr1-4_0141",
            "standard_gene_symbol": "THI20",
            "display_name": "Hydroxymethylpyrimidine phosphate kinase",
            "protein_name": "Protein with similarity to hydroxymethylpyrimidine phosphate kinases",
            "function_annotation": "thiamine biosynthetic process",
            "external_ids": {"uniprot": "C4QXJ6", "ncbi_gene": "8197125", "kegg": "ppa:PAS_chr1-4_0141"},
            "ec_numbers": ["2.7.4.7"],
            "go_terms": ["thiamine biosynthetic process [GO:0009228]"],
            "wet_lab_readiness": "database_supported_experiment_candidate",
            "evidence_sources": ["UniProt", "NCBI Gene", "KEGG"],
            "evidence_confidence": "high_exact_locus_tag",
        }
    )

    assert record.model_gene_id == "PAS_chr1-4_0141"
    assert record.standard_gene_symbol == "THI20"
    assert record.display_name == "Hydroxymethylpyrimidine phosphate kinase"
    assert record.external_ids["uniprot"] == "C4QXJ6"
    assert record.wet_lab_readiness == "database_supported_experiment_candidate"
    assert record.to_dict()["go_terms"] == ["thiamine biosynthetic process [GO:0009228]"]


def test_gene_evidence_builder_records_exact_and_model_only_genes(tmp_path: Path, monkeypatch) -> None:
    from pcsec_pichia.services import gene_evidence

    monkeypatch.setattr(
        gene_evidence,
        "_fetch_uniprot_proteome_by_gene",
        lambda *args, **kwargs: {
            "PAS_chr1-4_0141": {
                "Entry": "C4QXJ6",
                "Gene Names": "PAS_chr1-4_0141",
                "Gene Names (primary)": "",
                "Protein names": "Protein with similarity to hydroxymethylpyrimidine phosphate kinases",
                "GeneID": "8197125;",
                "KEGG": "ppa:PAS_chr1-4_0141;",
                "RefSeq": "XP_002490250.1;",
                "Gene Ontology (biological process)": "thiamine biosynthetic process [GO:0009228]",
                "EC number": "",
            }
        },
    )
    monkeypatch.setattr(
        gene_evidence,
        "_fetch_kegg_gene_descriptions",
        lambda *args, **kwargs: {
            "PAS_chr1-4_0141": "Protein with similarity to hydroxymethylpyrimidine phosphate kinases"
        },
    )

    output_path = tmp_path / "gene_evidence.json"
    summary_path = tmp_path / "summary.json"
    summary = gene_evidence.build_gene_evidence_cache(
        ["PAS_chr1-4_0141", "AT250_GQ_6803479"],
        output_path=output_path,
        summary_path=summary_path,
    )
    loaded = gene_evidence.load_gene_evidence_cache(output_path)

    assert summary["total_genes"] == 2
    assert summary["database_supported_count"] == 1
    assert summary["model_only_count"] == 1
    assert loaded["PAS_chr1-4_0141"].external_ids["uniprot"] == "C4QXJ6"
    assert loaded["PAS_chr1-4_0141"].wet_lab_readiness == "database_supported_experiment_candidate"
    assert loaded["AT250_GQ_6803479"].wet_lab_readiness == "model_only_not_experiment_ready"


def test_gene_knockout_does_not_disable_text_only_gr_rule() -> None:
    class TinyModel:
        rxns = ["R_TEXT_ONLY"]
        rules = ["[]"]
        gr_rules = ["G1"]
        genes = ["G1"]
        gene_index = {"G1": 0}
        reaction_index = {"R_TEXT_ONLY": 0}

    plan = plan_gene_knockout(TinyModel(), "G1")

    assert plan.resolved is True
    assert plan.affected_reactions == ("R_TEXT_ONLY",)
    assert plan.inactive_reactions == ()
    assert plan.ko_support_status == "ko_no_reaction_disabled"
    assert "model_rule_token_mapping" in plan.missing_information


def test_reactions_for_gene_handles_missing_parallel_rule_lists() -> None:
    class TinyModel:
        rxns = ["R_NUMERIC", "R_TEXT_ONLY", "R_NO_RULE"]
        rules = ["x(1)"]
        gr_rules = ["", "G1"]
        genes = ["G1"]
        gene_index = {"G1": 0}
        reaction_index = {"R_NUMERIC": 0, "R_TEXT_ONLY": 1, "R_NO_RULE": 2}

    assert reactions_for_gene(TinyModel(), "G1") == ["R_NUMERIC", "R_TEXT_ONLY"]


def test_gene_overexpression_complex_subunit_is_explain_only() -> None:
    class TinyModel:
        rxns = ["sec_Kar2p_complex_formation"]
        rules = ["x(1)"]
        gr_rules = ["Kar2p"]
        gene_index = {"Kar2p": 0}
        reaction_index = {"sec_Kar2p_complex_formation": 0}

    plan = plan_gene_overexpression(
        TinyModel(),
        "Kar2p",
        complex_subunits={"sec_Kar2p_complex": [{"subunit_id": "Kar2p", "stoichiometry": 1.0}]},
    )
    rows = explain_only_gene_overexpression_rows(
        "OPN_ALPHA_FULL_PROJECT",
        (plan,),
        baseline_objective_value=0.1,
        complex_subunits={"sec_Kar2p_complex": [{"subunit_id": "Kar2p", "stoichiometry": 1.0}]},
    )

    assert plan.executable_reactions == ()
    assert plan.explain_only_reactions == ("sec_Kar2p_complex_formation",)
    assert plan.capacity_effect == "complex_subunit_limited"
    assert rows[0]["status"] == "not_run_complex_subunit_limited"
    assert rows[0]["effect_label"] == "未运行"
    assert rows[0]["simulation_basis"] == "explain_only"
    assert rows[0]["capacity_effect"] == "complex_subunit_limited"


def test_legacy_oe_gene_resolver_is_compat_export_only() -> None:
    import pcsec_pichia.screens as screens

    preview_source = (REPO_ROOT / "app" / "services" / "pichia_screen_preview_service.py").read_text(encoding="utf-8")
    assert hasattr(screens, "resolve_oe_gene_reactions")
    assert "resolve_oe_gene_reactions" not in preview_source


def test_screen_solve_tests_are_slow_gated() -> None:
    module_ast = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    solve_calls = {"run_knockout_screen", "run_overexpression_screen", "run_reaction_knockout_screen"}
    ungated: list[str] = []
    for node in module_ast.body:
        if not isinstance(node, ast.FunctionDef) or not node.name.startswith("test_"):
            continue
        calls_screen_solve = any(
            isinstance(child, ast.Call)
            and isinstance(child.func, ast.Name)
            and child.func.id in solve_calls
            for child in ast.walk(node)
        )
        has_slow_marker = any(
            isinstance(decorator, ast.Name) and decorator.id == "slow_screen"
            for decorator in node.decorator_list
        )
        if calls_screen_solve and not has_slow_marker:
            ungated.append(node.name)

    assert ungated == []


@lru_cache(maxsize=1)
def _inputs() -> PcSecPichiaInputs:
    return load_pcsec_pichia_inputs(REPO_ROOT)


@lru_cache(maxsize=1)
def _builtin_targets() -> dict[str, TargetSpec]:
    return {target.target_id: target for target in load_builtin_targets(REPO_ROOT)}


def _assert_common_screen_result(result: ScreenResult, target_id: str, screen_type: str) -> None:
    summary = summarize_screen_result(result)

    assert isinstance(result, ScreenResult)
    assert result.target_id == target_id
    assert result.screen_type == screen_type
    assert result.success is True
    assert result.candidate_count == 2
    assert len(result.rows) == 2
    assert result.constraint_counts["eq_total"] > 0
    assert result.constraint_counts["ub_total"] == 1
    assert result.baseline_objective_value is not None
    assert result.result_status == "draft"
    assert result.matlab_alignment_status == "pending"
    assert summary["candidate_count"] == 2

    for row in result.rows:
        assert REQUIRED_ROW_FIELDS.issubset(row)
        assert row["objective_value"] is not None
        assert row["baseline_objective_value"] == pytest.approx(result.baseline_objective_value)


@slow_screen
@pytest.mark.parametrize("target_id", ("OPN_ALPHA_FULL_PROJECT", "hLF"))
def test_builtin_targets_run_knockout_screen_smoke(target_id: str) -> None:
    inputs = _inputs()
    target = _builtin_targets()[target_id]

    result = run_knockout_screen(
        inputs.prepared_model,
        target,
        inputs.amino_acids,
        inputs.metabolic,
        inputs.secretory,
        inputs.combined,
        genes=list(KO_GENES),
    )

    _assert_common_screen_result(result, target_id, "knockout")
    for row in result.rows:
        assert row["intervention_type"] == "KO"
        assert row["gene_id"] in KO_GENES
        assert row["input_gene_id"] == row["gene_id"]
        assert row["resolved_reaction_id"] == row["reaction_id"]
        assert row["effect_label"] in {"提升分泌", "降低分泌", "无明显变化", "求解失败", "未解析"}
        assert row["secretory_process"]
        assert row["reaction_id"]


@slow_screen
@pytest.mark.parametrize("target_id", ("OPN_ALPHA_FULL_PROJECT", "hLF"))
def test_builtin_targets_run_overexpression_screen_smoke(target_id: str) -> None:
    inputs = _inputs()
    target = _builtin_targets()[target_id]

    result = run_overexpression_screen(
        inputs.prepared_model,
        target,
        inputs.amino_acids,
        inputs.metabolic,
        inputs.secretory,
        inputs.combined,
        reactions=list(OE_REACTIONS),
    )

    _assert_common_screen_result(result, target_id, "overexpression")
    for row in result.rows:
        assert row["intervention_type"] == "OE_reaction"
        assert row["gene_id"] is None
        assert row["input_gene_id"] is None
        assert row["resolved_reaction_id"] == row["reaction_id"]
        assert row["effect_label"] in {"提升分泌", "降低分泌", "无明显变化", "求解失败", "未解析"}
        assert row["secretory_process"]
        assert row["reaction_id"] in OE_REACTIONS
        assert row["complex_subunit_ids"]
        assert row["complex_subunit_stoichiometry"]


@slow_screen
def test_builtin_opn_runs_reaction_knockout_screen_smoke() -> None:
    inputs = _inputs()
    target = _builtin_targets()["OPN_ALPHA_FULL_PROJECT"]

    result = run_reaction_knockout_screen(
        inputs.prepared_model,
        target,
        inputs.amino_acids,
        inputs.metabolic,
        inputs.secretory,
        inputs.combined,
        reactions=list(KO_REACTIONS),
        write_ribosome_translation_constraint=True,
        write_misfolding_constraints=True,
    )

    summary = summarize_screen_result(result)

    assert isinstance(result, ScreenResult)
    assert result.target_id == "OPN_ALPHA_FULL_PROJECT"
    assert result.screen_type == "knockout"
    assert result.candidate_count == 1
    assert len(result.rows) == 1
    assert summary["candidate_count"] == 1
    assert result.constraint_counts["ribosome_translation"] == 1
    assert result.constraint_counts["misfolding"] == 1418
    row = result.rows[0]
    assert REQUIRED_ROW_FIELDS.issubset(row)
    assert row["intervention_type"] == "KO_reaction"
    assert row["gene_id"] is None
    assert row["reaction_id"] == KO_REACTIONS[0]
    assert row["resolved_reaction_id"] == KO_REACTIONS[0]
    assert row["baseline_objective_value"] == pytest.approx(result.baseline_objective_value)
    assert row["complex_subunit_ids"]
    assert row["status"] == "2"
    assert row["success"] is False
    assert row["objective_value"] is None
    assert row["delta_objective"] is None
    assert row["effect_label"] == "约束不可行"
    assert row["solver_status_label"] == "约束不可行"
    assert row["failure_reason"] == "约束不可行"
