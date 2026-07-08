from __future__ import annotations

import ast
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_hlf_opn_gpr_overlay_review_facade_uses_cache_without_applying_overlay(tmp_path, monkeypatch) -> None:
    from pcsec_pichia.core.paths import ProjectPaths
    from pcsec_pichia.services.gene_rule_overlay import GeneRuleEvidence, HIGH_CONFIDENCE

    from app.services import pichia_gene_catalog_service as service

    calls = {"candidates": 0, "evidence": 0, "reactions": 0}

    def fake_candidates(*, paths=None, **kwargs) -> list[dict[str, object]]:
        calls["candidates"] += 1
        return [
            {
                "target_context": "hLF",
                "source_common_name": "PDI1",
                "gene_id": "PAS_PDI1_HOMOLOGY",
                "candidate_role": "disulfide_bond_folding",
                "recommended_intervention": "OE",
                "model_operable": False,
            }
        ]

    def fake_evidence(paths=None) -> dict[str, GeneRuleEvidence]:
        calls["evidence"] += 1
        return {
            "PDI1": GeneRuleEvidence(
                common_name="PDI1",
                candidate_locus_tag="PAS_PDI1_EXTERNAL",
                external_ids={"uniprot": "U-PDI1"},
                evidence_sources=("UniProt GS115 proteome exact locus", "KEGG ppa exact locus"),
                confidence=HIGH_CONFIDENCE,
                target_reaction_ids=("sec_PDI1_ERV2_Ero1p_complex_formation",),
                rule_status="high_confidence_locus_candidate",
                recommended_action="eligible_for_overlay_if_all_complex_subunits_are_confirmed",
            )
        }

    def fake_reactions(paths=None) -> set[str]:
        calls["reactions"] += 1
        return {"sec_PDI1_ERV2_Ero1p_complex_formation"}

    monkeypatch.setattr(service, "load_hlf_opn_candidate_genes", fake_candidates)
    monkeypatch.setattr(service, "_load_offline_gene_rule_evidence", fake_evidence)
    monkeypatch.setattr(service, "_load_model_reaction_ids", fake_reactions)
    paths = ProjectPaths(repo_root=tmp_path)

    first_rows = service.load_hlf_opn_gpr_overlay_review(target_context="hLF", paths=paths)
    second_rows = service.load_hlf_opn_gpr_overlay_review(target_context="hLF", paths=paths)
    summary = service.hlf_opn_gpr_overlay_review_summary(paths=paths)

    assert calls == {"candidates": 1, "evidence": 1, "reactions": 1}
    assert first_rows == second_rows
    assert first_rows[0]["review_status"] == "candidate_gpr_overlay_review"
    assert first_rows[0]["gene_id"] == "PAS_PDI1_EXTERNAL"
    assert summary["candidate_gpr_overlay_review_count"] == 1
    assert service.hlf_opn_gpr_overlay_review_cache_path(paths).exists()


def test_hlf_opn_gpr_overlay_review_facade_does_not_import_runtime_overlay_or_online_clients() -> None:
    service_path = REPO_ROOT / "app" / "services" / "pichia_gene_catalog_service.py"
    module_ast = ast.parse(service_path.read_text(encoding="utf-8"))
    target_functions = {
        "load_hlf_opn_gpr_overlay_review",
        "hlf_opn_gpr_overlay_review_summary",
        "hlf_opn_gpr_overlay_review_cache_path",
        "_load_offline_gene_rule_evidence",
        "_load_model_reaction_ids",
    }
    function_nodes = [
        node
        for node in module_ast.body
        if isinstance(node, ast.FunctionDef) and node.name in target_functions
    ]
    imports: set[str] = set()
    called_names: set[str] = set()
    for function_node in function_nodes:
        for node in ast.walk(function_node):
            if isinstance(node, ast.ImportFrom) and node.module:
                imports.add(node.module)
            elif isinstance(node, ast.Import):
                imports.update(alias.name for alias in node.names)
            elif isinstance(node, ast.Call):
                if isinstance(node.func, ast.Name):
                    called_names.add(node.func.id)
                elif isinstance(node.func, ast.Attribute):
                    called_names.add(node.func.attr)

    assert {node.name for node in function_nodes} == target_functions
    assert "pcsec_pichia.homology.blast_runner" not in imports
    assert "pcsec_pichia.homology.external_fetch" not in imports
    assert "urllib.request" not in imports
    assert not {
        "run_blastp",
        "make_blast_db",
        "fetch_external_name_references",
        "build_gene_rule_evidence_cache",
        "build_gpr_overlay",
        "apply_gpr_overlay_for_analysis",
        "urlopen",
    } & called_names
