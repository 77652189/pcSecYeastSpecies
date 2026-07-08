from __future__ import annotations

import ast
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_hlf_opn_candidate_gene_facade_uses_persistent_cache(tmp_path, monkeypatch) -> None:
    from pcsec_pichia.core.paths import ProjectPaths
    from pcsec_pichia.services.gene_catalog import SecretionGeneEntry
    from pcsec_pichia.services.homology_evidence import GeneHomologyEvidence

    from app.services import pichia_gene_catalog_service as service
    import pcsec_pichia.services.gene_catalog as gene_catalog
    import pcsec_pichia.services.homology_evidence as homology_evidence

    calls = {"full": 0, "standard": 0}

    def fake_full_model_catalog(*, force_refresh: bool = False, paths=None) -> list[dict[str, object]]:
        calls["full"] += 1
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
            }
        ]

    def fake_standard_names(*, force_refresh: bool = False, paths=None) -> list[dict[str, object]]:
        calls["standard"] += 1
        return [
            {
                "gene_id": "PAS_KAR2",
                "display_name": "KAR2",
                "standard_symbol": "KAR2",
                "protein_name": "ER chaperone",
                "external_ids": {"uniprot": "U-KAR2"},
                "annotation_sources": ["UniProt"],
                "annotation_confidence": "high_exact_locus_tag",
                "model_operable": True,
                "gpr_status": "ko_and_oe_model_executable",
            }
        ]

    monkeypatch.setattr(service, "load_pichia_full_model_gene_catalog", fake_full_model_catalog)
    monkeypatch.setattr(service, "load_pichia_gene_id_standardization", fake_standard_names)
    monkeypatch.setattr(
        gene_catalog,
        "SECRETION_GENE_CATALOG",
        (
            SecretionGeneEntry(
                category="ER 折叠与分子伴侣",
                common_name="KAR2 / BiP",
                description="ER chaperone",
                intervention="OE",
                oe_reaction_id="sec_Kar2p_complex_formation",
                evidence="已报道 Kar2 过表达可提升毕赤酵母外源蛋白分泌",
            ),
        ),
    )
    monkeypatch.setattr(
        homology_evidence,
        "load_homology_evidence_cache",
        lambda *args, **kwargs: {
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
            )
        },
    )
    paths = ProjectPaths(repo_root=tmp_path)

    first_rows = service.load_hlf_opn_candidate_genes(target_context="hLF", paths=paths)
    second_rows = service.load_hlf_opn_candidate_genes(target_context="hLF", paths=paths)
    summary = service.hlf_opn_candidate_gene_summary(paths=paths)
    executable = service.hlf_opn_executable_candidate_inputs(target_context="hLF", paths=paths)

    assert calls == {"full": 1, "standard": 1}
    assert first_rows == second_rows
    assert first_rows[0]["gene_id"] == "PAS_KAR2"
    assert first_rows[0]["operability_status"] == "model_oe_proxy_executable"
    assert summary["total_candidates"] == 1
    assert executable["oe_gene_ids"] == ["PAS_KAR2"]
    assert service.hlf_opn_candidate_gene_cache_path(paths).exists()
    assert "local_runs" in str(service.hlf_opn_candidate_gene_cache_path(paths))


def test_hlf_opn_candidate_gene_facade_does_not_import_blast_or_online_clients() -> None:
    service_path = REPO_ROOT / "app" / "services" / "pichia_gene_catalog_service.py"
    module_ast = ast.parse(service_path.read_text(encoding="utf-8"))
    target_functions = {
        "load_hlf_opn_candidate_genes",
        "hlf_opn_candidate_gene_summary",
        "hlf_opn_candidate_gene_cache_path",
        "hlf_opn_executable_candidate_inputs",
        "_load_offline_homology_evidence",
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
    assert "pcsec_pichia.services.gene_evidence" not in imports
    assert not {"run_blastp", "make_blast_db", "fetch_external_name_references", "build_gene_evidence_cache"} & called_names
