from __future__ import annotations

import ast
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_app_gene_id_standardization_facade_uses_persistent_cache(tmp_path, monkeypatch) -> None:
    from pcsec_pichia.core.paths import ProjectPaths

    from app.services import pichia_gene_catalog_service as service

    calls = {"count": 0}

    def fake_full_model_catalog(*, force_refresh: bool = False, paths=None) -> list[dict[str, object]]:
        calls["count"] += 1
        return [
            {
                "gene_id": "PAS_chr2-1_0140",
                "display_name": "KAR2",
                "standard_gene_symbol": "KAR2",
                "protein_name": "BiP molecular chaperone",
                "external_ids": {"uniprot": "C4R8K4"},
                "evidence_sources": ["UniProt"],
                "evidence_confidence": "high_exact_locus_tag",
                "ko_support_status": "ko_runnable_gpr_gene_deletion",
                "oe_support_status": "oe_runnable_reaction_proxy",
                "affected_reactions": ["sec_Kar2p_complex_formation"],
                "gpr_role": "single_gene",
            },
            {
                "gene_id": "AT250_GQ_6803479",
                "display_name": "AT250_GQ_6803479",
                "external_ids": {},
                "evidence_sources": [],
                "evidence_confidence": "low_model_only",
                "ko_support_status": "ko_no_gpr_effect",
                "oe_support_status": "oe_no_gpr_effect",
                "affected_reactions": [],
            },
        ]

    monkeypatch.setattr(service, "load_pichia_full_model_gene_catalog", fake_full_model_catalog)
    paths = ProjectPaths(repo_root=tmp_path)

    first_rows = service.load_pichia_gene_id_standardization(paths=paths)
    second_rows = service.load_pichia_gene_id_standardization(paths=paths)
    summary = service.pichia_gene_id_standardization_summary(paths=paths)

    assert calls["count"] == 1
    assert first_rows == second_rows
    assert first_rows[0]["gene_id"] == "PAS_chr2-1_0140"
    assert first_rows[0]["standard_symbol"] == "KAR2"
    assert first_rows[1]["annotation_sources"] == ["model_only"]
    assert summary["total_genes"] == 2
    assert summary["annotated_gene_count"] == 1
    assert summary["model_only_gene_ids"] == ["AT250_GQ_6803479"]
    assert service.pichia_gene_id_standardization_cache_path(paths).exists()
    assert "local_runs" in str(service.pichia_gene_id_standardization_cache_path(paths))
    assert "gene_catalog_cache" in str(service.pichia_gene_id_standardization_cache_path(paths))

    service.load_pichia_gene_id_standardization(force_refresh=True, paths=paths)

    assert calls["count"] == 2


def test_gene_id_standardization_summary_rebuilds_invalid_cache(tmp_path, monkeypatch) -> None:
    from pcsec_pichia.core.paths import ProjectPaths

    from app.services import pichia_gene_catalog_service as service

    def fake_full_model_catalog(*, force_refresh: bool = False, paths=None) -> list[dict[str, object]]:
        return [
            {
                "gene_id": "PAS_chr2-1_0140",
                "display_name": "KAR2",
                "external_ids": {"uniprot": "C4R8K4"},
                "evidence_sources": ["UniProt"],
                "ko_support_status": "ko_runnable_gpr_gene_deletion",
                "oe_support_status": "oe_no_gpr_effect",
            }
        ]

    monkeypatch.setattr(service, "load_pichia_full_model_gene_catalog", fake_full_model_catalog)
    paths = ProjectPaths(repo_root=tmp_path)
    cache_path = service.pichia_gene_id_standardization_cache_path(paths)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_text('{"schema_version": -1, "rows": []}', encoding="utf-8")

    summary = service.pichia_gene_id_standardization_summary(paths=paths)

    assert summary["total_genes"] == 1
    assert summary["annotated_gene_count"] == 1
    assert summary["model_only_count"] == 0


def test_gene_id_standardization_service_does_not_import_online_or_blast_clients() -> None:
    service_path = REPO_ROOT / "app" / "services" / "pichia_gene_catalog_service.py"
    module_ast = ast.parse(service_path.read_text(encoding="utf-8"))
    target_functions = {
        "load_pichia_gene_id_standardization",
        "pichia_gene_id_standardization_summary",
        "pichia_gene_id_standardization_cache_path",
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
