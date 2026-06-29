from __future__ import annotations

import ast
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]

REVIEW_PACKAGES = {
    "formal_engine": [
        "python_pichia/src/pcsec_pichia/loading/__init__.py",
        "python_pichia/src/pcsec_pichia/targets/__init__.py",
        "python_pichia/src/pcsec_pichia/secretion_plan/__init__.py",
        "python_pichia/src/pcsec_pichia/constraints/_prototype_adapter.py",
        "python_pichia/src/pcsec_pichia/simulation/__init__.py",
        "python_pichia/src/pcsec_pichia/screens/__init__.py",
        "python_pichia/src/pcsec_pichia/screens/gene_perturbation_map.py",
        "python_pichia/src/pcsec_pichia/reports/__init__.py",
        "python_pichia/src/pcsec_pichia/alignment/__init__.py",
        "python_pichia/src/pcsec_pichia/analysis/__init__.py",
        "python_pichia/src/pcsec_pichia/pipeline.py",
    ],
    "app_service_facade": [
        "app/services/pichia_secretion_service.py",
        "app/services/pichia_secretion_schema.py",
        "app/services/pichia_secretion_runner.py",
        "app/services/pichia_request_mapping_service.py",
        "app/services/pichia_background_tasks.py",
        "app/services/pichia_gene_catalog_service.py",
        "app/services/pichia_screen_preview_service.py",
        "app/services/pichia_target_catalog_service.py",
    ],
    "streamlit_draft_ui": [
        "app/ui/streamlit_app.py",
        "app/ui/views/simulation.py",
        "app/ui/views/simulation_builder.py",
        "app/ui/views/simulation_display.py",
        "app/ui/views/simulation_gene_catalog.py",
        "app/ui/views/simulation_gene_inputs.py",
        "app/ui/views/simulation_gene_text.py",
        "app/ui/views/simulation_results.py",
        "app/ui/views/candidate_path_graph.py",
        "app/ui/views/simulation_matlab_reference.py",
    ],
    "experimental_api": [
        "app/api/__init__.py",
        "app/api/pichia_secretion_api.py",
    ],
    "startup_scripts": [
        "run_streamlit.ps1",
        "scripts/start_pcSecYeastSpecies_lan.ps1",
        "scripts/repair_pcSecYeastSpecies_desktop_shortcut.ps1",
    ],
    "active_docs": [
        "docs/pichia_python_next_development_slices_2026-06-26.md",
        "docs/pichia_python_architecture.md",
        "docs/pichia_python_release_validation_2026-06-25.md",
    ],
    "boundary_tests": [
        "tests/test_pichia_secretion_service_contract.py",
        "tests/test_pichia_fastapi_entrypoints.py",
        "tests/test_streamlit_startup_scripts.py",
        "tests/test_docs_active_boundary.py",
        "tests/test_review_package_boundaries.py",
        "tests/test_slow_test_gates.py",
    ],
}


def test_review_packages_have_their_expected_anchor_files() -> None:
    missing_by_package: dict[str, list[str]] = {}

    for package_name, relative_paths in REVIEW_PACKAGES.items():
        missing = [
            relative_path
            for relative_path in relative_paths
            if not (REPO_ROOT / relative_path).exists()
        ]
        if missing:
            missing_by_package[package_name] = missing

    assert missing_by_package == {}


def test_pichia_app_orchestration_files_stay_outside_formal_engine_package() -> None:
    pichia_app_files = [
        path
        for relative_paths in (
            REVIEW_PACKAGES["app_service_facade"],
            REVIEW_PACKAGES["streamlit_draft_ui"],
            REVIEW_PACKAGES["experimental_api"],
        )
        for path in relative_paths
    ]

    assert all(not path.startswith("python_pichia/") for path in pichia_app_files)


def test_all_pichia_app_service_modules_are_in_review_package() -> None:
    discovered = {
        str(path.relative_to(REPO_ROOT)).replace("\\", "/")
        for path in (REPO_ROOT / "app" / "services").glob("pichia_*.py")
    }
    reviewed = {
        path.replace("\\", "/")
        for path in REVIEW_PACKAGES["app_service_facade"]
    }

    assert discovered <= reviewed


def test_formal_engine_review_package_stays_under_python_pichia() -> None:
    assert all(
        path.startswith("python_pichia/src/pcsec_pichia/")
        for path in REVIEW_PACKAGES["formal_engine"]
    )


def test_formal_gene_catalog_service_does_not_import_probe_modules_directly() -> None:
    source_path = REPO_ROOT / "python_pichia" / "src" / "pcsec_pichia" / "services" / "gene_catalog.py"
    source = source_path.read_text(encoding="utf-8")
    module_ast = ast.parse(source)
    imported_modules: list[str] = []
    for node in ast.walk(module_ast):
        if isinstance(node, ast.ImportFrom) and node.module:
            imported_modules.append(node.module)
        elif isinstance(node, ast.Import):
            imported_modules.extend(alias.name for alias in node.names)

    assert not any(module_name.startswith("pcsec_pichia.probe") for module_name in imported_modules)


def test_direct_probe_imports_in_formal_engine_are_explicit_migration_debt() -> None:
    expected_probe_backed_modules = {
        "python_pichia/src/pcsec_pichia/constraints/_prototype_adapter.py",
        "python_pichia/src/pcsec_pichia/loading/__init__.py",
        "python_pichia/src/pcsec_pichia/reports/_prototype_adapter.py",
        "python_pichia/src/pcsec_pichia/screens/_prototype_adapter.py",
        "python_pichia/src/pcsec_pichia/secretion_plan/_prototype_adapter.py",
        "python_pichia/src/pcsec_pichia/simulation/__init__.py",
        "python_pichia/src/pcsec_pichia/targets/__init__.py",
    }
    discovered_probe_backed_modules: set[str] = set()
    for source_path in (REPO_ROOT / "python_pichia" / "src" / "pcsec_pichia").rglob("*.py"):
        relative_path = str(source_path.relative_to(REPO_ROOT)).replace("\\", "/")
        if "/probe/" in relative_path:
            continue
        module_ast = ast.parse(source_path.read_text(encoding="utf-8-sig"))
        imported_modules: list[str] = []
        for node in ast.walk(module_ast):
            if isinstance(node, ast.ImportFrom) and node.module:
                imported_modules.append(node.module)
            elif isinstance(node, ast.Import):
                imported_modules.extend(alias.name for alias in node.names)
        if any(module_name.startswith("pcsec_pichia.probe") for module_name in imported_modules):
            discovered_probe_backed_modules.add(relative_path)

    assert discovered_probe_backed_modules == expected_probe_backed_modules
