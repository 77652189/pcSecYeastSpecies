from __future__ import annotations

import ast
from pathlib import Path

import pandas as pd

from pcsec_pichia.core.paths import ProjectPaths
from app.services.results import PlotBuilder, ResultCatalog, ResultLoader


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_app_no_longer_imports_deleted_app_core_paths_module() -> None:
    deleted_path_imports: list[str] = []
    for path in (REPO_ROOT / "app").rglob("*.py"):
        module_ast = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(module_ast):
            imported: list[str] = []
            if isinstance(node, ast.ImportFrom) and node.module:
                imported.append(node.module)
            elif isinstance(node, ast.Import):
                imported.extend(alias.name for alias in node.names)
            for module_name in imported:
                if module_name == "app.core.paths":
                    deleted_path_imports.append(str(path.relative_to(REPO_ROOT)))

    assert deleted_path_imports == []


def test_project_paths_locate_required_directories() -> None:
    paths = ProjectPaths.discover()

    assert paths.repo_root.name == "pcSecYeastSpecies"
    assert paths.results_dir.exists()
    assert paths.models_dir.joinpath("pcSecYeast.mat").exists()
    assert paths.local_runs_dir.exists()


def test_result_catalog_finds_xlsx_and_mat_datasets() -> None:
    paths = ProjectPaths.discover()
    datasets = ResultCatalog(paths).list_datasets()

    assert any(dataset.id == "Results/CSource/CSource_res.xlsx" for dataset in datasets)
    assert any(dataset.suffix == ".mat" for dataset in datasets)
    assert any(dataset.species == "SCE" for dataset in datasets)


def test_loader_reads_csource_excel_and_builds_chart() -> None:
    paths = ProjectPaths.discover()
    catalog = ResultCatalog(paths)
    dataset = catalog.get_dataset("Results/CSource/CSource_res.xlsx")

    loaded = ResultLoader().load_dataset(dataset)

    assert loaded.kind == "table"
    assert {"Kma", "Sce", "Ppa"}.issubset(loaded.tables)
    assert pd.DataFrame(loaded.tables["Sce"]).shape[1] >= 4
    assert PlotBuilder().build_chart(loaded, "Sce") is not None


def test_loader_summarizes_mat_dataset() -> None:
    paths = ProjectPaths.discover()
    catalog = ResultCatalog(paths)
    dataset = catalog.get_dataset("Results/Growth_rate_TP/fluxesSCE_ref.mat")

    loaded = ResultLoader().load_dataset(dataset)

    assert loaded.kind == "mat"
    assert any(item["变量"] == "fluxes" for item in loaded.variable_summary)
    assert "fluxes" in loaded.tables
