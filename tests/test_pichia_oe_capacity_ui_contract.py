from __future__ import annotations

import ast
from pathlib import Path
from types import SimpleNamespace

from pcsec_pichia.oe_capacity import (
    ConfidenceLevel,
    OECapacityOutputs,
    OECapacityScreenResult,
    OECapacityScreenRow,
    OEDoseMode,
    OEExecutionMode,
    OEExecutionStatus,
    ParameterScenario,
)

from app.services import pichia_oe_capacity_service as service
from app.ui.common import OE_CAPACITY_PAGE
from pcsec_pichia.screens import prepare_screen_inputs


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_oe_capacity_page_is_registered_in_navigation_and_entrypoint() -> None:
    common = (REPO_ROOT / "app" / "ui" / "common.py").read_text(encoding="utf-8")
    entrypoint = (REPO_ROOT / "app" / "ui" / "streamlit_app.py").read_text(
        encoding="utf-8"
    )

    assert OE_CAPACITY_PAGE in common
    assert "render_oe_capacity" in entrypoint
    assert "elif page == OE_CAPACITY_PAGE" in entrypoint


def test_oe_capacity_view_uses_service_only_and_target_scoped_session_keys() -> None:
    view_path = REPO_ROOT / "app" / "ui" / "views" / "oe_capacity.py"
    source = view_path.read_text(encoding="utf-8")
    module_ast = ast.parse(source)
    imported_modules: set[str] = set()
    for node in ast.walk(module_ast):
        if isinstance(node, ast.ImportFrom) and node.module:
            imported_modules.add(node.module)
        elif isinstance(node, ast.Import):
            imported_modules.update(alias.name for alias in node.names)

    assert "app.services.pichia_oe_capacity_service" in imported_modules
    assert not any(
        module.startswith(("pcsec_pichia", "python_pichia"))
        for module in imported_modules
    )
    for text in (
        "oe_capacity_target_id",
        "oe_capacity_last_previews_by_target",
        "oe_capacity_last_runs_by_target",
        'f"oe_capacity_gene_id_{target_id}"',
        'f"oe_capacity_scenarios_{target_id}"',
        "nonzero_baseline_formation_flux",
        "baseline / proxy / gene-capacity",
        "不会自动修改 recommendation tier",
    ):
        assert text in source
    assert "st.cache_data" not in source
    assert "st.cache_resource" not in source


def test_service_cache_key_includes_target_context_and_uncertainty() -> None:
    source = (
        REPO_ROOT / "app" / "services" / "pichia_oe_capacity_service.py"
    ).read_text(encoding="utf-8")

    assert "@lru_cache(maxsize=8)" in source
    assert "def _prepare_runtime(\n    target_id: str," in source
    assert "growth_rate: float" in source
    assert "carbon_source_id: str" in source
    assert "relative_uncertainty: float" in source
    assert "DEFAULT_TARGET_IDS = (\"hLF\", \"OPN_ALPHA_FULL_PROJECT\")" in source
    assert "from pcsec_pichia.screens import prepare_screen_inputs" in source
    assert "_prepare_screen_inputs" not in source


def test_canonical_screen_preparation_has_a_public_service_adapter() -> None:
    assert callable(prepare_screen_inputs)


def test_service_submits_core_screen_and_writes_target_scoped_summary(
    tmp_path: Path,
    monkeypatch,
) -> None:
    runtime = SimpleNamespace(target_id="hLF")
    captured: dict[str, object] = {}

    def fake_screen(prepared, requests, config):
        captured["prepared"] = prepared
        captured["requests"] = requests
        captured["config"] = config
        row = OECapacityScreenRow(
            gene_id="G1",
            target_id="hLF",
            context_id="glucose_mu_0.1",
            execution_mode=OEExecutionMode.COMPARISON,
            execution_status=OEExecutionStatus.GENE_LEVEL_EXECUTABLE,
            dose_id="2x",
            dose_mode=OEDoseMode.EXPLICIT_MULTIPLIER,
            expression_multiplier=2.0,
            mapping_ids=("map-G1",),
            parameter_sources=("local_enzyme_data:test",),
            parameter_confidence=ConfidenceLevel.HIGH,
            uncertainty_scenarios=(ParameterScenario.NOMINAL,),
            baseline_objective=1.0,
            proxy_objective=1.1,
            gene_capacity_objective=1.2,
            gene_capacity_vs_baseline_delta=0.2,
            gene_capacity_vs_proxy_delta=0.1,
            protein_resource_cost_delta=0.01,
        )
        return OECapacityScreenResult(
            model_fingerprint="model-v1",
            config=config,
            rows=(row,),
        )

    def fake_write(result, output_dir):
        root = Path(output_dir)
        root.mkdir(parents=True)
        rows_path = root / "oe_capacity_rows.jsonl"
        manifest_path = root / "oe_capacity_manifest.json"
        report_path = root / "oe_capacity_report.md"
        rows_path.write_text("{}\n", encoding="utf-8")
        manifest_path.write_text("{}\n", encoding="utf-8")
        report_path.write_text("# report\n", encoding="utf-8")
        return OECapacityOutputs(
            output_dir=str(root),
            rows_path=str(rows_path),
            manifest_path=str(manifest_path),
            report_path=str(report_path),
        )

    monkeypatch.setattr(service, "_prepare_runtime", lambda *args: runtime)
    monkeypatch.setattr(service, "run_gene_level_oe_screen", fake_screen)
    monkeypatch.setattr(service, "write_oe_capacity_outputs", fake_write)
    summary = service.submit_oe_capacity_screen(
        target_id="hLF",
        gene_ids=("G1",),
        dose_payload={
            "dose_id": "2x",
            "dose_mode": "explicit_multiplier",
            "expression_multiplier": 2.0,
        },
        parameter_scenarios=("nominal",),
        run_name="hlf-ui-test",
        output_root=tmp_path,
    )

    assert captured["prepared"] is runtime
    assert captured["requests"][0].target_id == "hLF"  # type: ignore[index]
    assert captured["config"].parameter_scenarios == (  # type: ignore[union-attr]
        ParameterScenario.NOMINAL,
    )
    assert summary["target_id"] == "hLF"
    assert summary["completed_count"] == 1
    assert summary["failure_count"] == 0
    assert summary["model_relative_only"] is True
    assert summary["mutates_recommendation_tier"] is False
    assert (tmp_path / "hlf-ui-test" / "ui_run_summary.json").is_file()


def test_service_rejects_duplicate_run_name_without_overwriting(
    tmp_path: Path,
    monkeypatch,
) -> None:
    runtime = SimpleNamespace(target_id="hLF")
    existing = tmp_path / "duplicate-run"
    existing.mkdir()
    sentinel = existing / "keep.txt"
    sentinel.write_text("user-owned\n", encoding="utf-8")

    monkeypatch.setattr(service, "_prepare_runtime", lambda *args: runtime)
    monkeypatch.setattr(
        service,
        "run_gene_level_oe_screen",
        lambda prepared, requests, config: (_ for _ in ()).throw(
            AssertionError("duplicate run must be rejected before solving")
        ),
    )

    try:
        service.submit_oe_capacity_screen(
            target_id="hLF",
            gene_ids=("G1",),
            dose_payload={
                "dose_id": "categorical",
                "dose_mode": "categorical_only",
                "promoter": "unspecified",
            },
            parameter_scenarios=("nominal",),
            run_name="duplicate-run",
            output_root=tmp_path,
        )
    except FileExistsError as exc:
        assert "duplicate-run" in str(exc)
    else:
        raise AssertionError("duplicate run names must not overwrite existing output")

    assert sentinel.read_text(encoding="utf-8") == "user-owned\n"


def test_run_history_is_filtered_by_target(tmp_path: Path) -> None:
    for name, target_id in (("a-hlf", "hLF"), ("b-opn", "OPN_ALPHA_FULL_PROJECT")):
        run_dir = tmp_path / name
        run_dir.mkdir()
        (run_dir / "ui_run_summary.json").write_text(
            '{"run_name": "' + name + '", "target_id": "' + target_id + '"}\n',
            encoding="utf-8",
        )

    rows = service.list_oe_capacity_runs(tmp_path, target_id="hLF")

    assert [row["run_name"] for row in rows] == ["a-hlf"]
    assert all(row["target_id"] == "hLF" for row in rows)
