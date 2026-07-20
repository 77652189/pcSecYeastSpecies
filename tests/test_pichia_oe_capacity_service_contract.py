from __future__ import annotations

from pathlib import Path
from threading import Event, Thread
from types import SimpleNamespace

from pcsec_pichia.oe_capacity import (
    AbsoluteCapacityAvailability,
    ConfidenceLevel,
    EvidenceSourceType,
    GeneCapacityCatalog,
    GeneEnzymeReactionMapping,
    GPRRole,
    OECapacityOutputs,
    OECapacityScreenResult,
    OECapacityScreenRow,
    OEDoseMode,
    OEExecutionMode,
    OEExecutionStatus,
    OECalibrationStatus,
    OEProductMode,
    OEProductState,
    ParameterScenario,
    ParameterPolicy,
)

from app.services import pichia_oe_capacity_service as service
from pcsec_pichia.screens import prepare_screen_inputs


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_preview_uses_core_product_summary_without_mutating_mapping_status(
    monkeypatch,
) -> None:
    mapping = GeneEnzymeReactionMapping(
        mapping_id="map-G1",
        model_fingerprint="model-v1",
        gene_id="G1",
        enzyme_id="R1_complex",
        reaction_id="R1",
        gpr_rule="G1",
        gpr_role=GPRRole.SINGLE_GENE,
        enzyme_variable_id="R1_complex_formation",
        formation_or_dilution_reaction_id="R1_complex_formation",
        mapping_source=EvidenceSourceType.CURRENT_MODEL,
        mapping_confidence=ConfidenceLevel.HIGH,
        execution_status=OEExecutionStatus.GENE_LEVEL_EXECUTABLE,
    )
    catalog = GeneCapacityCatalog(
        model_fingerprint="model-v1",
        mappings=(mapping,),
    )
    runtime = SimpleNamespace(
        fixed_model=object(),
        gene_capacity_catalog=catalog,
        parameter_policy=ParameterPolicy(parameter_sets=()),
    )
    monkeypatch.setattr(service, "_prepare_runtime", lambda *args: runtime)
    core_plan = SimpleNamespace(executable_capacity_specs=())
    monkeypatch.setattr(
        service,
        "plan_gene_level_overexpression",
        lambda *args: core_plan,
    )
    monkeypatch.setattr(
        service,
        "resolve_oe_product_plan",
        lambda plan, **kwargs: plan,
    )
    monkeypatch.setattr(
        service,
        "summarize_oe_product_candidate",
        lambda plan: {
            "product_state": "absolute_unavailable",
            "absolute_capacity_availability": "unavailable_missing_reviewed_anchor",
            "missing_information": ["reviewed_baseline_capacity"],
        },
    )

    preview = service.preview_oe_capacity_candidate(target_id="hLF", gene_id="G1")

    assert preview["parameter_set_count"] == 0
    assert preview["executable_mapping_count"] == 0
    assert preview["mappings"][0]["execution_status"] == "gene_level_executable"
    assert preview["product"]["product_state"] == "absolute_unavailable"
    assert "reviewed_baseline_capacity" in preview["product"]["missing_information"]


def test_service_cache_key_includes_target_context_and_uncertainty() -> None:
    source = (
        REPO_ROOT / "app" / "services" / "pichia_oe_capacity_service.py"
    ).read_text(encoding="utf-8")

    assert "@lru_cache(maxsize=8)" in source
    assert "def _prepare_runtime(\n    target_id: str," in source
    assert "growth_rate: float" in source
    assert "carbon_source_id: str" in source
    assert "relative_uncertainty: float" in source
    assert "capacity_asset_version: str" in source
    assert "_capacity_asset_version()" in source
    assert "DEFAULT_TARGET_IDS = (\"hLF\", \"OPN_ALPHA_FULL_PROJECT\")" in source
    assert "from pcsec_pichia.oe_capacity.external_candidate_audit import (" in source
    assert "prepare_external_candidate_runtime" in source
    assert "from pcsec_pichia.screens" not in source
    assert "from pcsec_pichia.loading" not in source
    assert "_prepare_screen_inputs" not in source


def test_candidate_review_missing_cache_is_safe_and_promotion_requires_explicit_approval(
    tmp_path: Path,
) -> None:
    review = service.load_oe_capacity_candidate_review(
        tmp_path / "missing",
        target_id="hLF",
    )
    assert review["available"] is False
    assert review["candidates"] == []
    assert "正式容量资产保持不变" in review["message"]

    try:
        service.promote_oe_capacity_candidate_selection(
            candidate_root=tmp_path / "missing",
            candidate_ids=("candidate-1",),
            reviewer="reviewer",
            expected_candidate_manifest_sha256="a" * 64,
            expected_asset_sha256="b" * 64,
            explicit_approval=False,
        )
    except ValueError as exc:
        assert "explicit_approval=True" in str(exc)
    else:
        raise AssertionError("formal promotion must require explicit approval")


def test_candidate_promotion_keeps_the_cached_runtime_resolver_seam(
    tmp_path: Path,
    monkeypatch,
) -> None:
    captured: dict[str, object] = {}

    def fake_promote(repo_root, **kwargs):
        captured["repo_root"] = repo_root
        captured.update(kwargs)
        return {"decision": "approved", "acceptance_started": False}

    monkeypatch.setattr(service, "promote_external_candidate_selection", fake_promote)

    result = service.promote_oe_capacity_candidate_selection(
        candidate_root=tmp_path / "candidates",
        candidate_ids=("candidate-1",),
        reviewer="reviewer",
        expected_candidate_manifest_sha256="a" * 64,
        expected_asset_sha256="b" * 64,
        explicit_approval=True,
    )

    assert result == {"decision": "approved", "acceptance_started": False}
    assert captured["runtime_resolver"] is service._prepare_runtime
    assert captured["capacity_asset_path"] == service.OE_CAPACITY_ASSET_PATH


def test_canonical_screen_preparation_has_a_public_service_adapter() -> None:
    assert callable(prepare_screen_inputs)


def test_service_submits_core_screen_and_writes_target_scoped_summary(
    tmp_path: Path,
    monkeypatch,
) -> None:
    runtime = SimpleNamespace(
        target_id="hLF",
        capacity_asset_metadata={
            "path": "asset.json",
            "version": "pending-rd-review",
            "sha256": "a" * 64,
            "reviewed": False,
        },
    )
    captured: dict[str, object] = {}

    def fake_screen(prepared, requests, config):
        status = service._load_run_status(tmp_path / "hlf-ui-test")
        assert status["status"] == "running"
        captured["prepared"] = prepared
        captured["requests"] = requests
        captured["config"] = config
        row = OECapacityScreenRow(
            gene_id="G1",
            target_id="hLF",
            context_id="glucose_mu_0.1",
            execution_mode=OEExecutionMode.COMPARISON,
            execution_status=OEExecutionStatus.GENE_LEVEL_EXECUTABLE,
            product_mode=OEProductMode.ABSOLUTE_CAPACITY,
            product_state=OEProductState.ABSOLUTE_AVAILABLE,
            absolute_capacity_availability=AbsoluteCapacityAvailability.AVAILABLE_REVIEWED,
            calibration_status=OECalibrationStatus.REVIEWED_ABSOLUTE,
            absolute_solver_allowed=True,
            model_fingerprint="model-v1",
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
            nominal_capacity=2.0,
            limitations=("model_relative_only", "no_mg_per_litre_prediction"),
        )
        return OECapacityScreenResult(
            model_fingerprint="model-v1",
            config=config,
            rows=(row,),
        )

    def fake_write(result, output_dir, **metadata):
        root = Path(output_dir)
        assert root.is_dir()
        assert metadata["run_identity"]["case_kind"] == "screen"
        assert metadata["run_identity"]["gene_ids"] == ["G1"]
        assert metadata["capacity_asset"]["reviewed"] is False
        assert metadata["capacity_asset"]["sha256"] == "a" * 64
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
    assert summary["status"] == "completed"
    assert summary["model_relative_only"] is False
    assert summary["absolute_capacity_available"] is True
    assert summary["product_states"] == ["absolute_available"]
    assert summary["mutates_recommendation_tier"] is False
    assert (tmp_path / "hlf-ui-test" / "ui_run_summary.json").is_file()
    status = service._load_run_status(tmp_path / "hlf-ui-test")
    assert status["status"] == "completed"
    assert status["completed_count"] == 1


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


def test_service_preserves_failed_run_status_when_solver_raises(
    tmp_path: Path,
    monkeypatch,
) -> None:
    monkeypatch.setattr(
        service,
        "_prepare_runtime",
        lambda *args: (_ for _ in ()).throw(RuntimeError("baseline failed")),
    )

    try:
        service.submit_oe_capacity_screen(
            target_id="hLF",
            gene_ids=("G1",),
            dose_payload={
                "dose_id": "2x",
                "dose_mode": "explicit_multiplier",
                "expression_multiplier": 2.0,
            },
            parameter_scenarios=("nominal",),
            run_name="failed-run",
            output_root=tmp_path,
        )
    except RuntimeError as exc:
        assert "baseline failed" in str(exc)
    else:
        raise AssertionError("runtime failure must propagate")

    status = service._load_run_status(tmp_path / "failed-run")
    assert status["status"] == "failed"
    assert status["error_type"] == "RuntimeError"
    assert status["error_message"] == "baseline failed"
    loaded = service.load_oe_capacity_run(tmp_path / "failed-run")
    assert loaded["available"] is False
    assert loaded["target_id"] == "hLF"


def test_run_directory_is_reserved_before_runtime_preparation(
    tmp_path: Path,
    monkeypatch,
) -> None:
    entered_runtime = Event()
    release_runtime = Event()
    first_errors: list[Exception] = []

    def blocked_runtime(*args):
        entered_runtime.set()
        assert release_runtime.wait(timeout=5)
        raise RuntimeError("stop first run")

    monkeypatch.setattr(service, "_prepare_runtime", blocked_runtime)
    kwargs = {
        "target_id": "hLF",
        "gene_ids": ("G1",),
        "dose_payload": {
            "dose_id": "2x",
            "dose_mode": "explicit_multiplier",
            "expression_multiplier": 2.0,
        },
        "parameter_scenarios": ("nominal",),
        "run_name": "concurrent-run",
        "output_root": tmp_path,
    }

    def first_run() -> None:
        try:
            service.submit_oe_capacity_screen(**kwargs)
        except Exception as exc:
            first_errors.append(exc)

    thread = Thread(target=first_run)
    thread.start()
    assert entered_runtime.wait(timeout=5)
    try:
        service.submit_oe_capacity_screen(**kwargs)
    except FileExistsError:
        pass
    else:
        raise AssertionError("a concurrent run must not reuse a reserved directory")
    finally:
        release_runtime.set()
        thread.join(timeout=5)

    assert not thread.is_alive()
    assert len(first_errors) == 1
    assert isinstance(first_errors[0], RuntimeError)
