from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pcsec_pichia.oe_capacity.acceptance as acceptance_module
from pcsec_pichia.oe_capacity import (
    run_phase2_oe_capacity_acceptance,
    verify_oe_capacity_run,
)


def test_verify_run_accepts_only_hash_bound_v2_artifacts(tmp_path: Path) -> None:
    run_dir = _write_run(tmp_path, target_id="hLF", case_kind="executable")

    verification = verify_oe_capacity_run(
        run_dir,
        {"target_id": "hLF", "case_kind": "executable"},
    )

    assert verification["passed"] is True
    assert verification["coverage"]["by_execution_status"]["gene_level_executable"] == 1
    assert verification["manifest_sha256"]
    assert verification["errors"] == []


def test_verify_run_rejects_tampered_rows_even_when_manifest_claims_completed(
    tmp_path: Path,
) -> None:
    run_dir = _write_run(tmp_path, target_id="hLF", case_kind="executable")
    rows_path = run_dir / "oe_capacity_rows.jsonl"
    rows_path.write_text(
        rows_path.read_text(encoding="utf-8").replace("1.2", "999.0"),
        encoding="utf-8",
    )

    verification = verify_oe_capacity_run(
        run_dir,
        {"target_id": "hLF", "case_kind": "executable"},
    )

    assert verification["passed"] is False
    assert "files.rows.sha256_mismatch" in verification["errors"]


def test_verify_run_rejects_scientifically_inconsistent_product_state(
    tmp_path: Path,
) -> None:
    run_dir = _write_run(tmp_path, target_id="hLF", case_kind="executable")
    rows_path = run_dir / "oe_capacity_rows.jsonl"
    row = json.loads(rows_path.read_text(encoding="utf-8"))
    row["product_state"] = "relative_uncalibrated"
    row["execution_mode"] = "relative_gene_capacity"
    rows_path.write_text(json.dumps(row, sort_keys=True) + "\n", encoding="utf-8")
    manifest_path = run_dir / "oe_capacity_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["files"]["rows"]["sha256"] = _sha256(rows_path)
    manifest_path.write_text(json.dumps(manifest, sort_keys=True), encoding="utf-8")

    verification = verify_oe_capacity_run(
        run_dir,
        {"target_id": "hLF", "case_kind": "executable"},
    )

    assert verification["passed"] is False
    assert "rows.product_state.relative_uncalibrated" in verification["errors"]


def test_verify_run_rejects_v1_and_missing_scenario_evidence(tmp_path: Path) -> None:
    run_dir = _write_run(tmp_path, target_id="hLF", case_kind="executable")
    manifest_path = run_dir / "oe_capacity_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["schema_version"] = 1
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    verification = verify_oe_capacity_run(run_dir, {"target_id": "hLF"})

    assert verification["passed"] is False
    assert "manifest.schema_version" in verification["errors"]

    run_dir = _write_run(
        tmp_path / "second",
        target_id="hLF",
        case_kind="executable",
        scenario_results=[],
    )
    verification = verify_oe_capacity_run(run_dir, {"target_id": "hLF"})
    assert verification["passed"] is False
    assert "rows.scenario_evidence_incomplete" in verification["errors"]


def test_verify_run_recomputes_coverage_and_model_identity(tmp_path: Path) -> None:
    run_dir = _write_run(tmp_path, target_id="hLF", case_kind="executable")
    manifest_path = run_dir / "oe_capacity_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["coverage"] = {"junk": 0}
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    verification = verify_oe_capacity_run(
        run_dir,
        {"target_id": "hLF", "model_fingerprint": "c" * 64},
    )

    assert verification["passed"] is False
    assert "manifest.coverage" in verification["errors"]
    assert "identity.model_fingerprint_mismatch" in verification["errors"]


def test_verify_run_rejects_wrong_gene_even_when_artifact_hashes_are_valid(
    tmp_path: Path,
) -> None:
    expected_gene = "PAS_chr2-1_0308"
    run_dir = _write_run(
        tmp_path,
        target_id="hLF",
        case_kind="executable",
        gene_id=expected_gene,
    )
    rows_path = run_dir / "oe_capacity_rows.jsonl"
    row = json.loads(rows_path.read_text(encoding="utf-8"))
    row["gene_id"] = "WRONG_GENE"
    rows_path.write_text(json.dumps(row, sort_keys=True) + "\n", encoding="utf-8")
    manifest_path = run_dir / "oe_capacity_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["run_identity"]["gene_ids"] = ["WRONG_GENE"]
    manifest["files"]["rows"]["sha256"] = _sha256(rows_path)
    manifest_path.write_text(json.dumps(manifest, sort_keys=True), encoding="utf-8")

    verification = verify_oe_capacity_run(
        run_dir,
        {
            "target_id": "hLF",
            "case_kind": "executable",
            "gene_id": expected_gene,
        },
    )

    assert verification["passed"] is False
    assert "identity.gene_id_mismatch" in verification["errors"]


def test_verify_run_rejects_completed_row_with_failed_solver_evidence(
    tmp_path: Path,
) -> None:
    run_dir = _write_run(
        tmp_path,
        target_id="hLF",
        case_kind="executable",
        scenario_results=[
            {
                "scenario": "nominal",
                "baseline": {"success": True, "solver_status": "optimal"},
                "perturbed": {
                    "success": False,
                    "solver_status": "infeasible",
                    "message": "infeasible",
                },
                "failure_reason": "scenario_perturbation_failed",
            }
        ],
    )

    verification = verify_oe_capacity_run(run_dir, {"target_id": "hLF"})

    assert verification["passed"] is False
    assert "rows.completed_scenario_failure" in verification["errors"]


def test_phase2_runner_uses_only_fresh_controlled_smokes(
    tmp_path: Path,
    monkeypatch,
) -> None:
    repo_root = tmp_path / "repo"
    _write_run(
        repo_root / "local_runs" / "oe_capacity" / "round6" / "smokes" / "forged",
        "hLF",
        "executable",
    )
    monkeypatch.setattr(
        acceptance_module,
        "_validate_capacity_asset",
        lambda _path: {
            "path": "asset.json",
            "available": True,
            "valid": True,
            "version": "v1",
            "sha256": "b" * 64,
            "errors": [],
        },
    )

    def fake_smokes(_repo: Path, smoke_root: Path):
        for target_id, case_kind, gene_id in acceptance_module.ACCEPTANCE_CASE_MATRIX:
            _write_run(
                smoke_root / f"{target_id}-{case_kind}-{gene_id}",
                target_id,
                case_kind,
                gene_id=gene_id,
            )
        return [{"check_id": "controlled_smokes", "exit_code": 0}]

    monkeypatch.setattr(acceptance_module, "_run_smoke_matrix", fake_smokes)
    monkeypatch.setattr(
        acceptance_module,
        "_run_acceptance_commands",
        lambda _repo: [{"check_id": "focused", "exit_code": 0}],
    )

    summary = run_phase2_oe_capacity_acceptance(repo_root, tmp_path / "acceptance")

    assert summary["passed"] is True
    assert all(
        "acceptance_runs" in item["run_dir"]
        for item in summary["artifact_verifications"]
    )
    assert all(case["eligible_artifact_count"] == 1 for case in summary["required_cases"])
    assert {
        (case["target_id"], case["case_kind"], case["gene_id"])
        for case in summary["required_cases"]
    } == set(acceptance_module.ACCEPTANCE_CASE_MATRIX)


def test_phase2_runner_stably_fails_without_reviewed_capacity_asset(
    tmp_path: Path,
) -> None:
    repo_root = tmp_path / "repo"
    smoke_root = repo_root / "local_runs" / "oe_capacity" / "round6" / "smokes"
    for target_id in ("hLF", "OPN_ALPHA_FULL_PROJECT"):
        _write_run(smoke_root / f"{target_id}-exec", target_id, "executable")
        _write_run(smoke_root / f"{target_id}-boundary", target_id, "boundary")

    summary = run_phase2_oe_capacity_acceptance(repo_root, tmp_path / "acceptance")

    assert summary["passed"] is False
    assert summary["phase"] == "phase_2_gene_level_oe"
    assert "reviewed_baseline_capacity" in summary["missing_information"]
    assert summary["inputs"]["capacity_asset"]["available"] is False
    assert (tmp_path / "acceptance" / "phase2_acceptance.json").is_file()


def _write_run(
    root: Path,
    target_id: str,
    case_kind: str,
    *,
    gene_id: str | None = None,
    scenario_results: list[dict[str, object]] | None = None,
) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    rows_path = root / "oe_capacity_rows.jsonl"
    report_path = root / "oe_capacity_report.md"
    scenarios = scenario_results
    if scenarios is None:
        scenarios = [
            {
                "scenario": "nominal",
                "baseline": {"success": True, "solver_status": "optimal"},
                "perturbed": {"success": True, "solver_status": "optimal"},
            }
        ]
    missing = [] if case_kind == "executable" else ["reviewed_baseline_capacity"]
    row = {
        "target_id": target_id,
        "context_id": "glucose_mu_0.1",
        "gene_id": gene_id or ("G1" if case_kind == "executable" else "G_BOUNDARY"),
        "screen_status": "completed" if case_kind == "executable" else "failed",
        "execution_mode": "comparison" if case_kind == "executable" else "not_executable",
        "execution_status": (
            "gene_level_executable" if case_kind == "executable" else "partial_mapping"
        ),
        "product_mode": "absolute_capacity",
        "product_state": (
            "absolute_available" if case_kind == "executable" else "absolute_unavailable"
        ),
        "absolute_capacity_availability": (
            "available_reviewed"
            if case_kind == "executable"
            else "unavailable_missing_reviewed_anchor"
        ),
        "calibration_status": (
            "reviewed_absolute" if case_kind == "executable" else "unavailable"
        ),
        "absolute_solver_allowed": case_kind == "executable",
        "model_fingerprint": "a" * 64,
        "baseline_objective": 1.0 if case_kind == "executable" else None,
        "proxy_objective": 1.1 if case_kind == "executable" else None,
        "gene_capacity_objective": 1.2 if case_kind == "executable" else None,
        "missing_information": missing,
        "uncertainty_scenarios": ["nominal"],
        "scenario_results": scenarios if case_kind == "executable" else [],
        "proxy_attempts": (
            [{"attempt_id": "R1", "success": True, "solver_status": "optimal"}]
            if case_kind == "executable"
            else []
        ),
    }
    rows_path.write_text(json.dumps(row, sort_keys=True) + "\n", encoding="utf-8")
    report_path.write_text("# verified report\n", encoding="utf-8")
    manifest = {
        "schema_version": 2,
        "run_identity": {
            "run_id": root.name,
            "target_ids": [target_id],
            "context_ids": ["glucose_mu_0.1"],
            "gene_ids": [row["gene_id"]],
            "case_kind": case_kind,
        },
        "model": {"fingerprint": "a" * 64},
        "capacity_asset": {
            "path": "Enzymedata/oe_capacity_baseline_capacity.json",
            "version": "fixture-v1",
            "sha256": "b" * 64,
            "reviewed": True,
        },
        "config": {
            "feature_enabled": True,
            "compare_proxy": True,
            "parameter_scenarios": ["nominal"],
        },
        "files": {
            "rows": {"path": rows_path.name, "sha256": _sha256(rows_path)},
            "report": {"path": report_path.name, "sha256": _sha256(report_path)},
        },
        "status": {
            "state": "completed" if case_kind == "executable" else "failed",
            "completed_count": 1 if case_kind == "executable" else 0,
            "failure_count": 0 if case_kind == "executable" else 1,
        },
        "scenario_completeness": {
            "required": ["nominal"],
            "complete": bool(scenarios) if case_kind == "executable" else False,
        },
        "coverage": {
            "total_rows": 1,
            "by_execution_status": {row["execution_status"]: 1},
            "by_product_state": {row["product_state"]: 1},
        },
    }
    (root / "oe_capacity_manifest.json").write_text(
        json.dumps(manifest, sort_keys=True), encoding="utf-8"
    )
    return root


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()
