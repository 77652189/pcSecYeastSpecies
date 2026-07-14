from __future__ import annotations

import hashlib
import json
import math
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence

from pcsec_pichia.oe_capacity.schema import OECapacityValidationError


EXPECTED_TARGET_IDS = ("hLF", "OPN_ALPHA_FULL_PROJECT")
EXPECTED_CASE_KINDS = ("executable", "boundary")
ACCEPTANCE_CASE_MATRIX = (
    ("hLF", "executable", "PAS_chr2-1_0308"),
    ("hLF", "boundary", "PAS_chr1-4_0458"),
    ("OPN_ALPHA_FULL_PROJECT", "executable", "PAS_chr2-1_0308"),
    ("OPN_ALPHA_FULL_PROJECT", "boundary", "PAS_chr1-4_0458"),
)
CAPACITY_ASSET_RELATIVE_PATH = Path("Enzymedata") / "oe_capacity_baseline_capacity.json"
MANIFEST_NAME = "oe_capacity_manifest.json"


def verify_oe_capacity_run(
    run_dir: str | Path,
    expected_identity: Mapping[str, object] | None = None,
) -> dict[str, object]:
    """Verify one run from files on disk; caller-supplied pass/fail is ignored."""

    root = Path(run_dir).resolve()
    errors: list[str] = []
    warnings: list[str] = []
    manifest_path = root / MANIFEST_NAME
    manifest = _read_json_object(manifest_path, errors, "manifest")
    if not manifest:
        return _verification_payload(root, {}, [], errors, warnings, "")

    if manifest.get("schema_version") != 2:
        errors.append("manifest.schema_version")
    identity = manifest.get("run_identity")
    if not isinstance(identity, Mapping):
        errors.append("manifest.run_identity")
        identity = {}
    _verify_identity(identity, expected_identity or {}, errors)

    model = manifest.get("model")
    fingerprint = model.get("fingerprint") if isinstance(model, Mapping) else None
    if not _is_sha256(fingerprint):
        errors.append("manifest.model.fingerprint")
    expected_model = str((expected_identity or {}).get("model_fingerprint") or "")
    if expected_model and expected_model != str(fingerprint or ""):
        errors.append("identity.model_fingerprint_mismatch")

    capacity_asset = manifest.get("capacity_asset")
    if not isinstance(capacity_asset, Mapping):
        errors.append("manifest.capacity_asset")
    else:
        if capacity_asset.get("reviewed") is not True:
            errors.append("manifest.capacity_asset.reviewed")
        if not str(capacity_asset.get("version") or "").strip():
            errors.append("manifest.capacity_asset.version")
        if not _is_sha256(capacity_asset.get("sha256")):
            errors.append("manifest.capacity_asset.sha256")

    files = manifest.get("files")
    rows_path = _verified_file(root, files, "rows", errors)
    _verified_file(root, files, "report", errors)
    rows = _read_jsonl(rows_path, errors) if rows_path is not None else []
    if not rows:
        errors.append("rows.empty")

    _verify_status_and_coverage(manifest, rows, errors)
    _verify_row_identity(identity, rows, errors)
    _verify_scenario_and_proxy_evidence(manifest, rows, errors)
    _verify_case_contract(identity, rows, errors)

    manifest_sha256 = _sha256(manifest_path) if manifest_path.is_file() else ""
    coverage = _coverage_from_rows(rows)
    return _verification_payload(
        root,
        dict(identity),
        rows,
        errors,
        warnings,
        manifest_sha256,
        coverage,
        manifest,
    )


def run_phase2_oe_capacity_acceptance(
    repo_root: str | Path,
    output_root: str | Path,
) -> dict[str, object]:
    """Build the Phase 2 gate from trusted artifacts and real repository checks."""

    repo = Path(repo_root).resolve()
    output = Path(output_root).resolve()
    output.mkdir(parents=True, exist_ok=True)
    capacity_asset = _validate_capacity_asset(repo / CAPACITY_ASSET_RELATIVE_PATH)
    acceptance_run_id = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    smoke_root = (
        repo
        / "local_runs"
        / "oe_capacity"
        / "round6"
        / "acceptance_runs"
        / acceptance_run_id
        / "smokes"
    )
    smoke_commands: list[dict[str, object]] = []
    if capacity_asset["valid"]:
        smoke_commands = _run_smoke_matrix(repo, smoke_root)
    verifications: list[dict[str, object]] = []
    if smoke_root.is_dir():
        expected_by_run_name = {
            f"{target_id}-{case_kind}-{gene_id}": {
                "target_id": target_id,
                "case_kind": case_kind,
                "gene_id": gene_id,
            }
            for target_id, case_kind, gene_id in ACCEPTANCE_CASE_MATRIX
        }
        for manifest_path in sorted(smoke_root.glob(f"*/{MANIFEST_NAME}")):
            verifications.append(
                verify_oe_capacity_run(
                    manifest_path.parent,
                    expected_by_run_name.get(manifest_path.parent.name),
                )
            )

    required_cases: list[dict[str, object]] = []
    asset_hash = str(capacity_asset.get("sha256") or "")
    for target_id, case_kind, gene_id in ACCEPTANCE_CASE_MATRIX:
        matches = [
            item
            for item in verifications
            if _identity_matches(item.get("identity"), target_id, case_kind, gene_id)
        ]
        eligible = [
            item
            for item in matches
            if bool(item.get("passed"))
            and bool(capacity_asset["valid"])
            and _artifact_capacity_hash(item) == asset_hash
        ]
        passed = bool(eligible)
        required_cases.append(
            {
                "target_id": target_id,
                "case_kind": case_kind,
                "gene_id": gene_id,
                "passed": passed,
                "artifact_count": len(matches),
                "eligible_artifact_count": len(eligible),
                "artifacts": [str(item.get("run_dir") or "") for item in matches],
            }
        )

    commands: list[dict[str, object]] = []
    missing_information: list[str] = []
    if not capacity_asset["available"] or not capacity_asset["valid"]:
        missing_information.append("reviewed_baseline_capacity")
    else:
        commands = [*smoke_commands, *_run_acceptance_commands(repo)]

    artifact_asset_match = bool(required_cases) and all(
        bool(case["passed"]) for case in required_cases
    )
    if capacity_asset["valid"] and not artifact_asset_match:
        missing_information.append("capacity_asset_artifact_match")

    passed = (
        bool(required_cases)
        and all(bool(case["passed"]) for case in required_cases)
        and bool(capacity_asset["valid"])
        and artifact_asset_match
        and bool(commands)
        and all(int(command["exit_code"]) == 0 for command in commands)
    )
    summary: dict[str, object] = {
        "schema_version": 2,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "phase": "phase_2_gene_level_oe",
        "round": "round_6_hlf_opn_acceptance",
        "passed": passed,
        "expected_target_ids": list(EXPECTED_TARGET_IDS),
        "required_cases": required_cases,
        "artifact_verifications": verifications,
        "commands": commands,
        "inputs": {
            "repo_root": str(repo),
            "smoke_root": str(smoke_root),
            "acceptance_run_id": acceptance_run_id,
            "capacity_asset": capacity_asset,
        },
        "missing_information": sorted(set(missing_information)),
        "model_relative_only": True,
        "predicts_absolute_yield": False,
        "mutates_recommendation_tier": False,
    }
    _write_acceptance_outputs(summary, output)
    return summary


def _verify_identity(
    actual: Mapping[str, object],
    expected: Mapping[str, object],
    errors: list[str],
) -> None:
    if not str(actual.get("run_id") or "").strip():
        errors.append("manifest.run_identity.run_id")
    target_ids = actual.get("target_ids")
    context_ids = actual.get("context_ids")
    gene_ids = actual.get("gene_ids")
    if not _nonempty_strings(target_ids):
        errors.append("manifest.run_identity.target_ids")
    if not _nonempty_strings(context_ids):
        errors.append("manifest.run_identity.context_ids")
    if not _nonempty_strings(gene_ids):
        errors.append("manifest.run_identity.gene_ids")
    case_kind = str(actual.get("case_kind") or "")
    if case_kind not in {*EXPECTED_CASE_KINDS, "screen"}:
        errors.append("manifest.run_identity.case_kind")
    expected_target = str(expected.get("target_id") or "")
    if expected_target and expected_target not in (target_ids or []):
        errors.append("identity.target_id_mismatch")
    expected_case = str(expected.get("case_kind") or "")
    if expected_case and expected_case != case_kind:
        errors.append("identity.case_kind_mismatch")
    expected_gene = str(expected.get("gene_id") or "")
    if expected_gene and expected_gene not in (gene_ids or []):
        errors.append("identity.gene_id_mismatch")


def _verified_file(
    root: Path,
    files: object,
    key: str,
    errors: list[str],
) -> Path | None:
    if not isinstance(files, Mapping) or not isinstance(files.get(key), Mapping):
        errors.append(f"manifest.files.{key}")
        return None
    entry = files[key]
    relative = str(entry.get("path") or "")
    expected_hash = entry.get("sha256")
    path = (root / relative).resolve()
    try:
        path.relative_to(root)
    except ValueError:
        errors.append(f"files.{key}.path_escape")
        return None
    if not path.is_file():
        errors.append(f"files.{key}.missing")
        return None
    if not _is_sha256(expected_hash) or _sha256(path) != expected_hash:
        errors.append(f"files.{key}.sha256_mismatch")
    return path


def _verify_status_and_coverage(
    manifest: Mapping[str, object],
    rows: Sequence[Mapping[str, object]],
    errors: list[str],
) -> None:
    completed = sum(row.get("screen_status") == "completed" for row in rows)
    failures = len(rows) - completed
    status = manifest.get("status")
    if not isinstance(status, Mapping):
        errors.append("manifest.status")
    else:
        if status.get("completed_count") != completed:
            errors.append("manifest.status.completed_count")
        if status.get("failure_count") != failures:
            errors.append("manifest.status.failure_count")
        expected_state = (
            "partial_failure"
            if completed and failures
            else "failed"
            if failures
            else "completed"
        )
        if status.get("state") != expected_state:
            errors.append("manifest.status.state")
    if manifest.get("coverage") != _coverage_from_rows(rows):
        errors.append("manifest.coverage")


def _verify_row_identity(
    identity: Mapping[str, object],
    rows: Sequence[Mapping[str, object]],
    errors: list[str],
) -> None:
    targets = set(identity.get("target_ids") or [])
    contexts = set(identity.get("context_ids") or [])
    genes = set(identity.get("gene_ids") or [])
    if any(row.get("target_id") not in targets for row in rows):
        errors.append("rows.target_identity_mismatch")
    if any(row.get("context_id") not in contexts for row in rows):
        errors.append("rows.context_identity_mismatch")
    if any(row.get("gene_id") not in genes for row in rows):
        errors.append("rows.gene_identity_mismatch")


def _verify_scenario_and_proxy_evidence(
    manifest: Mapping[str, object],
    rows: Sequence[Mapping[str, object]],
    errors: list[str],
) -> None:
    config = manifest.get("config")
    config = config if isinstance(config, Mapping) else {}
    required = [str(item) for item in config.get("parameter_scenarios") or []]
    feature_enabled = config.get("feature_enabled") is True
    compare_proxy = config.get("compare_proxy") is True
    executable_rows = [
        row
        for row in rows
        if row.get("screen_status") == "completed"
        and row.get("execution_status") == "gene_level_executable"
    ]
    actual_scenarios_complete = not feature_enabled or (
        bool(required)
        and bool(executable_rows)
        and all(_has_scenario_pairs(row, required) for row in executable_rows)
    )
    if feature_enabled and executable_rows:
        if not actual_scenarios_complete:
            errors.append("rows.scenario_evidence_incomplete")
        if any(
            not _scenario_pairs_successful(row, required)
            for row in executable_rows
        ):
            errors.append("rows.completed_scenario_failure")
    scenario_manifest = manifest.get("scenario_completeness")
    if not isinstance(scenario_manifest, Mapping):
        errors.append("manifest.scenario_completeness")
    elif list(scenario_manifest.get("required") or []) != required:
        errors.append("manifest.scenario_completeness.required")
    elif scenario_manifest.get("complete") is not actual_scenarios_complete:
        errors.append("manifest.scenario_completeness.complete")
    if compare_proxy and executable_rows:
        if any(not _proxy_evidence_complete(row) for row in executable_rows):
            errors.append("rows.proxy_evidence_incomplete")
        if any(
            any(attempt.get("success") is not True for attempt in _proxy_attempts(row))
            for row in executable_rows
        ):
            errors.append("rows.completed_proxy_failure")


def _verify_case_contract(
    identity: Mapping[str, object],
    rows: Sequence[Mapping[str, object]],
    errors: list[str],
) -> None:
    case_kind = identity.get("case_kind")
    if case_kind == "executable":
        clean = bool(rows) and all(
            row.get("screen_status") == "completed"
            and row.get("execution_status") == "gene_level_executable"
            and not row.get("missing_information")
            and all(
                _finite(row.get(name))
                for name in (
                    "baseline_objective",
                    "proxy_objective",
                    "gene_capacity_objective",
                )
            )
            and _scenario_pairs_successful(
                row,
                [str(item) for item in row.get("uncertainty_scenarios") or []],
            )
            and bool(_proxy_attempts(row))
            and all(
                attempt.get("success") is True
                for attempt in _proxy_attempts(row)
            )
            for row in rows
        )
        if not clean:
            errors.append("case.executable_not_clean")
    elif case_kind == "boundary":
        preserved = bool(rows) and all(
            bool(row.get("missing_information"))
            and (
                row.get("screen_status") != "completed"
                or row.get("execution_status") != "gene_level_executable"
            )
            for row in rows
        )
        if not preserved:
            errors.append("case.boundary_not_preserved")


def _validate_capacity_asset(path: Path) -> dict[str, object]:
    result: dict[str, object] = {
        "path": str(path),
        "available": path.is_file(),
        "valid": False,
        "version": "",
        "sha256": "",
        "errors": [],
    }
    if not path.is_file():
        result["errors"] = ["capacity_asset.missing"]
        return result
    errors: list[str] = []
    payload = _read_json_object(path, errors, "capacity_asset")
    version = str(payload.get("asset_version") or "") if isinstance(payload, Mapping) else ""
    catalog = None
    if not errors:
        try:
            from pcsec_pichia.oe_capacity.parameters import (
                load_capacity_anchor_catalog,
            )
            catalog = load_capacity_anchor_catalog(path)
        except OECapacityValidationError as exc:
            errors.append(f"capacity_asset.invalid:{exc}")
    if not version:
        errors.append("capacity_asset.asset_version")
    if catalog is not None and not catalog.anchors:
        errors.append("capacity_asset.anchors_empty")
    result.update(
        {
            "valid": not errors,
            "version": version,
            "sha256": _sha256(path),
            "errors": errors,
        }
    )
    return result


def _run_acceptance_commands(repo: Path) -> list[dict[str, object]]:
    commands = (
        (
            "compileall",
            [sys.executable, "-m", "compileall", "-q", "python_pichia/src/pcsec_pichia/oe_capacity"],
        ),
        (
            "focused_tests",
            [
                sys.executable,
                "-m",
                "pytest",
                "-q",
                "python_pichia/tests/test_oe_capacity_acceptance.py",
                "python_pichia/tests/test_oe_capacity_contracts.py",
                "python_pichia/tests/test_oe_capacity_parameters.py",
                "python_pichia/tests/test_oe_capacity_constraints.py",
                "python_pichia/tests/test_oe_capacity_screen_reports.py",
            ],
        ),
        (
            "dependency_files",
            [
                "git",
                "diff",
                "--exit-code",
                "--",
                "requirements.txt",
                "python_pichia/pyproject.toml",
            ],
        ),
        (
            "ignored_outputs",
            ["git", "check-ignore", "local_runs/oe_capacity"],
        ),
    )
    results = [_run_command(repo, check_id, argv) for check_id, argv in commands]
    results.append(_run_protected_paths_check(repo))
    results.append(_run_sensitive_files_check(repo))
    return results


def _run_smoke_matrix(repo: Path, smoke_root: Path) -> list[dict[str, object]]:
    script = repo / "python_pichia" / "tools" / "run_oe_capacity_smoke_case.py"
    results: list[dict[str, object]] = []
    for target_id, case_kind, gene_id in ACCEPTANCE_CASE_MATRIX:
        run_name = f"{target_id}-{case_kind}-{gene_id}"
        results.append(
            _run_command(
                repo,
                f"smoke:{target_id}:{case_kind}",
                [
                    sys.executable,
                    str(script),
                    "--target-id",
                    target_id,
                    "--gene-id",
                    gene_id,
                    "--case-kind",
                    case_kind,
                    "--output-root",
                    str(smoke_root),
                    "--run-name",
                    run_name,
                ],
            )
        )
    return results


def _run_command(repo: Path, check_id: str, argv: list[str]) -> dict[str, object]:
    try:
        completed = subprocess.run(
            argv,
            cwd=repo,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            check=False,
        )
        output = (completed.stdout + completed.stderr).strip()
        return {
            "check_id": check_id,
            "argv": argv,
            "exit_code": completed.returncode,
            "output_tail": output[-4000:],
        }
    except OSError as exc:
        return {
            "check_id": check_id,
            "argv": argv,
            "exit_code": 127,
            "output_tail": str(exc),
        }


def _run_protected_paths_check(repo: Path) -> dict[str, object]:
    result = _run_command(
        repo,
        "protected_paths",
        ["git", "status", "--porcelain", "--", "Code", "Model", "Enzymedata", "Results"],
    )
    changed = [
        line[3:].strip().replace("\\", "/")
        for line in str(result.get("output_tail") or "").splitlines()
        if len(line) > 3
    ]
    disallowed = [
        path
        for path in changed
        if path != CAPACITY_ASSET_RELATIVE_PATH.as_posix()
    ]
    if int(result["exit_code"]) == 0 and disallowed:
        result["exit_code"] = 1
        result["output_tail"] = "unexpected protected path changes: " + ", ".join(disallowed)
    return result


def _run_sensitive_files_check(repo: Path) -> dict[str, object]:
    result = _run_command(
        repo,
        "sensitive_files",
        [
            "git",
            "ls-files",
            "--",
            ":(glob)**/.env",
            ":(glob)**/*.pem",
            ":(glob)**/*credential*",
        ],
    )
    if int(result["exit_code"]) == 0 and str(result.get("output_tail") or "").strip():
        result["exit_code"] = 1
    return result


def _write_acceptance_outputs(summary: Mapping[str, object], root: Path) -> None:
    json_path = root / "phase2_acceptance.json"
    markdown_path = root / "phase2_acceptance.md"
    json_path.write_text(
        json.dumps(dict(summary), ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    lines = [
        "# Phase 2 gene-level OE capacity acceptance",
        "",
        f"- Passed: `{bool(summary.get('passed'))}`",
        f"- Missing information: `{', '.join(summary.get('missing_information') or []) or 'none'}`",
        "",
        "## Required artifact cases",
        "",
    ]
    for case in summary.get("required_cases") or []:
        if isinstance(case, Mapping):
            lines.append(
                f"- `{case.get('target_id')}/{case.get('case_kind')}`: "
                f"`{case.get('passed')}` ({case.get('artifact_count')} artifacts)"
            )
    lines.extend(("", "## Executed checks", ""))
    for command in summary.get("commands") or []:
        if isinstance(command, Mapping):
            lines.append(
                f"- `{command.get('check_id')}`: exit `{command.get('exit_code')}`"
            )
    markdown_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _read_json_object(path: Path, errors: list[str], prefix: str) -> dict[str, object]:
    if not path.is_file():
        errors.append(f"{prefix}.missing")
        return {}
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError):
        errors.append(f"{prefix}.invalid_json")
        return {}
    if not isinstance(payload, dict):
        errors.append(f"{prefix}.not_object")
        return {}
    return payload


def _read_jsonl(path: Path, errors: list[str]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    try:
        for line in path.read_text(encoding="utf-8").splitlines():
            if not line.strip():
                continue
            payload = json.loads(line)
            if not isinstance(payload, dict):
                raise ValueError
            rows.append(payload)
    except (OSError, UnicodeError, json.JSONDecodeError, ValueError):
        errors.append("rows.invalid_jsonl")
    return rows


def _coverage_from_rows(rows: Sequence[Mapping[str, object]]) -> dict[str, object]:
    by_status: dict[str, int] = {}
    for row in rows:
        status = str(row.get("execution_status") or "unknown")
        by_status[status] = by_status.get(status, 0) + 1
    return {"total_rows": len(rows), "by_execution_status": by_status}


def _has_scenario_pairs(row: Mapping[str, object], required: Sequence[str]) -> bool:
    raw = row.get("scenario_results") or row.get("gene_capacity_scenario_results")
    if not isinstance(raw, list):
        return False
    found: set[str] = set()
    for item in raw:
        if not isinstance(item, Mapping):
            continue
        scenario = item.get("scenario") or item.get("parameter_scenario")
        if isinstance(scenario, Mapping):
            scenario = scenario.get("value")
        baseline = item.get("baseline") or item.get("baseline_snapshot")
        perturbed = item.get("perturbed") or item.get("perturbed_snapshot")
        failed = (
            isinstance(baseline, Mapping)
            and baseline.get("success") is not True
        ) or (
            isinstance(perturbed, Mapping)
            and perturbed.get("success") is not True
        )
        failure_preserved = not failed or bool(
            item.get("failure_reason")
            or (baseline.get("message") if isinstance(baseline, Mapping) else "")
            or (perturbed.get("message") if isinstance(perturbed, Mapping) else "")
        )
        if (
            isinstance(baseline, Mapping)
            and isinstance(perturbed, Mapping)
            and failure_preserved
        ):
            found.add(str(scenario))
    return set(required).issubset(found)


def _proxy_attempts(row: Mapping[str, object]) -> list[Mapping[str, object]]:
    raw = row.get("proxy_attempts") or row.get("proxy_results")
    return (
        [item for item in raw if isinstance(item, Mapping)]
        if isinstance(raw, list)
        else []
    )


def _scenario_pairs_successful(
    row: Mapping[str, object],
    required: Sequence[str],
) -> bool:
    raw = row.get("scenario_results") or row.get("gene_capacity_scenario_results")
    if not isinstance(raw, list) or not required:
        return False
    by_scenario: dict[str, Mapping[str, object]] = {}
    for item in raw:
        if not isinstance(item, Mapping):
            continue
        scenario = item.get("scenario") or item.get("parameter_scenario")
        if isinstance(scenario, Mapping):
            scenario = scenario.get("value")
        by_scenario[str(scenario)] = item
    for scenario in required:
        item = by_scenario.get(str(scenario))
        if item is None:
            return False
        baseline = item.get("baseline") or item.get("baseline_snapshot")
        perturbed = item.get("perturbed") or item.get("perturbed_snapshot")
        if not isinstance(baseline, Mapping) or not isinstance(perturbed, Mapping):
            return False
        if baseline.get("success") is not True or perturbed.get("success") is not True:
            return False
        if str(item.get("failure_reason") or ""):
            return False
    return True


def _proxy_evidence_complete(row: Mapping[str, object]) -> bool:
    attempts = _proxy_attempts(row)
    if not attempts:
        return False
    return all(
        isinstance(attempt, Mapping)
        and str(attempt.get("solver_status") or "").strip()
        and (
            attempt.get("success") is True
            or bool(str(attempt.get("message") or "").strip())
        )
        for attempt in attempts
    )


def _verification_payload(
    root: Path,
    identity: Mapping[str, object],
    rows: Sequence[Mapping[str, object]],
    errors: Sequence[str],
    warnings: Sequence[str],
    manifest_sha256: str,
    coverage: Mapping[str, object] | None = None,
    manifest: Mapping[str, object] | None = None,
) -> dict[str, object]:
    return {
        "passed": not errors,
        "run_dir": str(root),
        "identity": dict(identity),
        "coverage": dict(coverage or _coverage_from_rows(rows)),
        "manifest_sha256": manifest_sha256,
        "capacity_asset": dict((manifest or {}).get("capacity_asset") or {}),
        "errors": sorted(set(errors)),
        "warnings": list(warnings),
    }


def _identity_matches(
    value: object,
    target_id: str,
    case_kind: str,
    gene_id: str,
) -> bool:
    return (
        isinstance(value, Mapping)
        and target_id in (value.get("target_ids") or [])
        and gene_id in (value.get("gene_ids") or [])
        and value.get("case_kind") == case_kind
    )


def _artifact_capacity_hash(item: Mapping[str, object]) -> str:
    capacity = item.get("capacity_asset")
    return str(capacity.get("sha256") or "") if isinstance(capacity, Mapping) else ""


def _nonempty_strings(value: object) -> bool:
    return isinstance(value, list) and bool(value) and all(str(item).strip() for item in value)


def _is_sha256(value: object) -> bool:
    text = str(value or "")
    return len(text) == 64 and all(character in "0123456789abcdef" for character in text.lower())


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _finite(value: object) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(float(value))


def _finite_positive(value: object) -> bool:
    return _finite(value) and float(value) > 0


__all__ = [
    "CAPACITY_ASSET_RELATIVE_PATH",
    "EXPECTED_TARGET_IDS",
    "run_phase2_oe_capacity_acceptance",
    "verify_oe_capacity_run",
]
