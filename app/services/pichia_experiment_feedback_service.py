from __future__ import annotations

import json
import re
from dataclasses import asdict, replace
from enum import Enum
from pathlib import Path
from typing import Any, Mapping

from app import ensure_python_pichia_on_path

ensure_python_pichia_on_path()

from pcsec_pichia.experimental_feedback import (
    CalibrationConfig,
    PredictionLinkageResult,
    build_calibration_summary,
    build_prediction_index,
    link_experiments_to_predictions,
    load_experiment_bundle,
    validate_experiment_bundle,
    write_calibration_outputs,
    write_experiment_feedback_cache,
    write_prediction_experiment_report_outputs,
)


DEFAULT_EXPERIMENT_FEEDBACK_ROOT = Path("local_runs") / "experiment_feedback" / "ui_runs"
_EXPERIMENT_SUFFIXES = {".csv", ".xlsx", ".jsonl"}


def submit_experiment_feedback_import(
    *,
    experiment_filename: str,
    experiment_bytes: bytes,
    prediction_filename: str = "",
    prediction_bytes: bytes | None = None,
    prediction_path: str | Path | None = None,
    run_name: str = "",
    output_root: str | Path = DEFAULT_EXPERIMENT_FEEDBACK_ROOT,
    calibration_config: CalibrationConfig | Mapping[str, Any] | None = None,
    experiment_metadata: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    safe_run_name = _safe_run_name(run_name)
    run_dir = Path(output_root) / safe_run_name
    if run_dir.exists():
        raise FileExistsError(f"experiment feedback run already exists: {run_dir}")
    inbox_dir = run_dir / "inbox"
    inbox_dir.mkdir(parents=True, exist_ok=False)
    experiment_path = inbox_dir / _safe_upload_name(experiment_filename, _EXPERIMENT_SUFFIXES)
    experiment_path.write_bytes(experiment_bytes)
    prediction_payload: dict[str, Any] | None = None
    prediction_source = ""
    if prediction_bytes is not None:
        prediction_upload = inbox_dir / _safe_upload_name(prediction_filename or "fact_pack.json", {".json"})
        prediction_upload.write_bytes(prediction_bytes)
        prediction_payload = _read_json_object(prediction_upload)
        prediction_source = str(prediction_upload)
    elif prediction_path:
        resolved_prediction = Path(prediction_path)
        prediction_payload = _read_json_object(resolved_prediction)
        prediction_source = str(resolved_prediction)

    bundle = load_experiment_bundle(experiment_path, metadata=experiment_metadata)
    initial_validation = validate_experiment_bundle(bundle)
    prediction_index = build_prediction_index((prediction_payload,)) if prediction_payload else build_prediction_index(())
    linkage = link_experiments_to_predictions(bundle, prediction_index)
    linked_bundle = replace(bundle, prediction_links=linkage.links)
    validation = validate_experiment_bundle(linked_bundle)
    cache_outputs = write_experiment_feedback_cache(linked_bundle, run_dir / "validated")
    linkage_dir = run_dir / "linkage"
    linkage_dir.mkdir(parents=True, exist_ok=True)
    linkage_path = linkage_dir / "linkage_summary.json"
    linkage_path.write_text(
        json.dumps(_linkage_payload(linkage), ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    calibration_payload: dict[str, Any]
    calibration_paths: dict[str, str] = {}
    if validation.is_valid:
        config = _calibration_config(calibration_config)
        calibration = build_calibration_summary(validation, linkage, config)
        calibration_outputs = write_calibration_outputs(calibration, run_dir / "calibration")
        calibration_payload = _calibration_payload(calibration)
        calibration_paths = {
            "records_path": str(calibration_outputs.records_path),
            "summary_path": str(calibration_outputs.summary_path),
            "manifest_path": str(calibration_outputs.manifest_path),
        }
        report_outputs = write_prediction_experiment_report_outputs(
            calibration=calibration,
            validation=validation,
            linkage=linkage,
            output_dir=run_dir / "report",
            experiment_path=experiment_path,
            prediction_runs=(prediction_payload,) if prediction_payload else (),
            prediction_sources=(prediction_source,) if prediction_source else (),
            source_classification="local_unreviewed_input",
            supporting_files={
                "validated_records": cache_outputs.validated_records_path.name,
                "import_manifest": cache_outputs.manifest_path.name,
                "linkage": linkage_path.name,
                "calibration_records": calibration_outputs.records_path.name,
                "calibration_summary": calibration_outputs.summary_path.name,
                "calibration_manifest": calibration_outputs.manifest_path.name,
            },
        )
        calibration_paths.update(
            {
                "report_manifest_path": str(report_outputs.manifest_path),
                "report_summary_path": str(report_outputs.summary_path),
                "report_path": str(report_outputs.report_path),
            }
        )
    else:
        calibration_payload = {
            "available": False,
            "reason": "validation_failed",
            "targets": [],
            "records": [],
        }

    summary = {
        "run_name": safe_run_name,
        "run_dir": str(run_dir),
        "experiment_source": str(experiment_path),
        "prediction_source": prediction_source,
        "validation": _validation_payload(validation),
        "initial_validation": _validation_payload(initial_validation),
        "linkage": _linkage_payload(linkage),
        "calibration": calibration_payload,
        "paths": {
            "validated_records_path": str(cache_outputs.validated_records_path),
            "conflicts_path": str(cache_outputs.conflicts_path),
            "warnings_path": str(cache_outputs.warnings_path),
            "manifest_path": str(cache_outputs.manifest_path),
            "linkage_path": str(linkage_path),
            **calibration_paths,
        },
        "import_is_calibration": False,
    }
    (run_dir / "run_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return summary


def load_experiment_feedback_run(run_dir: str | Path) -> dict[str, Any]:
    root = Path(run_dir)
    summary_path = root / "run_summary.json"
    if not summary_path.exists():
        return {
            "run_name": root.name,
            "run_dir": str(root),
            "available": False,
            "validation": {"is_valid": False, "errors": [], "warnings": []},
            "linkage": {},
            "calibration": {"available": False, "targets": [], "records": []},
            "paths": {},
        }
    payload = _read_json_object(summary_path)
    payload["available"] = True
    return payload


def list_experiment_feedback_runs(
    output_root: str | Path = DEFAULT_EXPERIMENT_FEEDBACK_ROOT,
) -> list[dict[str, Any]]:
    root = Path(output_root)
    if not root.exists():
        return []
    rows = [load_experiment_feedback_run(path) for path in root.iterdir() if path.is_dir()]
    return sorted(rows, key=lambda row: str(row.get("run_name") or ""), reverse=True)


def list_prediction_fact_packs(search_root: str | Path = "local_runs") -> list[str]:
    root = Path(search_root)
    if not root.exists():
        return []
    return [str(path) for path in sorted(root.rglob("fact_pack.json"), reverse=True)]


def export_experiment_feedback_issues(
    run_dir: str | Path,
    *,
    issue_kind: str,
) -> bytes:
    filename = {"conflicts": "conflicts.jsonl", "warnings": "warnings.jsonl"}.get(issue_kind)
    if filename is None:
        raise ValueError(f"unsupported issue_kind: {issue_kind}")
    path = Path(run_dir) / "validated" / filename
    return path.read_bytes() if path.exists() else b""


def export_experiment_feedback_report(run_dir: str | Path) -> bytes:
    path = Path(run_dir) / "report" / "prediction_experiment_report.md"
    return path.read_bytes() if path.is_file() else b""


def _validation_payload(validation: object) -> dict[str, Any]:
    return {
        "is_valid": bool(validation.is_valid),  # type: ignore[attr-defined]
        "errors": [_json_ready(asdict(issue)) for issue in validation.errors],  # type: ignore[attr-defined]
        "warnings": [_json_ready(asdict(issue)) for issue in validation.warnings],  # type: ignore[attr-defined]
    }


def _linkage_payload(linkage: PredictionLinkageResult) -> dict[str, Any]:
    return {
        "matched_count": linkage.matched_count,
        "ambiguous_count": linkage.ambiguous_count,
        "missing_prediction_count": linkage.missing_prediction_count,
        "context_mismatch_count": linkage.context_mismatch_count,
        "control_count": linkage.control_count,
        "links": [_json_ready(asdict(link)) for link in linkage.links],
    }


def _calibration_payload(calibration: object) -> dict[str, Any]:
    return {
        "available": True,
        "config": _json_ready(asdict(calibration.config)),  # type: ignore[attr-defined]
        "targets": [_json_ready(asdict(item)) for item in calibration.targets],  # type: ignore[attr-defined]
        "records": [_json_ready(asdict(item)) for item in calibration.records],  # type: ignore[attr-defined]
    }


def _calibration_config(
    value: CalibrationConfig | Mapping[str, Any] | None,
) -> CalibrationConfig:
    if value is None:
        return CalibrationConfig()
    if isinstance(value, CalibrationConfig):
        return value
    return CalibrationConfig(**dict(value))


def _safe_run_name(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "-", str(value or "").strip()).strip(".-")
    if not cleaned:
        raise ValueError("run_name must contain letters, numbers, dot, underscore, or hyphen.")
    return cleaned


def _safe_upload_name(value: str, allowed_suffixes: set[str]) -> str:
    name = Path(str(value or "")).name
    if not name or name in {".", ".."}:
        raise ValueError("uploaded file name is invalid.")
    if Path(name).suffix.lower() not in allowed_suffixes:
        raise ValueError(f"unsupported uploaded file suffix: {Path(name).suffix}")
    return name


def _read_json_object(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(payload, dict):
        raise ValueError(f"JSON file must contain an object: {path}")
    return payload


def _json_ready(value: object) -> object:
    if isinstance(value, Enum):
        return value.value
    if isinstance(value, tuple):
        return [_json_ready(item) for item in value]
    if isinstance(value, list):
        return [_json_ready(item) for item in value]
    if isinstance(value, dict):
        return {str(key): _json_ready(item) for key, item in value.items()}
    return value


__all__ = [
    "DEFAULT_EXPERIMENT_FEEDBACK_ROOT",
    "export_experiment_feedback_issues",
    "export_experiment_feedback_report",
    "list_experiment_feedback_runs",
    "list_prediction_fact_packs",
    "load_experiment_feedback_run",
    "submit_experiment_feedback_import",
]
