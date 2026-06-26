from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class AlignmentSummary:
    target_id: str
    python_result_status: str = "draft"
    matlab_alignment_status: str = "pending"
    baseline_available: bool = False
    matlab_success: bool | None = None
    rows_diff: int | None = None
    cols_diff: int | None = None
    objective_relative_diff: float | None = None
    constraint_diff_status: str = "pending"
    diagnostic_message: str = ""
    artifact_path: Path | None = None
    success: bool = False
    baseline_source: str = "matlab"
    rows_python: int | None = None
    rows_matlab: int | None = None
    cols_python: int | None = None
    cols_matlab: int | None = None
    objective_python: float | None = None
    objective_matlab: float | None = None
    diagnostics: dict[str, Any] = field(default_factory=dict)
    compatibility_exceptions: tuple[dict[str, Any], ...] = ()
    report_path: Path | None = None

    @property
    def relative_objective_diff(self) -> float | None:
        return self.objective_relative_diff


ALIGNMENT_STATUSES = {
    "pending",
    "baseline_missing",
    "matlab_failed",
    "python_draft",
    "aligned",
    "not_aligned",
    "aligned_except_known_matlab_compatibility_differences",
}


KNOWN_OPN_MATLAB_COMPATIBILITY_EXCEPTIONS: tuple[dict[str, Any], ...] = (
    {
        "id": "corrected_medium_exchange_bounds",
        "category": "bound_difference",
        "count": 9,
        "description": "Corrected Python medium opens minimal exchange bounds that are closed in the existing MATLAB artifact.",
    },
    {
        "id": "misfolding_dilution_bounds",
        "category": "bound_difference",
        "count": 1418,
        "description": "Existing MATLAB artifact fixes most dilution_misfolding variables to zero; corrected Python keeps them open.",
    },
    {
        "id": "target_secretory_coupling_missing_in_matlab_artifact",
        "category": "row_coefficient_difference",
        "count": 18,
        "description": "MATLAB artifact omits OPN target secretory reaction coupling coefficients that corrected Python includes.",
    },
    {
        "id": "ribosome_optional_row_order_and_proteins_term",
        "category": "row_coefficient_difference",
        "count": 2,
        "description": "Ribosome optional rows differ by row body order and a Python PROTEINS term in the alignment probe.",
    },
)


KNOWN_HLF_PROJECT_710_MATLAB_COMPATIBILITY_EXCEPTIONS: tuple[dict[str, Any], ...] = (
    {
        "id": "corrected_medium_exchange_bounds",
        "category": "bound_difference",
        "count": 9,
        "description": "Corrected Python medium opens minimal exchange bounds that are closed in the MATLAB hLF_PROJECT_710 artifact.",
    },
    {
        "id": "misfolding_dilution_bounds",
        "category": "bound_difference",
        "count": 1418,
        "description": "MATLAB hLF_PROJECT_710 artifact fixes most dilution_misfolding variables to zero; corrected Python keeps them open.",
    },
    {
        "id": "ribosome_optional_row_mapping",
        "category": "row_coefficient_difference",
        "count": 2,
        "description": "Ribosome optional rows differ by row body/order in the alignment probe; probe-only replacement reduces row coefficient differences to zero.",
    },
)


def opn_known_matlab_compatibility_exceptions() -> tuple[dict[str, Any], ...]:
    return tuple(dict(item) for item in KNOWN_OPN_MATLAB_COMPATIBILITY_EXCEPTIONS)


def hlf_project_710_known_matlab_compatibility_exceptions() -> tuple[dict[str, Any], ...]:
    return tuple(dict(item) for item in KNOWN_HLF_PROJECT_710_MATLAB_COMPATIBILITY_EXCEPTIONS)


def load_matlab_alignment_artifact(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"artifact_path": str(path), "baseline_available": False, "targets": []}
    payload = json.loads(path.read_text(encoding="utf-8"))
    if isinstance(payload, list):
        return {"artifact_path": str(path), "baseline_available": True, "targets": payload}
    if isinstance(payload, dict):
        if "targets" not in payload and payload.get("target_id") is not None:
            return {"artifact_path": str(path), "baseline_available": True, **payload, "targets": [payload]}
        return {"artifact_path": str(path), "baseline_available": True, **payload}
    return {"artifact_path": str(path), "baseline_available": False, "targets": []}


def build_alignment_summary(
    target_id: str,
    python_result_status: str = "draft",
    artifact_path: Path | None = None,
    artifact: dict[str, Any] | None = None,
    rows_python: int | None = None,
    cols_python: int | None = None,
    objective_python: float | None = None,
    constraint_diff_status: str | None = None,
    compatibility_exceptions: tuple[dict[str, Any], ...] | list[dict[str, Any]] | None = None,
) -> AlignmentSummary:
    loaded = artifact or (load_matlab_alignment_artifact(artifact_path) if artifact_path else None)
    baseline_available = bool(loaded and loaded.get("baseline_available"))
    record = _find_target_record(loaded, target_id) if baseline_available else None
    matlab_success = bool(record.get("success")) if record else None

    rows_matlab = _matlab_constraint_row_count(record)
    cols_matlab = _int_or_none((record or {}).get("lp_stats", {}).get("bound_lines"))
    objective_matlab = _float_or_none((record or {}).get("production_ratio_fixed"))
    rows_diff = rows_python - rows_matlab if rows_python is not None and rows_matlab is not None else None
    cols_diff = cols_python - cols_matlab if cols_python is not None and cols_matlab is not None else None
    objective_relative_diff = _relative_diff(objective_python, objective_matlab)
    normalized_exceptions = tuple(dict(item) for item in (compatibility_exceptions or ()))
    inferred_constraint_status = constraint_diff_status or _constraint_diff_status(rows_diff, cols_diff, record)
    diagnostic_message = _diagnostic_message(record, baseline_available, inferred_constraint_status)
    matlab_alignment_status = classify_alignment_status(
        python_result_status=python_result_status,
        baseline_available=baseline_available,
        matlab_success=matlab_success,
        objective_relative_diff=objective_relative_diff,
        constraint_diff_status=inferred_constraint_status,
        has_compatibility_exceptions=bool(normalized_exceptions),
    )
    return AlignmentSummary(
        target_id=target_id,
        python_result_status=python_result_status,
        matlab_alignment_status=matlab_alignment_status,
        baseline_available=baseline_available,
        matlab_success=matlab_success,
        rows_diff=rows_diff,
        cols_diff=cols_diff,
        objective_relative_diff=objective_relative_diff,
        constraint_diff_status=inferred_constraint_status,
        diagnostic_message=diagnostic_message,
        artifact_path=artifact_path or (Path(str(loaded["artifact_path"])) if loaded and loaded.get("artifact_path") else None),
        success=matlab_alignment_status == "aligned",
        baseline_source="matlab",
        rows_python=rows_python,
        rows_matlab=rows_matlab,
        cols_python=cols_python,
        cols_matlab=cols_matlab,
        objective_python=objective_python,
        objective_matlab=objective_matlab,
        diagnostics=dict(record or {}),
        compatibility_exceptions=normalized_exceptions,
    )


def classify_alignment_status(
    python_result_status: str,
    baseline_available: bool,
    matlab_success: bool | None,
    objective_relative_diff: float | None = None,
    constraint_diff_status: str = "pending",
    objective_tolerance: float = 0.01,
    has_compatibility_exceptions: bool = False,
) -> str:
    if not baseline_available:
        return "baseline_missing"
    if matlab_success is False:
        return "matlab_failed"
    if (
        has_compatibility_exceptions
        and matlab_success is True
        and constraint_diff_status == "known_matlab_compatibility_differences"
    ):
        return "aligned_except_known_matlab_compatibility_differences"
    if python_result_status == "corrected_condition":
        return "pending"
    if python_result_status == "draft":
        return "python_draft"
    if matlab_success is True:
        objective_ok = objective_relative_diff is not None and objective_relative_diff <= objective_tolerance
        constraints_ok = constraint_diff_status == "matched"
        return "aligned" if objective_ok and constraints_ok else "not_aligned"
    return "pending"


def summarize_alignment(summary: AlignmentSummary) -> dict[str, Any]:
    return {
        "target_id": summary.target_id,
        "python_result_status": summary.python_result_status,
        "matlab_alignment_status": summary.matlab_alignment_status,
        "baseline_available": summary.baseline_available,
        "matlab_success": summary.matlab_success,
        "rows_diff": summary.rows_diff,
        "cols_diff": summary.cols_diff,
        "objective_relative_diff": summary.objective_relative_diff,
        "constraint_diff_status": summary.constraint_diff_status,
        "diagnostic_message": summary.diagnostic_message,
        "artifact_path": str(summary.artifact_path) if summary.artifact_path else None,
        "rows_python": summary.rows_python,
        "rows_matlab": summary.rows_matlab,
        "cols_python": summary.cols_python,
        "cols_matlab": summary.cols_matlab,
        "objective_python": summary.objective_python,
        "objective_matlab": summary.objective_matlab,
        "compatibility_exceptions": list(summary.compatibility_exceptions),
        "is_fully_aligned": summary.matlab_alignment_status == "aligned",
        "is_aligned_except_known_matlab_compatibility_differences": (
            summary.matlab_alignment_status == "aligned_except_known_matlab_compatibility_differences"
        ),
    }


def _find_target_record(artifact: dict[str, Any] | None, target_id: str) -> dict[str, Any] | None:
    for record in (artifact or {}).get("targets", []):
        if record.get("target_id") == target_id:
            return record
    return None


def _int_or_none(value: Any) -> int | None:
    try:
        return int(value) if value is not None else None
    except (TypeError, ValueError):
        return None


def _float_or_none(value: Any) -> float | None:
    try:
        return float(value) if value is not None else None
    except (TypeError, ValueError):
        return None


def _relative_diff(left: float | None, right: float | None) -> float | None:
    if left is None or right is None or right == 0:
        return None
    return abs(left - right) / abs(right)


def _constraint_diff_status(rows_diff: int | None, cols_diff: int | None, record: dict[str, Any] | None) -> str:
    if not record:
        return "pending"
    if not record.get("success"):
        return "not_available"
    if rows_diff is None or cols_diff is None:
        return "row_level_diff_missing"
    return "matched" if rows_diff == 0 and cols_diff == 0 else "shape_diff"


def _matlab_constraint_row_count(record: dict[str, Any] | None) -> int | None:
    if not record:
        return None
    label_count = _lp_constraint_label_count(record.get("lp_file"))
    if label_count is not None:
        return label_count
    return _int_or_none(record.get("lp_stats", {}).get("constraint_lines"))


def _lp_constraint_label_count(path_value: Any) -> int | None:
    if not path_value:
        return None
    path = Path(str(path_value))
    if not path.exists():
        return None
    try:
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError:
        return None
    try:
        subject_index = next(index for index, line in enumerate(lines) if line.strip() == "Subject To")
        bounds_index = next(index for index, line in enumerate(lines) if line.strip() == "Bounds")
    except StopIteration:
        return None
    return sum(1 for line in lines[subject_index + 1 : bounds_index] if re.match(r"\s*[^:\s]+\s*:", line))


def _diagnostic_message(record: dict[str, Any] | None, baseline_available: bool, constraint_diff_status: str) -> str:
    if not baseline_available:
        return "MATLAB alignment artifact is missing."
    if not record:
        return "Target is missing from MATLAB alignment artifact."
    if not record.get("success"):
        stack = record.get("error_stack") or []
        stack_text = f" Stack: {'; '.join(str(item) for item in stack)}" if stack else ""
        return f"MATLAB failed: {record.get('error_identifier') or 'unknown'} {record.get('error_message') or ''}.{stack_text}"
    if constraint_diff_status == "row_level_diff_missing":
        return "MATLAB baseline exists, but row-level LP diff has not been generated."
    if constraint_diff_status == "shape_diff":
        return "MATLAB baseline shape differs from Python result."
    return "MATLAB baseline artifact was loaded."


__all__ = [
    "ALIGNMENT_STATUSES",
    "AlignmentSummary",
    "KNOWN_HLF_PROJECT_710_MATLAB_COMPATIBILITY_EXCEPTIONS",
    "KNOWN_OPN_MATLAB_COMPATIBILITY_EXCEPTIONS",
    "build_alignment_summary",
    "classify_alignment_status",
    "hlf_project_710_known_matlab_compatibility_exceptions",
    "load_matlab_alignment_artifact",
    "opn_known_matlab_compatibility_exceptions",
    "summarize_alignment",
]
