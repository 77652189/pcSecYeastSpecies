from __future__ import annotations

from pathlib import Path

from app import ensure_python_pichia_on_path

ensure_python_pichia_on_path()

from pcsec_pichia.core.paths import ProjectPaths

from app.services.pichia_background_tasks import (
    poll_background_simulation,
    status_path_for_background_task,
    submit_background_simulation,
)
from app.services.pichia_secretion_schema import (
    BuiltinTargetTemplate,
    NormalizationMode,
    SecretionRunRequest,
    SecretionRunResponse,
    SequenceRole,
    TargetSource,
    TerminalStopPolicy,
)
from app.services.pichia_request_mapping_service import (
    engine_target_id as _engine_target_id,
    request_warnings as _request_warnings,
    resolve_output_dir as _resolve_output_dir,
    sequence_contract_for_engine as _sequence_contract_for_engine,
    target_input_payload as _target_input_payload,
)
from app.services.pichia_secretion_runner import run_pichia_pipeline_draft


def discover_project_paths(start: Path | None = None) -> ProjectPaths:
    return ProjectPaths.discover(start or Path(__file__))


def run_pichia_secretion_draft(request: SecretionRunRequest, paths: ProjectPaths | None = None) -> SecretionRunResponse:
    resolved_paths = paths or discover_project_paths()
    ensure_python_pichia_on_path()

    warnings = _request_warnings(request)
    output_dir = _resolve_output_dir(request, resolved_paths)
    engine_target_id = _engine_target_id(request)
    return run_pichia_pipeline_draft(
        request,
        output_dir=output_dir,
        warnings=warnings,
        engine_target_id=engine_target_id,
        target_input=_target_input_payload(request),
        sequence_contract=_sequence_contract_for_engine(request),
    )


__all__ = [
    "BuiltinTargetTemplate",
    "NormalizationMode",
    "SecretionRunRequest",
    "SecretionRunResponse",
    "SequenceRole",
    "TargetSource",
    "TerminalStopPolicy",
    "discover_project_paths",
    "poll_background_simulation",
    "run_pichia_secretion_draft",
    "status_path_for_background_task",
    "submit_background_simulation",
]
