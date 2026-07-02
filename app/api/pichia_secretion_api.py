"""Experimental FastAPI facade for pcSecPichia draft runs.

This module intentionally stays at the app-service boundary. It maps HTTP
payloads to ``SecretionRunRequest`` and delegates execution/status handling to
``app.services.pichia_secretion_service``.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Literal

from fastapi import FastAPI
from pydantic import BaseModel, Field

from app.services.pichia_secretion_service import (
    SecretionRunRequest,
    poll_background_simulation,
    status_path_for_background_task,
    submit_background_simulation,
)


TargetSource = Literal["builtin", "custom_json", "custom_sequence"]
SequenceRole = Literal["preprotein", "mature_secreted", "custom_user_sequence", "unknown"]
NormalizationMode = Literal["as_provided", "remove_terminal_stop", "none"]
TerminalStopPolicy = Literal["reject", "strip", "allow_for_record_only"]


class PichiaSecretionRunRequest(BaseModel):
    target_source: TargetSource = "builtin"
    target_id: str = "OPN_ALPHA_FULL_PROJECT"
    target_name: str = ""
    sequence: str = ""
    leader_sequence: str = ""
    signal_peptide_sequence: str = ""
    sequence_role: SequenceRole = "unknown"
    normalization_mode: NormalizationMode = "as_provided"
    contains_signal_peptide: bool | None = None
    contains_leader: bool | None = None
    terminal_stop_policy: TerminalStopPolicy = "allow_for_record_only"
    disulfide_sites: int = 0
    n_glycosylation_sites: int = 0
    o_glycosylation_sites: int = 0
    custom_json_path: str | None = None
    enable_ribosome_translation_constraint: bool = False
    enable_misfolding_constraint: bool = False
    ko_gene_ids: list[str] = Field(default_factory=list)
    oe_gene_ids: list[str] = Field(default_factory=list)
    ko_reaction_ids: list[str] = Field(default_factory=list)
    oe_reaction_ids: list[str] = Field(default_factory=list)
    screen_candidate_limit: int = 20
    growth_points: list[float] = Field(default_factory=lambda: [0.10])
    mu: float = 0.10
    media_type: int = 4
    carbon_source_id: str = "glucose"
    output_dir: str | None = None

    def to_service_request(self) -> SecretionRunRequest:
        return SecretionRunRequest(
            target_source=self.target_source,
            target_id=self.target_id,
            target_name=self.target_name,
            sequence=self.sequence,
            leader_sequence=self.leader_sequence,
            signal_peptide_sequence=self.signal_peptide_sequence,
            sequence_role=self.sequence_role,
            normalization_mode=self.normalization_mode,
            contains_signal_peptide=self.contains_signal_peptide,
            contains_leader=self.contains_leader,
            terminal_stop_policy=self.terminal_stop_policy,
            disulfide_sites=self.disulfide_sites,
            n_glycosylation_sites=self.n_glycosylation_sites,
            o_glycosylation_sites=self.o_glycosylation_sites,
            custom_json_path=Path(self.custom_json_path) if self.custom_json_path else None,
            enable_ribosome_translation_constraint=self.enable_ribosome_translation_constraint,
            enable_misfolding_constraint=self.enable_misfolding_constraint,
            ko_gene_ids=tuple(_compact_ids(self.ko_gene_ids)),
            oe_gene_ids=tuple(_compact_ids(self.oe_gene_ids)),
            ko_reaction_ids=tuple(_compact_ids(self.ko_reaction_ids)),
            oe_reaction_ids=tuple(_compact_ids(self.oe_reaction_ids)),
            screen_candidate_limit=self.screen_candidate_limit,
            growth_points=tuple(self.growth_points or [0.10]),
            mu=self.mu,
            media_type=self.media_type,
            carbon_source_id=self.carbon_source_id,
            output_dir=Path(self.output_dir) if self.output_dir else None,
        )


def create_app() -> FastAPI:
    app = FastAPI(
        title="pcSecPichia Python Draft API (Experimental)",
        version="0.1.0",
        description=(
            "Experimental minimal API facade for Python corrected-condition "
            "draft Pichia secretion simulations. Core modeling remains inside "
            "python_pichia and is reached only through the app service facade."
        ),
    )

    @app.post("/pichia/secretion/runs")
    def create_secretion_run(request: PichiaSecretionRunRequest) -> dict[str, Any]:
        service_request = request.to_service_request()
        task_id, status_path = submit_background_simulation(service_request)
        return {
            "task_id": task_id,
            "status": "pending",
            "message": "任务已提交，结果仍是 Python draft；默认使用 corrected medium。",
            "status_path": str(status_path),
            "status_url": app.url_path_for("get_secretion_run_status", task_id=task_id),
            "result_url": app.url_path_for("get_secretion_run_result", task_id=task_id),
        }

    @app.get("/pichia/secretion/runs/{task_id}/status", name="get_secretion_run_status")
    def get_secretion_run_status(task_id: str) -> dict[str, Any]:
        status, message, result = poll_background_simulation(_status_path_for_task(task_id))
        return {
            "task_id": task_id,
            "status": status,
            "message": message,
            "has_result": result is not None,
        }

    @app.get("/pichia/secretion/runs/{task_id}/result", name="get_secretion_run_result")
    def get_secretion_run_result(task_id: str) -> dict[str, Any]:
        status, message, result = poll_background_simulation(_status_path_for_task(task_id))
        return {
            "task_id": task_id,
            "status": status,
            "message": message,
            "result": result,
        }

    return app


def _status_path_for_task(task_id: str) -> Path:
    return status_path_for_background_task(task_id)


def _compact_ids(values: list[str]) -> list[str]:
    ids: list[str] = []
    for value in values:
        text = str(value or "").strip()
        if text:
            ids.append(text)
    return ids


app = create_app()
