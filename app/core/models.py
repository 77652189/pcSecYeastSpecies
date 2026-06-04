from __future__ import annotations

from pathlib import Path
from typing import Any, Literal

from pydantic import BaseModel, Field


SpeciesCode = Literal["SCE", "PPA", "KMX", "Unknown"]


class DatasetInfo(BaseModel):
    id: str
    name: str
    path: Path
    category: str
    suffix: str
    species: SpeciesCode = "Unknown"
    size_bytes: int
    modified_at: str


class LoadedDataset(BaseModel):
    info: DatasetInfo
    kind: Literal["table", "mat"]
    tables: dict[str, list[dict[str, Any]]] = Field(default_factory=dict)
    variable_summary: list[dict[str, Any]] = Field(default_factory=list)


class HealthItem(BaseModel):
    name: str
    status: Literal["ok", "warning", "missing", "error"]
    detail: str = ""


class HealthReport(BaseModel):
    items: list[HealthItem]
    preflight_output: str = ""


class SoplexSummary(BaseModel):
    optimal: bool
    objective_value: str | None = None
    status_line: str | None = None


class SimulationResult(BaseModel):
    success: bool
    mu: float
    message: str
    lp_file: Path | None = None
    output_file: Path | None = None
    objective_value: str | None = None
    command_output: str = ""
