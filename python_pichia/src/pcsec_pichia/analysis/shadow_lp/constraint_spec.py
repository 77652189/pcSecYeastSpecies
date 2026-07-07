from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Literal, Mapping


ConstraintSense = Literal["eq", "le", "ge"]


@dataclass(frozen=True)
class ConstraintSpec:
    """A symbolic LP constraint before reaction IDs are assembled into columns."""

    name: str
    layer: str
    sense: ConstraintSense
    terms: Mapping[str, float]
    rhs: float
    source: str
    enabled_by_default: bool
    metadata: Mapping[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class ConstraintBlock:
    """A layer-level group of constraints plus mapping diagnostics."""

    layer_id: str
    constraints: tuple[ConstraintSpec, ...]
    counts: Mapping[str, int] = field(default_factory=dict)
    mapped_reaction_count: int = 0
    missing_mapping_count: int = 0
    warnings: tuple[str, ...] = ()
    metadata: Mapping[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class ShadowConstraintConfig:
    """Explicit defaults for pcSec shadow resource constraints."""

    growth_rate: float = 0.10
    total_protein_content: float = 0.37
    unmodeled_er_protein_fraction: float = 0.040
    mitochondrial_protein_fraction: float = 0.05
    enable_ribosome_translation: bool = False
    enable_misfolding: bool = False
    metadata: Mapping[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class LPProblem:
    """Backend-neutral LP payload for shadow pcSec solves."""

    reaction_ids: tuple[str, ...]
    bounds: tuple[tuple[float | None, float | None], ...]
    objective: Mapping[str, float]
    stoichiometric_matrix: Any
    rhs: Any
    constraint_blocks: tuple[ConstraintBlock, ...]
    metadata: Mapping[str, Any] = field(default_factory=dict)
