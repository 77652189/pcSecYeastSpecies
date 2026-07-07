from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping, Protocol, runtime_checkable

from pcsec_pichia.analysis.shadow_lp.constraint_spec import LPProblem


@dataclass(frozen=True)
class SolverResult:
    """Backend-neutral result for a shadow LP solve."""

    success: bool
    status: str
    objective: float | None
    fluxes: Mapping[str, float]
    message: str
    timings: Mapping[str, float] = field(default_factory=dict)
    duals: Mapping[str, Any] | None = None
    backend_metadata: Mapping[str, Any] = field(default_factory=dict)


@runtime_checkable
class SolverBackend(Protocol):
    """Protocol implemented by replaceable shadow LP solver backends."""

    name: str
    supports_duals: bool
    supports_time_limit: bool

    def available(self) -> bool:
        """Return whether this backend can solve in the current environment."""

    def solve(self, problem: LPProblem, options: Mapping[str, Any] | None = None) -> SolverResult:
        """Solve a backend-neutral LP problem."""
