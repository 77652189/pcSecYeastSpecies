from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

import numpy as np
from scipy.optimize import linprog

from pcsec_pichia.core.pichia_model import PichiaModel


OptimizationSense = Literal["maximize", "minimize"]


@dataclass(frozen=True)
class FluxValue:
    reaction_id: str
    flux: float
    lower_bound: float
    upper_bound: float


@dataclass(frozen=True)
class LinearProgrammingResult:
    success: bool
    status: str
    message: str
    objective_reaction: str
    objective_value: float | None = None
    sense: OptimizationSense = "maximize"
    fluxes: dict[str, float] = field(default_factory=dict)
    key_fluxes: list[FluxValue] = field(default_factory=list)

    def to_dict(self) -> dict[str, object]:
        return {
            "success": self.success,
            "status": self.status,
            "message": self.message,
            "objective_reaction": self.objective_reaction,
            "objective_value": self.objective_value,
            "sense": self.sense,
            "key_fluxes": [item.__dict__ for item in self.key_fluxes],
        }


@dataclass(frozen=True)
class ScipyHiGHSSolver:
    time_limit_seconds: float | None = 60

    def solve(
        self,
        model: PichiaModel,
        objective_reaction: str,
        sense: OptimizationSense = "maximize",
        key_reactions: list[str] | tuple[str, ...] = (),
    ) -> LinearProgrammingResult:
        if objective_reaction not in model.reaction_index:
            raise KeyError(f"Reaction not found: {objective_reaction}")
        model_with_objective = model.set_objective(objective_reaction, coefficient=1.0)
        objective = np.array(model_with_objective.c, dtype=float)
        scipy_objective = -objective if sense == "maximize" else objective
        options = {}
        if self.time_limit_seconds is not None:
            options["time_limit"] = self.time_limit_seconds
        result = linprog(
            scipy_objective,
            A_eq=model_with_objective.s_matrix,
            b_eq=np.array(model_with_objective.b, dtype=float),
            bounds=list(zip(model_with_objective.lb.astype(float), model_with_objective.ub.astype(float))),
            method="highs",
            options=options,
        )
        objective_value = None
        fluxes: dict[str, float] = {}
        key_fluxes: list[FluxValue] = []
        if result.success and result.x is not None:
            objective_value = float(result.x[model_with_objective.reaction_index[objective_reaction]])
            fluxes = {
                reaction_id: float(result.x[index])
                for reaction_id, index in model_with_objective.reaction_index.items()
                if abs(float(result.x[index])) > 1e-12
            }
            for reaction_id in key_reactions:
                if reaction_id not in model_with_objective.reaction_index:
                    continue
                index = model_with_objective.reaction_index[reaction_id]
                key_fluxes.append(
                    FluxValue(
                        reaction_id=reaction_id,
                        flux=float(result.x[index]),
                        lower_bound=float(model_with_objective.lb[index]),
                        upper_bound=float(model_with_objective.ub[index]),
                    )
                )
        return LinearProgrammingResult(
            success=bool(result.success),
            status=str(result.status),
            message=str(result.message),
            objective_reaction=objective_reaction,
            objective_value=objective_value,
            sense=sense,
            fluxes=fluxes,
            key_fluxes=key_fluxes,
        )


def write_lp_result_summary(result: LinearProgrammingResult, output_path: Path) -> Path:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(result.to_dict(), ensure_ascii=False, indent=2), encoding="utf-8")
    return output_path
