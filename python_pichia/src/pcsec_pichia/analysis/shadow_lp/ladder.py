from __future__ import annotations

import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Mapping

from pcsec_pichia.analysis.shadow_lp.backends import ScipyHighsBackend, SolverBackend
from pcsec_pichia.analysis.shadow_lp.constraint_builders import build_shadow_constraint_blocks
from pcsec_pichia.analysis.shadow_lp.constraint_spec import ConstraintBlock, ShadowConstraintConfig
from pcsec_pichia.analysis.shadow_lp.lp_problem import assemble_lp_problem, lp_problem_from_model
from pcsec_pichia.analysis.shadow_lp.model_adapter import ShadowTargetPreparation, prepare_builtin_shadow_target


FORMAL_SHADOW_LADDER_ORDER: tuple[str, ...] = (
    "target_extension",
    "base_cobrapy_fba",
    "fixed_growth",
    "metabolic_coupling",
    "secretory_coupling",
    "protein_mass",
    "proteasome",
    "ribosome_assembly",
    "ribosome_translation",
    "misfolding",
    "mitochondrial",
)

_LAYER_BLOCK_IDS: Mapping[str, tuple[str, ...]] = {
    "base_cobrapy_fba": (),
    "fixed_growth": (),
    "metabolic_coupling": ("metabolic_coupling",),
    "secretory_coupling": ("metabolic_coupling", "secretory_coupling"),
    "protein_mass": ("metabolic_coupling", "secretory_coupling", "protein_mass"),
    "proteasome": ("metabolic_coupling", "secretory_coupling", "protein_mass", "proteasome"),
    "ribosome_assembly": (
        "metabolic_coupling",
        "secretory_coupling",
        "protein_mass",
        "proteasome",
        "ribosome_assembly",
    ),
    "mitochondrial": (
        "metabolic_coupling",
        "secretory_coupling",
        "protein_mass",
        "proteasome",
        "ribosome_assembly",
        "mitochondrial",
    ),
}


@dataclass(frozen=True)
class ShadowLadderLayerResult:
    target_id: str
    layer_id: str
    success: bool
    status: str
    message: str
    objective: float | None
    key_fluxes: Mapping[str, float | None]
    variable_count: int
    constraint_count: int
    eq_constraint_count: int
    ub_constraint_count: int
    backend_metadata: Mapping[str, Any]
    timings: Mapping[str, float]
    enabled_layers: tuple[str, ...]
    skipped_layers: Mapping[str, str]
    warnings: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class ShadowLadderResult:
    target_id: str
    exchange_reaction_id: str
    backend_name: str
    layers: tuple[ShadowLadderLayerResult, ...]
    final_layer_id: str = "mitochondrial"
    warnings: tuple[str, ...] = ()
    reference_validation: Mapping[str, Any] | None = None

    @property
    def final_layer(self) -> ShadowLadderLayerResult:
        for layer in reversed(self.layers):
            if layer.layer_id == self.final_layer_id:
                return layer
        raise KeyError(f"Final layer not found: {self.final_layer_id}")

    def with_reference_validation(self, validation: Mapping[str, Any]) -> "ShadowLadderResult":
        return ShadowLadderResult(
            target_id=self.target_id,
            exchange_reaction_id=self.exchange_reaction_id,
            backend_name=self.backend_name,
            layers=self.layers,
            final_layer_id=self.final_layer_id,
            warnings=self.warnings,
            reference_validation=validation,
        )

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def run_shadow_ladder(
    target_id: str,
    root: Path | None = None,
    config: ShadowConstraintConfig | None = None,
    backend: SolverBackend | None = None,
    solver_options: Mapping[str, Any] | None = None,
) -> ShadowLadderResult:
    """Run the formal isolated shadow LP ladder for one built-in target."""

    prep = prepare_builtin_shadow_target(target_id, root=root, config=config)
    return run_shadow_ladder_for_prepared_target(
        prep,
        config=config,
        backend=backend,
        solver_options=solver_options,
    )


def run_shadow_ladder_for_prepared_target(
    prep: ShadowTargetPreparation,
    config: ShadowConstraintConfig | None = None,
    backend: SolverBackend | None = None,
    solver_options: Mapping[str, Any] | None = None,
) -> ShadowLadderResult:
    resolved_backend = backend or ScipyHighsBackend()
    resolved_options = dict(solver_options or {"time_limit": 600.0, "presolve": True})
    blocks = build_shadow_constraint_blocks(prep, config=config)
    block_by_id = {block.layer_id: block for block in blocks}
    layers: list[ShadowLadderLayerResult] = [
        _target_extension_layer(prep, resolved_backend.name),
    ]
    for layer_id in FORMAL_SHADOW_LADDER_ORDER[1:]:
        if layer_id in {"ribosome_translation", "misfolding"}:
            layers.append(_skipped_layer(prep, layer_id, block_by_id[layer_id], resolved_backend.name))
            continue
        selected_blocks = tuple(block_by_id[block_id] for block_id in _LAYER_BLOCK_IDS[layer_id])
        model = prep.model if layer_id == "base_cobrapy_fba" else prep.fixed_model
        layers.append(_solve_layer(prep, model, layer_id, selected_blocks, resolved_backend, resolved_options))

    warnings = tuple(
        warning
        for block in blocks
        for warning in block.warnings
        if block.layer_id not in {"ribosome_translation", "misfolding"}
    )
    return ShadowLadderResult(
        target_id=prep.target_id,
        exchange_reaction_id=prep.exchange_reaction_id,
        backend_name=resolved_backend.name,
        layers=tuple(layers),
        warnings=warnings,
    )


def _solve_layer(
    prep: ShadowTargetPreparation,
    model: Any,
    layer_id: str,
    blocks: tuple[ConstraintBlock, ...],
    backend: SolverBackend,
    solver_options: Mapping[str, Any],
) -> ShadowLadderLayerResult:
    problem = lp_problem_from_model(
        model,
        prep.exchange_reaction_id,
        blocks,
        key_reaction_ids=("BIOMASS", "Ex_glc_D", "Ex_o2", prep.exchange_reaction_id),
        metadata={
            "target_id": prep.target_id,
            "layer_id": layer_id,
            "canonical_final_layer": "mitochondrial",
        },
    )
    assembled = assemble_lp_problem(problem)
    result = backend.solve(problem, options=solver_options)
    selected_layer_ids = tuple(block.layer_id for block in blocks if block.constraints)
    warnings = tuple(warning for block in blocks for warning in block.warnings)
    return ShadowLadderLayerResult(
        target_id=prep.target_id,
        layer_id=layer_id,
        success=result.success,
        status=result.status,
        message=result.message,
        objective=result.objective,
        key_fluxes=result.key_fluxes,
        variable_count=assembled.diagnostics.variable_count,
        constraint_count=assembled.diagnostics.constraint_count,
        eq_constraint_count=assembled.diagnostics.eq_constraint_count,
        ub_constraint_count=assembled.diagnostics.ub_constraint_count,
        backend_metadata=result.backend_metadata,
        timings=result.timings,
        enabled_layers=selected_layer_ids,
        skipped_layers={},
        warnings=warnings,
    )


def _target_extension_layer(prep: ShadowTargetPreparation, backend_name: str) -> ShadowLadderLayerResult:
    started = time.perf_counter()
    constraint_count = int(prep.model.s_matrix.shape[0])
    return ShadowLadderLayerResult(
        target_id=prep.target_id,
        layer_id="target_extension",
        success=True,
        status="implemented_prerequisite",
        message="Target reactions and metabolites appended before shadow LP optimization.",
        objective=None,
        key_fluxes={
            "BIOMASS": None,
            "Ex_glc_D": None,
            "Ex_o2": None,
            prep.exchange_reaction_id: None,
        },
        variable_count=len(prep.model.rxns),
        constraint_count=constraint_count,
        eq_constraint_count=constraint_count,
        ub_constraint_count=0,
        backend_metadata={
            "backend": backend_name,
            "solver_backend": "not_run",
            "added_reaction_count": prep.added_reaction_count,
            "added_metabolite_count": prep.added_metabolite_count,
        },
        timings={"total_seconds": time.perf_counter() - started},
        enabled_layers=("target_extension",),
        skipped_layers={},
        warnings=prep.warnings,
    )


def _skipped_layer(
    prep: ShadowTargetPreparation,
    layer_id: str,
    block: ConstraintBlock,
    backend_name: str,
) -> ShadowLadderLayerResult:
    reason = block.warnings[0] if block.warnings else f"{layer_id} disabled by default."
    constraint_count = int(prep.fixed_model.s_matrix.shape[0])
    return ShadowLadderLayerResult(
        target_id=prep.target_id,
        layer_id=layer_id,
        success=False,
        status="skipped",
        message="Layer intentionally disabled to match current pcSec reference defaults.",
        objective=None,
        key_fluxes={
            "BIOMASS": None,
            "Ex_glc_D": None,
            "Ex_o2": None,
            prep.exchange_reaction_id: None,
        },
        variable_count=len(prep.fixed_model.rxns),
        constraint_count=constraint_count,
        eq_constraint_count=constraint_count,
        ub_constraint_count=0,
        backend_metadata={"backend": backend_name, "solver_backend": "not_run"},
        timings={"total_seconds": 0.0},
        enabled_layers=(),
        skipped_layers={layer_id: reason},
        warnings=(reason,),
    )
