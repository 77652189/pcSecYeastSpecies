from __future__ import annotations

from dataclasses import asdict
from typing import Any, Mapping

from pcsec_pichia.analysis.shadow_lp import (
    ConstraintBlock,
    ConstraintSpec,
    LPProblem,
    ShadowConstraintConfig,
    SolverBackend,
    SolverResult,
)


class FakeBackend:
    name = "fake"
    supports_duals = False
    supports_time_limit = True

    def available(self) -> bool:
        return True

    def solve(self, problem: LPProblem, options: Mapping[str, Any] | None = None) -> SolverResult:
        scale = float((options or {}).get("scale", 1.0))
        objective = sum(problem.objective.values()) * scale
        return SolverResult(
            success=True,
            status="optimal",
            objective=objective,
            fluxes={reaction_id: 0.0 for reaction_id in problem.reaction_ids},
            message="fake solve",
            timings={"solve_seconds": 0.0},
            backend_metadata={"backend": self.name},
        )


def test_constraint_spec_and_block_are_serializable_contracts() -> None:
    spec = ConstraintSpec(
        name="CM1",
        layer="metabolic_coupling",
        sense="eq",
        terms={"RXN1": 1.0, "ENZ1_formation": -2.5},
        rhs=0.0,
        source="metabolic enzymedata",
        enabled_by_default=True,
        metadata={"enzyme_id": "ENZ1"},
    )
    block = ConstraintBlock(
        layer_id="metabolic_coupling",
        constraints=(spec,),
        counts={"metabolic_coupling": 1},
        mapped_reaction_count=2,
        missing_mapping_count=0,
    )

    payload = asdict(block)

    assert payload["layer_id"] == "metabolic_coupling"
    assert payload["constraints"][0]["sense"] == "eq"
    assert payload["constraints"][0]["terms"]["ENZ1_formation"] == -2.5
    assert payload["mapped_reaction_count"] == 2


def test_shadow_constraint_config_keeps_reference_defaults_explicit() -> None:
    config = ShadowConstraintConfig()

    assert config.growth_rate == 0.10
    assert config.total_protein_content == 0.37
    assert config.unmodeled_er_protein_fraction == 0.040
    assert config.mitochondrial_protein_fraction == 0.05
    assert config.enable_ribosome_translation is False
    assert config.enable_misfolding is False


def test_fake_backend_satisfies_solver_backend_protocol_and_result_contract() -> None:
    spec = ConstraintSpec(
        name="protein_mass_total",
        layer="protein_mass",
        sense="eq",
        terms={"dilute_dummy": 40.0},
        rhs=0.004,
        source="shadow config and combined enzymedata",
        enabled_by_default=True,
    )
    problem = LPProblem(
        reaction_ids=("BIOMASS", "target_exchange"),
        bounds=((0.10, 0.10), (0.0, 1000.0)),
        objective={"target_exchange": 1.25},
        stoichiometric_matrix=None,
        rhs=None,
        constraint_blocks=(ConstraintBlock("protein_mass", (spec,)),),
        metadata={"target_id": "hLF"},
    )
    backend: SolverBackend = FakeBackend()

    result = backend.solve(problem, options={"scale": 2.0})

    assert isinstance(backend, SolverBackend)
    assert backend.available() is True
    assert result.success is True
    assert result.status == "optimal"
    assert result.objective == 2.5
    assert result.fluxes == {"BIOMASS": 0.0, "target_exchange": 0.0}
    assert result.backend_metadata["backend"] == "fake"
