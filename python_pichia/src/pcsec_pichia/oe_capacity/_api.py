from __future__ import annotations

from typing import Any, Mapping, Sequence

from pcsec_pichia.oe_capacity.schema import OECapacityError


class OECapacityPhaseError(OECapacityError):
    def __init__(self, api_name: str, required_round: int) -> None:
        self.api_name = api_name
        self.required_round = required_round
        super().__init__(
            f"{api_name} is defined by the Phase 2 contract and becomes executable in "
            f"Round {required_round}."
        )


def _phase_gate(api_name: str, required_round: int) -> None:
    raise OECapacityPhaseError(api_name, required_round)


def build_gene_enzyme_reaction_catalog(
    model: Any,
    metabolic: Any,
    combined: Any,
    external_evidence: Any = None,
) -> Any:
    from pcsec_pichia.oe_capacity.mapping import (
        build_gene_enzyme_reaction_catalog as _build,
    )

    return _build(model, metabolic, combined, external_evidence)


def validate_gene_capacity_catalog(catalog: Any) -> Any:
    from pcsec_pichia.oe_capacity.mapping import (
        validate_gene_capacity_catalog as _validate,
    )

    return _validate(catalog)


def build_oe_dose_spec(
    payload: Mapping[str, Any],
    dose_mapping: Mapping[str, Any] | None = None,
) -> Any:
    from pcsec_pichia.oe_capacity.parameters import build_oe_dose_spec as _build

    return _build(payload, dose_mapping)


def build_gene_capacity_specs(
    gene_id: str,
    catalog: Any,
    dose: Any,
    parameter_policy: Any,
) -> Any:
    from pcsec_pichia.oe_capacity.parameters import build_gene_capacity_specs as _build

    return _build(gene_id, catalog, dose, parameter_policy)


def plan_gene_level_overexpression(
    model: Any,
    gene_id: str,
    target_id: str,
    context_id: str,
    dose: Any,
    catalog: Any,
    parameter_policy: Any,
) -> Any:
    from pcsec_pichia.oe_capacity.parameters import (
        plan_gene_level_overexpression as _plan,
    )

    return _plan(
        model,
        gene_id,
        target_id,
        context_id,
        dose,
        catalog,
        parameter_policy,
    )


def build_oe_capacity_constraints(prepared_model: Any, plan: Any) -> Any:
    from pcsec_pichia.oe_capacity.constraints import (
        build_oe_capacity_constraints as _build,
    )

    return _build(prepared_model, plan)


def run_gene_level_oe_comparison(
    prepared_model: Any,
    plan: Any,
    solver_options: Mapping[str, Any] | None = None,
) -> Any:
    from pcsec_pichia.oe_capacity.simulation import (
        run_gene_level_oe_comparison as _run,
    )

    return _run(prepared_model, plan, solver_options)


def run_gene_level_oe_screen(
    prepared_model: Any,
    requests: Sequence[Any],
    screen_config: Any,
) -> Any:
    _phase_gate("run_gene_level_oe_screen", 4)


def write_oe_capacity_outputs(result: Any, output_dir: Any) -> Any:
    _phase_gate("write_oe_capacity_outputs", 4)


__all__ = [
    "OECapacityPhaseError",
    "build_gene_capacity_specs",
    "build_gene_enzyme_reaction_catalog",
    "build_oe_capacity_constraints",
    "build_oe_dose_spec",
    "plan_gene_level_overexpression",
    "run_gene_level_oe_comparison",
    "run_gene_level_oe_screen",
    "validate_gene_capacity_catalog",
    "write_oe_capacity_outputs",
]
