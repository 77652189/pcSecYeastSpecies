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


def resolve_oe_product_plan(
    plan: Any,
    *,
    requested_mode: Any,
    feature_enabled: bool,
    compare_proxy: bool,
) -> Any:
    from pcsec_pichia.oe_capacity.product_tiering import (
        resolve_oe_product_plan as _resolve,
    )

    return _resolve(
        plan,
        requested_mode=requested_mode,
        feature_enabled=feature_enabled,
        compare_proxy=compare_proxy,
    )


def summarize_oe_product_candidate(plan: Any) -> dict[str, object]:
    from pcsec_pichia.oe_capacity.product_tiering import (
        summarize_oe_product_candidate as _summarize,
    )

    return _summarize(plan)


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
    from pcsec_pichia.oe_capacity.simulation import run_gene_level_oe_screen as _run

    return _run(prepared_model, requests, screen_config)


def write_oe_capacity_outputs(
    result: Any,
    output_dir: Any,
    *,
    run_identity: Mapping[str, Any] | None = None,
    capacity_asset: Mapping[str, Any] | None = None,
) -> Any:
    from pcsec_pichia.oe_capacity.reports import write_oe_capacity_outputs as _write

    return _write(
        result,
        output_dir,
        run_identity=run_identity,
        capacity_asset=capacity_asset,
    )


def verify_oe_capacity_run(
    run_dir: Any,
    expected_identity: Mapping[str, object] | None = None,
) -> Any:
    from pcsec_pichia.oe_capacity.acceptance import verify_oe_capacity_run as _verify

    return _verify(run_dir, expected_identity)


def run_phase2_oe_capacity_acceptance(repo_root: Any, output_root: Any) -> Any:
    from pcsec_pichia.oe_capacity.acceptance import (
        run_phase2_oe_capacity_acceptance as _run,
    )

    return _run(repo_root, output_root)


__all__ = [
    "OECapacityPhaseError",
    "build_gene_capacity_specs",
    "build_gene_enzyme_reaction_catalog",
    "build_oe_capacity_constraints",
    "build_oe_dose_spec",
    "plan_gene_level_overexpression",
    "resolve_oe_product_plan",
    "run_gene_level_oe_comparison",
    "run_gene_level_oe_screen",
    "run_phase2_oe_capacity_acceptance",
    "summarize_oe_product_candidate",
    "validate_gene_capacity_catalog",
    "verify_oe_capacity_run",
    "write_oe_capacity_outputs",
]
