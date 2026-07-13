from __future__ import annotations

from dataclasses import replace
from typing import Any, Mapping

import numpy as np

from pcsec_pichia.oe_capacity.constraints import build_oe_capacity_constraints
from pcsec_pichia.oe_capacity.schema import (
    OECapacityComparisonResult,
    OECapacityPlan,
    OECapacityValidationError,
    OEExecutionMode,
    ParameterScenario,
    SolverSnapshot,
)
from pcsec_pichia.probe import run_pcsec_oe_screen, solve_pcsec_maximize


def run_gene_level_oe_comparison(
    prepared_model: Any,
    plan: OECapacityPlan,
    solver_options: Mapping[str, Any] | None = None,
) -> OECapacityComparisonResult:
    plan.validate()
    if str(getattr(prepared_model, "target_id", "")) != plan.target_id:
        raise OECapacityValidationError(
            "prepared target_id does not match OECapacityPlan target_id."
        )
    bundle = (
        build_oe_capacity_constraints(prepared_model, plan)
        if plan.executable_capacity_specs
        else None
    )
    model = prepared_model.fixed_model
    objective = str(prepared_model.exchange_reaction_id)
    mu = _fixed_growth_rate(model)
    options = dict(solver_options or {})
    time_limit = float(options.get("time_limit_seconds", 600.0))
    formation_ids = tuple(
        dict.fromkeys(
            spec.mapping.formation_or_dilution_reaction_id
            for spec in plan.executable_capacity_specs
        )
    )
    key_reactions = tuple(
        dict.fromkeys(("BIOMASS", objective, *formation_ids, *plan.proxy_reaction_ids))
    )
    baseline_result, baseline_counts = solve_pcsec_maximize(
        model,
        objective,
        metabolic=prepared_model.metabolic,
        secretory=prepared_model.secretory,
        combined=prepared_model.combined,
        mu=mu,
        key_reactions=key_reactions,
        time_limit_seconds=time_limit,
    )
    baseline_cost = _targeted_resource_cost(
        baseline_result.fluxes,
        plan,
        ParameterScenario.NOMINAL,
    )
    baseline = _snapshot(
        execution_mode=OEExecutionMode.NOT_EXECUTABLE,
        result=baseline_result,
        counts=baseline_counts,
        mu=mu,
        protein_resource_cost=baseline_cost,
    )
    if not baseline_result.success:
        comparison = OECapacityComparisonResult(
            gene_id=plan.gene_id,
            target_id=plan.target_id,
            context_id=plan.context_id,
            execution_status=plan.execution_status,
            baseline=baseline,
            proxy=None,
            gene_capacity_scenarios=(),
            gene_capacity_vs_baseline_delta=None,
            gene_capacity_vs_proxy_delta=None,
            protein_resource_cost_delta=None,
            skipped_reason="baseline_solve_failed",
            missing_information=plan.missing_information,
            warnings=plan.warnings,
        )
        comparison.validate()
        return comparison

    baseline_formation_fluxes = tuple(
        (formation_id, float(baseline_result.fluxes[formation_id]))
        for formation_id in formation_ids
        if formation_id in baseline_result.fluxes
    )
    zero_baseline_formations = tuple(
        formation_id
        for formation_id, flux in baseline_formation_fluxes
        if abs(flux) <= 1e-12
    )
    result_missing_information = plan.missing_information
    result_warnings = plan.warnings
    if zero_baseline_formations:
        result_missing_information = tuple(
            dict.fromkeys(
                (*result_missing_information, "nonzero_baseline_formation_flux")
            )
        )
        result_warnings = tuple(
            dict.fromkeys(
                (
                    *result_warnings,
                    "Relative gene-capacity multipliers cannot create nonzero capacity "
                    "for formation handles with zero baseline flux: "
                    + ", ".join(zero_baseline_formations),
                )
            )
        )

    proxy = _run_proxy_snapshot(
        prepared_model,
        plan,
        baseline_result,
        mu=mu,
        time_limit=time_limit,
    )
    scenarios: list[SolverSnapshot] = []
    for scenario in plan.uncertainty_scenarios:
        scenario_specs = tuple(
            spec
            for spec in plan.executable_capacity_specs
            if spec.parameter_scenario is scenario
        )
        if not scenario_specs:
            continue
        perturbed_model = _model_with_gene_capacity_bounds(
            model,
            baseline_result.fluxes,
            scenario_specs,
        )
        perturbed_metabolic, perturbed_secretory, perturbed_combined = (
            _enzyme_data_for_scenario(
                prepared_model.metabolic,
                prepared_model.secretory,
                prepared_model.combined,
                scenario_specs,
            )
        )
        solved, counts = solve_pcsec_maximize(
            perturbed_model,
            objective,
            metabolic=perturbed_metabolic,
            secretory=perturbed_secretory,
            combined=perturbed_combined,
            mu=mu,
            key_reactions=key_reactions,
            time_limit_seconds=time_limit,
        )
        scenarios.append(
            _snapshot(
                execution_mode=OEExecutionMode.GENE_CAPACITY,
                result=solved,
                counts=counts,
                mu=mu,
                protein_resource_cost=_targeted_resource_cost(
                    solved.fluxes,
                    plan,
                    scenario,
                ),
                parameter_scenario=scenario,
            )
        )

    nominal = next(
        (
            snapshot
            for snapshot in scenarios
            if snapshot.parameter_scenario is ParameterScenario.NOMINAL
        ),
        scenarios[0] if scenarios else None,
    )
    baseline_objective = baseline.secretion_objective
    nominal_objective = nominal.secretion_objective if nominal else None
    proxy_objective = proxy.secretion_objective if proxy else None
    comparison = OECapacityComparisonResult(
        gene_id=plan.gene_id,
        target_id=plan.target_id,
        context_id=plan.context_id,
        execution_status=plan.execution_status,
        baseline=baseline,
        proxy=proxy,
        gene_capacity_scenarios=tuple(scenarios),
        gene_capacity_vs_baseline_delta=_difference(
            nominal_objective,
            baseline_objective,
        ),
        gene_capacity_vs_proxy_delta=_difference(
            nominal_objective,
            proxy_objective,
        ),
        protein_resource_cost_delta=_difference(
            nominal.protein_resource_cost if nominal else None,
            baseline.protein_resource_cost,
        ),
        skipped_reason=(
            ""
            if scenarios
            else (
                "gene_capacity_not_executable"
                if not plan.executable_capacity_specs
                else "gene_capacity_scenarios_unavailable"
            )
        ),
        missing_information=result_missing_information,
        traceability=(
            (
                "constraint_change_count",
                str(len(bundle.changes) if bundle is not None else 0),
            ),
            ("solver_backend", "scipy_highs_reference"),
            *tuple(
                (f"baseline_formation_flux:{formation_id}", _format_flux(flux))
                for formation_id, flux in baseline_formation_fluxes
            ),
        ),
        warnings=result_warnings,
    )
    comparison.validate()
    return comparison


def _run_proxy_snapshot(
    prepared_model: Any,
    plan: OECapacityPlan,
    baseline_result: Any,
    *,
    mu: float,
    time_limit: float,
) -> SolverSnapshot | None:
    multiplier = plan.requested_dose.expression_multiplier
    if not plan.proxy_reaction_ids or multiplier is None:
        return None
    rows = run_pcsec_oe_screen(
        prepared_model.fixed_model,
        baseline_result,
        list(plan.proxy_reaction_ids),
        str(prepared_model.exchange_reaction_id),
        metabolic=prepared_model.metabolic,
        secretory=prepared_model.secretory,
        combined=prepared_model.combined,
        mu=mu,
        factor=float(multiplier),
        time_limit_seconds=time_limit,
    )
    successful = tuple(row for row in rows if bool(row.get("success")))
    selected = max(
        successful or tuple(rows),
        key=lambda row: (
            float(row.get("objective_value"))
            if row.get("objective_value") is not None
            else float("-inf")
        ),
        default=None,
    )
    if selected is None:
        return None
    selected_reaction = str(selected.get("reaction") or "")
    warnings = [
        "Legacy reaction proxy is reported independently from gene capacity.",
        f"Selected legacy proxy reaction: {selected_reaction}.",
    ]
    if len(plan.proxy_reaction_ids) > 1:
        warnings.append(
            "Multiple proxy reactions were evaluated independently; this snapshot "
            "reports the best single-reaction result, not a joint proxy perturbation."
        )
    return SolverSnapshot(
        execution_mode=OEExecutionMode.REACTION_PROXY,
        backend="scipy_highs_reference",
        solver_status=str(selected.get("status") or "unknown"),
        success=bool(selected.get("success")),
        secretion_objective=_optional_float(selected.get("objective_value")),
        growth_retention=None,
        max_feasible_growth_rate=None,
        protein_resource_cost=None,
        constraint_counts=tuple(
            sorted(
                (str(key), int(value))
                for key, value in dict(selected.get("constraint_counts") or {}).items()
            )
        ),
        key_fluxes=(),
        message=str(selected.get("message") or ""),
        warnings=tuple(warnings),
    )


def _model_with_gene_capacity_bounds(
    model: Any,
    baseline_fluxes: Mapping[str, float],
    specs: tuple[Any, ...],
) -> Any:
    changes: dict[str, tuple[float | None, float | None]] = {}
    for spec in specs:
        formation_id = spec.mapping.formation_or_dilution_reaction_id
        if formation_id not in baseline_fluxes:
            raise OECapacityValidationError(
                f"baseline solve did not return formation flux: {formation_id}"
            )
        baseline_flux = max(0.0, float(baseline_fluxes[formation_id]))
        baseline_amount = spec.baseline_enzyme_amount
        multiplier = spec.dose.expression_multiplier
        if baseline_amount is None or multiplier is None:
            raise OECapacityValidationError(
                "gene capacity bounds require baseline amount and numeric dose."
            )
        relative_capacity = (
            float(multiplier)
            * baseline_amount.value_for_scenario(spec.parameter_scenario)
            / baseline_amount.nominal_value
        )
        changes[formation_id] = (None, baseline_flux * relative_capacity)
    return model.with_bounds(changes)


def _enzyme_data_for_scenario(
    metabolic: Any,
    secretory: Any,
    combined: Any,
    specs: tuple[Any, ...],
) -> tuple[Any, Any, Any]:
    metabolic_enzymes = tuple(str(item) for item in metabolic.enzymes)
    metabolic_kcat = np.array(metabolic.kcat, dtype=float, copy=True)
    metabolic_index = {
        enzyme_id: position for position, enzyme_id in enumerate(metabolic_enzymes)
    }
    secretory_complexes = tuple(str(item) for item in secretory.complexes)
    secretory_kcat = np.array(secretory.kcat, dtype=float, copy=True)
    enzymes = tuple(str(item) for item in combined.enzymes)
    kcat = np.array(combined.kcat, dtype=float, copy=True)
    enzyme_mw = np.array(combined.enzyme_mw, dtype=float, copy=True)
    index = {enzyme_id: position for position, enzyme_id in enumerate(enzymes)}
    for spec in specs:
        position = index.get(spec.mapping.enzyme_id)
        if position is None:
            raise OECapacityValidationError(
                f"combined enzyme data missing {spec.mapping.enzyme_id}."
            )
        if spec.kcat is not None:
            scenario_kcat = spec.kcat.value_for_scenario(spec.parameter_scenario)
            kcat[position] = scenario_kcat
            coupling_updated = False
            metabolic_position = metabolic_index.get(spec.mapping.enzyme_id)
            if metabolic_position is not None:
                metabolic_kcat[metabolic_position] = scenario_kcat
                coupling_updated = True
            for secretory_position, complex_id in enumerate(secretory_complexes):
                if complex_id == spec.mapping.enzyme_id:
                    secretory_kcat[secretory_position] = scenario_kcat
                    coupling_updated = True
            if not coupling_updated:
                raise OECapacityValidationError(
                    "no active metabolic or secretory coupling found for "
                    f"{spec.mapping.enzyme_id}."
                )
        if spec.molecular_weight is not None:
            enzyme_mw[position] = spec.molecular_weight.value_for_scenario(
                spec.parameter_scenario
            )
    return (
        replace(metabolic, kcat=metabolic_kcat),
        replace(secretory, kcat=secretory_kcat),
        replace(combined, kcat=kcat, enzyme_mw=enzyme_mw),
    )


def _snapshot(
    *,
    execution_mode: OEExecutionMode,
    result: Any,
    counts: Mapping[str, int],
    mu: float,
    protein_resource_cost: float | None,
    parameter_scenario: ParameterScenario | None = None,
) -> SolverSnapshot:
    biomass = result.fluxes.get("BIOMASS") if result.success else None
    return SolverSnapshot(
        execution_mode=execution_mode,
        backend="scipy_highs_reference",
        solver_status=str(result.status),
        success=bool(result.success),
        secretion_objective=_optional_float(result.objective_value),
        growth_retention=(
            float(biomass) / mu if biomass is not None and mu > 0 else None
        ),
        max_feasible_growth_rate=None,
        protein_resource_cost=protein_resource_cost,
        parameter_scenario=parameter_scenario,
        constraint_counts=tuple(sorted((str(key), int(value)) for key, value in counts.items())),
        key_fluxes=tuple(sorted((str(key), float(value)) for key, value in result.fluxes.items())),
        message=str(result.message),
    )


def _targeted_resource_cost(
    fluxes: Mapping[str, float],
    plan: OECapacityPlan,
    scenario: ParameterScenario,
) -> float | None:
    costs: list[float] = []
    seen: set[str] = set()
    for spec in plan.executable_capacity_specs:
        if spec.parameter_scenario is not scenario:
            continue
        mapping_id = spec.mapping.mapping_id
        if mapping_id in seen or spec.molecular_weight is None:
            continue
        seen.add(mapping_id)
        formation_id = spec.mapping.formation_or_dilution_reaction_id
        if formation_id not in fluxes:
            continue
        costs.append(
            abs(float(fluxes[formation_id]))
            * spec.molecular_weight.value_for_scenario(scenario)
            / 1000.0
        )
    return sum(costs) if costs else None


def _fixed_growth_rate(model: Any) -> float:
    index = model.reaction_index.get("BIOMASS")
    if index is None:
        raise OECapacityValidationError("prepared model is missing BIOMASS.")
    lower = float(model.lb[index])
    upper = float(model.ub[index])
    if lower <= 0 or abs(lower - upper) > 1e-12:
        raise OECapacityValidationError(
            "prepared model must have a positive fixed BIOMASS bound."
        )
    return lower


def _difference(left: float | None, right: float | None) -> float | None:
    if left is None or right is None:
        return None
    return float(left) - float(right)


def _optional_float(value: Any) -> float | None:
    return None if value is None else float(value)


def _format_flux(value: float) -> str:
    return "0" if abs(value) <= 1e-12 else f"{value:.17g}"


__all__ = ["run_gene_level_oe_comparison"]
