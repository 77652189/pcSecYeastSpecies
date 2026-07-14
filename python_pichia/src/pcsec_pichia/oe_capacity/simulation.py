from __future__ import annotations

from dataclasses import replace
from typing import Any, Mapping

import numpy as np

from pcsec_pichia.oe_capacity.constraints import build_oe_capacity_constraints
from pcsec_pichia.oe_capacity.schema import (
    ConfidenceLevel,
    OECapacityComparisonResult,
    OECapacityScenarioResult,
    OECapacityPlan,
    OECapacityScreenConfig,
    OECapacityScreenRequest,
    OECapacityScreenResult,
    OECapacityScreenRow,
    OECapacityValidationError,
    OEExecutionMode,
    OEExecutionStatus,
    ParameterScenario,
    SolverSnapshot,
)
from pcsec_pichia.oe_capacity.mapping import (
    build_gene_enzyme_reaction_catalog,
    fingerprint_oe_capacity_model,
)
from pcsec_pichia.oe_capacity.parameters import (
    build_current_model_parameter_policy,
    plan_gene_level_overexpression,
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
        attempt_id="legacy_baseline",
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
            scenario_results=(),
            proxy_attempts=(),
        )
        comparison.validate()
        return comparison

    result_missing_information = plan.missing_information
    result_warnings = plan.warnings
    proxy, proxy_attempts = _run_proxy_snapshots(
        prepared_model,
        plan,
        baseline_result,
        mu=mu,
        time_limit=time_limit,
    )
    scenarios: list[SolverSnapshot] = []
    scenario_results: list[OECapacityScenarioResult] = []
    for scenario in plan.uncertainty_scenarios:
        scenario_specs = tuple(
            spec
            for spec in plan.executable_capacity_specs
            if spec.parameter_scenario is scenario
        )
        if not scenario_specs:
            continue
        scenario_model = _model_with_gene_capacity_bounds(
            model,
            scenario_specs,
            multiplier=1.0,
        )
        scenario_metabolic, scenario_secretory, scenario_combined = (
            _enzyme_data_for_scenario(
                prepared_model.metabolic,
                prepared_model.secretory,
                prepared_model.combined,
                scenario_specs,
            )
        )
        scenario_baseline_result, scenario_baseline_counts = solve_pcsec_maximize(
            scenario_model,
            objective,
            metabolic=scenario_metabolic,
            secretory=scenario_secretory,
            combined=scenario_combined,
            mu=mu,
            key_reactions=key_reactions,
            time_limit_seconds=time_limit,
        )
        scenario_baseline = _snapshot(
            execution_mode=OEExecutionMode.NOT_EXECUTABLE,
            result=scenario_baseline_result,
            counts=scenario_baseline_counts,
            mu=mu,
            protein_resource_cost=_targeted_resource_cost(
                scenario_baseline_result.fluxes,
                plan,
                scenario,
            ),
            parameter_scenario=scenario,
            attempt_id=f"{scenario.value}:capacity_baseline",
        )
        multiplier = scenario_specs[0].dose.expression_multiplier
        if multiplier is None:
            raise OECapacityValidationError(
                "gene capacity scenarios require a numeric expression multiplier."
            )
        if abs(float(multiplier) - 1.0) <= 1e-12:
            perturbed = replace(
                scenario_baseline,
                execution_mode=OEExecutionMode.GENE_CAPACITY,
                attempt_id=f"{scenario.value}:gene_capacity:1x_identity",
            )
        else:
            perturbed_model = _model_with_gene_capacity_bounds(
                model,
                scenario_specs,
                multiplier=float(multiplier),
            )
            solved, counts = solve_pcsec_maximize(
                perturbed_model,
                objective,
                metabolic=scenario_metabolic,
                secretory=scenario_secretory,
                combined=scenario_combined,
                mu=mu,
                key_reactions=key_reactions,
                time_limit_seconds=time_limit,
            )
            perturbed = _snapshot(
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
                attempt_id=f"{scenario.value}:gene_capacity:{float(multiplier):g}x",
            )
        failure_reason = _scenario_failure_reason(
            scenario,
            legacy_baseline=baseline,
            scenario_baseline=scenario_baseline,
            perturbed=perturbed,
        )
        scenario_result = OECapacityScenarioResult(
            parameter_scenario=scenario,
            baseline=scenario_baseline,
            perturbed=perturbed,
            objective_delta=_difference(
                perturbed.secretion_objective,
                scenario_baseline.secretion_objective,
            ),
            protein_resource_cost_delta=_difference(
                perturbed.protein_resource_cost,
                scenario_baseline.protein_resource_cost,
            ),
            failure_reason=failure_reason,
        )
        scenario_result.validate()
        scenario_results.append(scenario_result)
        scenarios.append(perturbed)

    nominal = next(
        (
            snapshot
            for snapshot in scenarios
            if snapshot.parameter_scenario is ParameterScenario.NOMINAL
        ),
        scenarios[0] if scenarios else None,
    )
    nominal_objective = nominal.secretion_objective if nominal else None
    proxy_objective = proxy.secretion_objective if proxy else None
    nominal_pair = next(
        (
            item
            for item in scenario_results
            if item.parameter_scenario is ParameterScenario.NOMINAL
        ),
        scenario_results[0] if scenario_results else None,
    )
    failed_scenarios = tuple(
        item.parameter_scenario.value
        for item in scenario_results
        if item.failure_reason
    )
    if failed_scenarios:
        result_missing_information = tuple(
            dict.fromkeys((*result_missing_information, "successful_capacity_scenarios"))
        )
        result_warnings = tuple(
            dict.fromkeys(
                (
                    *result_warnings,
                    "Capacity scenario failures were retained: "
                    + ", ".join(failed_scenarios),
                )
            )
        )
    comparison = OECapacityComparisonResult(
        gene_id=plan.gene_id,
        target_id=plan.target_id,
        context_id=plan.context_id,
        execution_status=plan.execution_status,
        baseline=baseline,
        proxy=proxy,
        gene_capacity_scenarios=tuple(scenarios),
        gene_capacity_vs_baseline_delta=_difference(
            nominal_pair.perturbed.secretion_objective if nominal_pair else None,
            nominal_pair.baseline.secretion_objective if nominal_pair else None,
        ),
        gene_capacity_vs_proxy_delta=_difference(
            nominal_objective,
            proxy_objective,
        ),
        protein_resource_cost_delta=_difference(
            nominal_pair.perturbed.protein_resource_cost if nominal_pair else None,
            nominal_pair.baseline.protein_resource_cost if nominal_pair else None,
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
            ("capacity_basis", "reviewed_absolute_model_flux_anchor"),
        ),
        warnings=result_warnings,
        scenario_results=tuple(scenario_results),
        proxy_attempts=proxy_attempts,
    )
    comparison.validate()
    return comparison


def run_gene_level_oe_screen(
    prepared_model: Any,
    requests: tuple[OECapacityScreenRequest, ...] | list[OECapacityScreenRequest],
    screen_config: OECapacityScreenConfig,
) -> OECapacityScreenResult:
    screen_config.validate()
    resolved_requests = tuple(requests)
    for request in resolved_requests:
        request.validate()
    request_identities = tuple(
        (request.gene_id, request.target_id, request.context_id, request.dose.dose_id)
        for request in resolved_requests
    )
    if len(set(request_identities)) != len(request_identities):
        raise OECapacityValidationError("screen requests must be unique.")

    model = getattr(prepared_model, "fixed_model", None)
    if model is None:
        raise OECapacityValidationError(
            "gene-level OE screen requires prepared_model.fixed_model."
        )
    fixed_growth_rate = _fixed_growth_rate(model)
    if abs(fixed_growth_rate - float(screen_config.growth_rate)) > 1e-12:
        raise OECapacityValidationError(
            "screen_config growth_rate does not match the prepared fixed BIOMASS bound."
        )
    catalog = getattr(prepared_model, "gene_capacity_catalog", None)
    if catalog is None:
        catalog = build_gene_enzyme_reaction_catalog(
            model,
            prepared_model.metabolic,
            prepared_model.combined,
        )
    catalog.validate()
    if catalog.model_fingerprint != fingerprint_oe_capacity_model(model):
        raise OECapacityValidationError(
            "prepared gene capacity catalog does not match the fixed model."
        )
    parameter_policy = getattr(prepared_model, "parameter_policy", None)
    if parameter_policy is None:
        parameter_policy = build_current_model_parameter_policy(
            catalog,
            prepared_model.combined,
        )
    parameter_policy = replace(
        parameter_policy,
        scenarios=screen_config.parameter_scenarios,
    )
    parameter_policy.validate()

    solver_options = dict(screen_config.solver_options)
    rows: list[OECapacityScreenRow] = []
    failures: list[OECapacityScreenRow] = []
    target_id = str(getattr(prepared_model, "target_id", ""))
    for request in resolved_requests:
        if request.target_id != target_id:
            failures.append(
                _request_failure_row(
                    request,
                    "request target does not match prepared target",
                )
            )
            continue
        try:
            plan = plan_gene_level_overexpression(
                model,
                request.gene_id,
                request.target_id,
                request.context_id,
                request.dose,
                catalog,
                parameter_policy,
            )
            plan = _plan_for_screen_mode(plan, request, screen_config)
            if plan.execution_mode is OEExecutionMode.NOT_EXECUTABLE:
                failures.append(_screen_row(plan, None))
                continue
            comparison = run_gene_level_oe_comparison(
                prepared_model,
                plan,
                solver_options,
            )
            row = _screen_row(plan, comparison)
            (rows if _comparison_succeeded(plan, comparison) else failures).append(row)
        except Exception as exc:
            failures.append(_request_failure_row(request, str(exc)))

    result_warnings = (
        ("Gene-capacity feature is disabled; only legacy reaction proxy was executed.",)
        if not screen_config.feature_enabled
        else ()
    )
    result = OECapacityScreenResult(
        model_fingerprint=catalog.model_fingerprint,
        config=screen_config,
        rows=tuple(rows),
        failures=tuple(failures),
        warnings=result_warnings,
    )
    result.validate()
    return result


def _plan_for_screen_mode(
    plan: OECapacityPlan,
    request: OECapacityScreenRequest,
    config: OECapacityScreenConfig,
) -> OECapacityPlan:
    requested_mode = request.execution_mode
    if requested_mode is OEExecutionMode.NOT_EXECUTABLE:
        return replace(
            plan,
            execution_mode=OEExecutionMode.NOT_EXECUTABLE,
            executable_capacity_specs=(),
            proxy_reaction_ids=(),
            uncertainty_scenarios=(),
            missing_information=tuple(
                dict.fromkeys((*plan.missing_information, "execution_not_requested"))
            ),
        )
    if not config.feature_enabled or requested_mode is OEExecutionMode.REACTION_PROXY:
        if not plan.proxy_reaction_ids:
            return replace(
                plan,
                execution_mode=OEExecutionMode.NOT_EXECUTABLE,
                executable_capacity_specs=(),
                uncertainty_scenarios=(),
                missing_information=tuple(
                    dict.fromkeys((*plan.missing_information, "proxy_reaction_mapping"))
                ),
            )
        return replace(
            plan,
            execution_mode=OEExecutionMode.REACTION_PROXY,
            execution_status=OEExecutionStatus.PROXY_ONLY,
            executable_capacity_specs=(),
            uncertainty_scenarios=(),
        )
    if requested_mode is OEExecutionMode.GENE_CAPACITY or not config.compare_proxy:
        if not plan.executable_capacity_specs:
            return replace(
                plan,
                execution_mode=OEExecutionMode.NOT_EXECUTABLE,
                proxy_reaction_ids=(),
                uncertainty_scenarios=(),
                missing_information=tuple(
                    dict.fromkeys((*plan.missing_information, "capacity_parameters"))
                ),
            )
        return replace(
            plan,
            execution_mode=OEExecutionMode.GENE_CAPACITY,
            proxy_reaction_ids=(),
        )
    return plan


def _comparison_succeeded(
    plan: OECapacityPlan,
    comparison: OECapacityComparisonResult,
) -> bool:
    if not comparison.baseline.success:
        return False
    proxy_attempts = comparison.proxy_attempts or (
        ((comparison.proxy,) if comparison.proxy is not None else ())
    )
    proxy_success = bool(proxy_attempts) and all(
        snapshot.success for snapshot in proxy_attempts
    )
    capacity_success = (
        bool(comparison.scenario_results)
        and all(
            item.baseline.success
            and item.perturbed.success
            and not item.failure_reason
            for item in comparison.scenario_results
        )
    ) or (
        not comparison.scenario_results
        and bool(comparison.gene_capacity_scenarios)
        and all(snapshot.success for snapshot in comparison.gene_capacity_scenarios)
    )
    if plan.execution_mode is OEExecutionMode.REACTION_PROXY:
        return proxy_success
    if plan.execution_mode is OEExecutionMode.GENE_CAPACITY:
        return capacity_success
    if plan.execution_mode is OEExecutionMode.COMPARISON:
        return capacity_success and proxy_success
    return False


def _screen_row(
    plan: OECapacityPlan,
    comparison: OECapacityComparisonResult | None,
) -> OECapacityScreenRow:
    specs = plan.executable_capacity_specs
    mappings = tuple(
        dict.fromkeys(
            (
                *(spec.mapping.mapping_id for spec in specs),
                *(mapping.mapping_id for mapping in plan.explain_only_mappings),
            )
        )
    )
    estimates = tuple(
        estimate
        for spec in specs
        for estimate in (
            spec.kcat,
            spec.molecular_weight,
            spec.baseline_enzyme_amount,
            spec.complex_stoichiometry,
        )
        if estimate is not None
    )
    parameter_sources = tuple(
        dict.fromkeys(
            f"{estimate.source_type.value}:{estimate.source_ref}"
            for estimate in estimates
        )
    )
    confidence = _lowest_confidence(
        tuple(
            (
                *(spec.mapping.mapping_confidence for spec in specs),
                *(estimate.confidence for estimate in estimates),
            )
        )
    )
    nominal = None
    if comparison is not None:
        nominal_pair = next(
            (
                item
                for item in comparison.scenario_results
                if item.parameter_scenario is ParameterScenario.NOMINAL
            ),
            comparison.scenario_results[0]
            if comparison.scenario_results
            else None,
        )
        nominal = nominal_pair.perturbed if nominal_pair is not None else None
    warnings = list(plan.warnings)
    missing = list(plan.missing_information)
    if comparison is not None:
        warnings.extend(comparison.warnings)
        missing.extend(comparison.missing_information)
        if comparison.proxy is not None:
            warnings.extend(comparison.proxy.warnings)
        warnings.extend(
            warning
            for snapshot in (*comparison.gene_capacity_scenarios, *comparison.proxy_attempts)
            for warning in snapshot.warnings
        )
    screen_status = _screen_status(plan, comparison)
    return OECapacityScreenRow(
        gene_id=plan.gene_id,
        target_id=plan.target_id,
        context_id=plan.context_id,
        execution_mode=plan.execution_mode,
        execution_status=plan.execution_status,
        dose_id=plan.requested_dose.dose_id,
        dose_mode=plan.requested_dose.dose_mode,
        expression_multiplier=plan.requested_dose.expression_multiplier,
        mapping_ids=mappings,
        parameter_sources=parameter_sources,
        parameter_confidence=confidence,
        uncertainty_scenarios=plan.uncertainty_scenarios,
        baseline_objective=(
            comparison.baseline.secretion_objective if comparison is not None else None
        ),
        proxy_objective=(
            comparison.proxy.secretion_objective
            if comparison is not None and comparison.proxy is not None
            else None
        ),
        gene_capacity_objective=(nominal.secretion_objective if nominal else None),
        gene_capacity_vs_baseline_delta=(
            comparison.gene_capacity_vs_baseline_delta
            if comparison is not None
            else None
        ),
        gene_capacity_vs_proxy_delta=(
            comparison.gene_capacity_vs_proxy_delta
            if comparison is not None
            else None
        ),
        protein_resource_cost_delta=(
            comparison.protein_resource_cost_delta
            if comparison is not None
            else None
        ),
        missing_information=tuple(dict.fromkeys(missing)),
        warnings=tuple(dict.fromkeys(warnings)),
        screen_status=screen_status,
        scenario_results=(comparison.scenario_results if comparison is not None else ()),
        proxy_attempts=(comparison.proxy_attempts if comparison is not None else ()),
        summary_source=(
            "nominal_scenario_perturbed"
            if nominal is not None
            else ("best_successful_proxy" if comparison and comparison.proxy else "none")
        ),
    )


def _request_failure_row(
    request: OECapacityScreenRequest,
    reason: str,
) -> OECapacityScreenRow:
    return OECapacityScreenRow(
        gene_id=request.gene_id,
        target_id=request.target_id,
        context_id=request.context_id,
        execution_mode=OEExecutionMode.NOT_EXECUTABLE,
        execution_status=OEExecutionStatus.UNRESOLVED,
        dose_id=request.dose.dose_id,
        dose_mode=request.dose.dose_mode,
        expression_multiplier=request.dose.expression_multiplier,
        mapping_ids=(),
        parameter_sources=(),
        parameter_confidence=None,
        uncertainty_scenarios=(),
        baseline_objective=None,
        proxy_objective=None,
        gene_capacity_objective=None,
        gene_capacity_vs_baseline_delta=None,
        gene_capacity_vs_proxy_delta=None,
        protein_resource_cost_delta=None,
        missing_information=("screen_execution",),
        warnings=(reason or "screen request failed",),
        screen_status="failed",
        summary_source="none",
    )


def _screen_status(
    plan: OECapacityPlan,
    comparison: OECapacityComparisonResult | None,
) -> str:
    if plan.execution_mode is OEExecutionMode.NOT_EXECUTABLE:
        return "not_executable"
    if comparison is None or not comparison.baseline.success:
        return "failed"
    scenario_states = tuple(
        item.baseline.success and item.perturbed.success and not item.failure_reason
        for item in comparison.scenario_results
    ) or tuple(item.success for item in comparison.gene_capacity_scenarios)
    proxy_states = tuple(item.success for item in comparison.proxy_attempts) or (
        ((comparison.proxy.success,) if comparison.proxy is not None else ())
    )
    required_states: tuple[bool, ...]
    if plan.execution_mode is OEExecutionMode.GENE_CAPACITY:
        required_states = scenario_states
    elif plan.execution_mode is OEExecutionMode.REACTION_PROXY:
        required_states = proxy_states
    else:
        required_states = (*scenario_states, *proxy_states)
    if required_states and all(required_states):
        return "completed"
    if any(required_states):
        return "partial_failure"
    return "failed"


def _lowest_confidence(
    values: tuple[ConfidenceLevel, ...],
) -> ConfidenceLevel | None:
    if not values:
        return None
    priority = {
        ConfidenceLevel.HIGH: 0,
        ConfidenceLevel.MEDIUM: 1,
        ConfidenceLevel.LOW: 2,
    }
    return max(values, key=priority.__getitem__)


def _run_proxy_snapshots(
    prepared_model: Any,
    plan: OECapacityPlan,
    baseline_result: Any,
    *,
    mu: float,
    time_limit: float,
) -> tuple[SolverSnapshot | None, tuple[SolverSnapshot, ...]]:
    multiplier = plan.requested_dose.expression_multiplier
    if not plan.proxy_reaction_ids or multiplier is None:
        return None, ()
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
    attempts = tuple(_proxy_snapshot(row) for row in rows)
    successful = tuple(snapshot for snapshot in attempts if snapshot.success)
    selected = max(
        successful or attempts,
        key=lambda snapshot: (
            snapshot.secretion_objective
            if snapshot.secretion_objective is not None
            else float("-inf")
        ),
        default=None,
    )
    if selected is None:
        return None, ()
    summary_warnings = [
        *selected.warnings,
        "Legacy reaction proxy is reported independently from gene capacity.",
        f"Selected legacy proxy reaction: {selected.attempt_id}.",
    ]
    if len(attempts) > 1:
        summary_warnings.append(
            "Multiple proxy reactions were evaluated independently; all attempts are "
            "retained and the summary selects the best successful objective."
        )
    return replace(selected, warnings=tuple(dict.fromkeys(summary_warnings))), attempts


def _proxy_snapshot(row: Mapping[str, Any]) -> SolverSnapshot:
    reaction_id = str(row.get("reaction") or "unknown_proxy")
    return SolverSnapshot(
        execution_mode=OEExecutionMode.REACTION_PROXY,
        backend="scipy_highs_reference",
        solver_status=str(row.get("status") or "unknown"),
        success=bool(row.get("success")),
        secretion_objective=_optional_float(row.get("objective_value")),
        growth_retention=None,
        max_feasible_growth_rate=None,
        protein_resource_cost=None,
        constraint_counts=tuple(
            sorted(
                (str(key), int(value))
                for key, value in dict(row.get("constraint_counts") or {}).items()
            )
        ),
        key_fluxes=(),
        message=str(row.get("message") or ""),
        warnings=(() if bool(row.get("success")) else (f"Proxy attempt failed: {reaction_id}.",)),
        attempt_id=reaction_id,
    )


def _model_with_gene_capacity_bounds(
    model: Any,
    specs: tuple[Any, ...],
    *,
    multiplier: float,
) -> Any:
    changes: dict[str, tuple[float | None, float | None]] = {}
    for spec in specs:
        formation_id = spec.mapping.formation_or_dilution_reaction_id
        baseline_amount = spec.baseline_enzyme_amount
        if baseline_amount is None:
            raise OECapacityValidationError(
                "gene capacity bounds require a reviewed baseline capacity anchor."
            )
        if baseline_amount.unit != "model_flux":
            raise OECapacityValidationError(
                "baseline capacity anchor must use canonical unit model_flux."
            )
        absolute_capacity = baseline_amount.value_for_scenario(
            spec.parameter_scenario
        )
        upper_bound = float(multiplier) * absolute_capacity
        existing = changes.get(formation_id)
        if existing is not None and not np.isclose(
            float(existing[1]), upper_bound, rtol=0.0, atol=1e-12
        ):
            raise OECapacityValidationError(
                "conflicting reviewed capacity anchors for formation handle: "
                f"{formation_id}"
            )
        changes[formation_id] = (None, upper_bound)
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
    attempt_id: str = "",
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
        attempt_id=attempt_id,
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
        formation_id = spec.mapping.formation_or_dilution_reaction_id
        if formation_id in seen or spec.molecular_weight is None:
            continue
        seen.add(formation_id)
        if formation_id not in fluxes:
            continue
        costs.append(
            abs(float(fluxes[formation_id]))
            * spec.molecular_weight.value_for_scenario(scenario)
            / 1000.0
        )
    return sum(costs) if costs else None


def _scenario_failure_reason(
    scenario: ParameterScenario,
    *,
    legacy_baseline: SolverSnapshot,
    scenario_baseline: SolverSnapshot,
    perturbed: SolverSnapshot,
) -> str:
    if not scenario_baseline.success:
        return "scenario_baseline_failed"
    if not perturbed.success:
        return "scenario_perturbation_failed"
    if scenario is ParameterScenario.NOMINAL and not _objectives_close(
        scenario_baseline.secretion_objective,
        legacy_baseline.secretion_objective,
    ):
        return "capacity_baseline_incompatible_with_legacy"
    return ""


def _objectives_close(left: float | None, right: float | None) -> bool:
    if left is None or right is None:
        return False
    return bool(np.isclose(float(left), float(right), rtol=1e-6, atol=1e-10))


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


__all__ = ["run_gene_level_oe_comparison", "run_gene_level_oe_screen"]
