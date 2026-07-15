from __future__ import annotations

from dataclasses import replace

from pcsec_pichia.oe_capacity.schema import (
    AbsoluteCapacityAvailability,
    OECapacityPlan,
    OECapacityValidationError,
    OECalibrationStatus,
    OEExecutionMode,
    OEExecutionStatus,
    OEProductMode,
    OEProductState,
)


def requested_product_mode(
    *,
    product_mode: OEProductMode | None,
    execution_mode: OEExecutionMode,
) -> OEProductMode:
    if product_mode is not None:
        return product_mode
    if execution_mode is OEExecutionMode.REACTION_PROXY:
        return OEProductMode.REACTION_PROXY
    if execution_mode is OEExecutionMode.GENE_CAPACITY:
        return OEProductMode.ABSOLUTE_CAPACITY
    if execution_mode is OEExecutionMode.NOT_EXECUTABLE:
        return OEProductMode.NOT_EXECUTABLE
    # Historical `comparison` requests become the explicitly labelled relative
    # decision product. They never imply that an absolute anchor exists.
    return OEProductMode.RELATIVE_UNCALIBRATED


def resolve_oe_product_plan(
    plan: OECapacityPlan,
    *,
    requested_mode: OEProductMode,
    feature_enabled: bool,
    compare_proxy: bool,
) -> OECapacityPlan:
    plan.validate()
    if not isinstance(requested_mode, OEProductMode):
        raise OECapacityValidationError("requested_mode must be an OEProductMode.")

    if not feature_enabled or requested_mode is OEProductMode.REACTION_PROXY:
        return _reaction_proxy_plan(plan)
    if requested_mode is OEProductMode.RELATIVE_UNCALIBRATED:
        return _relative_plan(plan, compare_proxy=compare_proxy)
    if requested_mode is OEProductMode.ABSOLUTE_CAPACITY:
        return _absolute_plan(plan, compare_proxy=compare_proxy)
    return _not_executable_plan(plan, "execution_not_requested")


def summarize_oe_product_candidate(plan: OECapacityPlan) -> dict[str, object]:
    plan.validate()
    return {
        "product_mode": plan.product_mode.value,
        "product_state": plan.product_state.value,
        "execution_mode": plan.execution_mode.value,
        "execution_status": plan.execution_status.value,
        "absolute_capacity_availability": plan.absolute_capacity_availability.value,
        "calibration_status": plan.calibration_status.value,
        "absolute_solver_allowed": plan.absolute_solver_allowed,
        "model_fingerprint": plan.model_fingerprint,
        "mapping_ids": [mapping.mapping_id for mapping in plan.structural_mappings],
        "relative_scenarios": [
            {
                "scenario": spec.parameter_scenario.value,
                "relative_capacity_factor": spec.relative_capacity_factor,
                "parameter_sources": list(spec.parameter_sources),
                "warnings": list(spec.warnings),
                "limitations": list(spec.limitations),
            }
            for spec in plan.relative_scenario_specs
        ],
        "missing_information": list(plan.missing_information),
        "warnings": list(plan.warnings),
        "limitations": list(plan.limitations),
    }


def _reaction_proxy_plan(plan: OECapacityPlan) -> OECapacityPlan:
    if not plan.proxy_reaction_ids:
        return _not_executable_plan(plan, "proxy_reaction_mapping")
    resolved = replace(
        plan,
        execution_mode=OEExecutionMode.REACTION_PROXY,
        execution_status=OEExecutionStatus.PROXY_ONLY,
        executable_capacity_specs=(),
        relative_scenario_specs=(),
        constraint_changes=(),
        uncertainty_scenarios=(),
        product_mode=OEProductMode.REACTION_PROXY,
        product_state=OEProductState.REACTION_PROXY,
        absolute_capacity_availability=_absolute_availability(plan),
        calibration_status=OECalibrationStatus.PROXY_ONLY,
        absolute_solver_allowed=False,
        limitations=tuple(
            dict.fromkeys(
                (
                    *plan.limitations,
                    "reaction_proxy_is_not_gene_level_capacity",
                    "no_absolute_capacity_claim",
                )
            )
        ),
    )
    resolved.validate()
    return resolved


def _relative_plan(plan: OECapacityPlan, *, compare_proxy: bool) -> OECapacityPlan:
    relative_specs = (
        plan.relative_scenario_specs or plan.available_relative_scenario_specs
    )
    if not relative_specs:
        return _not_executable_plan(plan, "relative_parameter_mapping")
    scenarios = tuple(
        dict.fromkeys(spec.parameter_scenario for spec in relative_specs)
    )
    resolved = replace(
        plan,
        execution_mode=OEExecutionMode.RELATIVE_GENE_CAPACITY,
        execution_status=(
            OEExecutionStatus.PARTIAL_MAPPING
            if plan.execution_status is OEExecutionStatus.PARTIAL_MAPPING
            else OEExecutionStatus.GENE_LEVEL_EXECUTABLE
        ),
        executable_capacity_specs=(),
        relative_scenario_specs=relative_specs,
        constraint_changes=(),
        proxy_reaction_ids=(plan.proxy_reaction_ids if compare_proxy else ()),
        uncertainty_scenarios=scenarios,
        product_mode=OEProductMode.RELATIVE_UNCALIBRATED,
        product_state=OEProductState.RELATIVE_UNCALIBRATED,
        absolute_capacity_availability=_absolute_availability(plan),
        calibration_status=OECalibrationStatus.RELATIVE_UNCALIBRATED,
        absolute_solver_allowed=False,
        warnings=tuple(
            dict.fromkeys(
                (
                    *plan.warnings,
                    "Relative OE applies a dimensionless factor to the targeted "
                    "current-model enzyme-capacity coupling; it does not represent "
                    "a measured kcat, true expression fold-change, or absolute capacity.",
                )
            )
        ),
        limitations=tuple(
            dict.fromkeys(
                (
                    *plan.limitations,
                    "relative_scenarios_are_uncalibrated",
                    *(
                        ("optional_proxy_comparison_is_directional_evidence_only",)
                        if compare_proxy and plan.proxy_reaction_ids
                        else ()
                    ),
                    "no_absolute_nominal_capacity",
                )
            )
        ),
    )
    resolved.validate()
    return resolved


def _absolute_plan(plan: OECapacityPlan, *, compare_proxy: bool) -> OECapacityPlan:
    if not plan.executable_capacity_specs or not plan.absolute_solver_allowed:
        resolved = replace(
            plan,
            execution_mode=OEExecutionMode.NOT_EXECUTABLE,
            execution_status=(
                OEExecutionStatus.PARTIAL_MAPPING
                if plan.structural_mappings
                else OEExecutionStatus.UNRESOLVED
            ),
            executable_capacity_specs=(),
            relative_scenario_specs=(),
            proxy_reaction_ids=(),
            constraint_changes=(),
            uncertainty_scenarios=(),
            product_mode=OEProductMode.ABSOLUTE_CAPACITY,
            product_state=OEProductState.ABSOLUTE_UNAVAILABLE,
            absolute_capacity_availability=(
                AbsoluteCapacityAvailability.UNAVAILABLE_MISSING_REVIEWED_ANCHOR
            ),
            calibration_status=OECalibrationStatus.UNAVAILABLE,
            absolute_solver_allowed=False,
            missing_information=tuple(
                dict.fromkeys(
                    (*plan.missing_information, "reviewed_baseline_capacity")
                )
            ),
            warnings=tuple(
                dict.fromkeys(
                    (
                        *plan.warnings,
                        "Absolute capacity was requested but no compatible reviewed "
                        "baseline anchor is available; no solver was called.",
                    )
                )
            ),
            limitations=tuple(
                dict.fromkeys((*plan.limitations, "absolute_capacity_unavailable"))
            ),
        )
        resolved.validate()
        return resolved
    mode = (
        OEExecutionMode.COMPARISON
        if compare_proxy and plan.proxy_reaction_ids
        else OEExecutionMode.GENE_CAPACITY
    )
    resolved = replace(
        plan,
        execution_mode=mode,
        relative_scenario_specs=(),
        proxy_reaction_ids=(plan.proxy_reaction_ids if mode is OEExecutionMode.COMPARISON else ()),
        product_mode=OEProductMode.ABSOLUTE_CAPACITY,
        product_state=OEProductState.ABSOLUTE_AVAILABLE,
        absolute_capacity_availability=AbsoluteCapacityAvailability.AVAILABLE_REVIEWED,
        calibration_status=OECalibrationStatus.REVIEWED_ABSOLUTE,
        absolute_solver_allowed=True,
    )
    resolved.validate()
    return resolved


def _not_executable_plan(plan: OECapacityPlan, reason: str) -> OECapacityPlan:
    resolved = replace(
        plan,
        execution_mode=OEExecutionMode.NOT_EXECUTABLE,
        executable_capacity_specs=(),
        relative_scenario_specs=(),
        proxy_reaction_ids=(),
        constraint_changes=(),
        uncertainty_scenarios=(),
        product_mode=OEProductMode.NOT_EXECUTABLE,
        product_state=OEProductState.NOT_EXECUTABLE,
        absolute_capacity_availability=_absolute_availability(plan),
        calibration_status=OECalibrationStatus.NOT_APPLICABLE,
        absolute_solver_allowed=False,
        missing_information=tuple(
            dict.fromkeys((*plan.missing_information, reason))
        ),
        limitations=tuple(dict.fromkeys((*plan.limitations, "not_executable"))),
    )
    resolved.validate()
    return resolved


def _absolute_availability(plan: OECapacityPlan) -> AbsoluteCapacityAvailability:
    return (
        AbsoluteCapacityAvailability.AVAILABLE_REVIEWED
        if plan.executable_capacity_specs and plan.absolute_solver_allowed
        else AbsoluteCapacityAvailability.UNAVAILABLE_MISSING_REVIEWED_ANCHOR
    )


__all__ = [
    "requested_product_mode",
    "resolve_oe_product_plan",
    "summarize_oe_product_candidate",
]
