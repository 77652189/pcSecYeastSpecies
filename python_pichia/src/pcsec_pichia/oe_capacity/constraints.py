from __future__ import annotations

from typing import Any

from pcsec_pichia.oe_capacity.schema import (
    CapacityConstraintChange,
    ConstraintChangeKind,
    OECapacityConstraintBundle,
    OECapacityPlan,
    OECapacityValidationError,
)


def build_oe_capacity_constraints(
    prepared_model: Any,
    plan: OECapacityPlan,
) -> OECapacityConstraintBundle:
    plan.validate()
    model = getattr(prepared_model, "fixed_model", None)
    reaction_index = getattr(model, "reaction_index", None)
    combined = getattr(prepared_model, "combined", None)
    combined_enzymes = {
        str(enzyme_id) for enzyme_id in getattr(combined, "enzymes", ())
    }
    if reaction_index is None:
        raise OECapacityValidationError(
            "capacity constraints require prepared_model.fixed_model."
        )
    if not plan.executable_capacity_specs:
        raise OECapacityValidationError(
            "gene-level capacity constraints require executable_capacity_specs."
        )
    model_fingerprints = {
        spec.mapping.model_fingerprint for spec in plan.executable_capacity_specs
    }
    if len(model_fingerprints) != 1:
        raise OECapacityValidationError(
            "all capacity specs must use one model_fingerprint."
        )
    changes: list[CapacityConstraintChange] = []
    seen_capacity_handles: dict[
        tuple[object, str], tuple[float, float, float, float]
    ] = {}
    for spec in plan.executable_capacity_specs:
        scenario = spec.parameter_scenario
        mapping = spec.mapping
        if mapping.reaction_id not in reaction_index:
            raise OECapacityValidationError(
                f"capacity mapping reaction is missing from prepared model: {mapping.reaction_id}"
            )
        if mapping.formation_or_dilution_reaction_id not in reaction_index:
            raise OECapacityValidationError(
                "capacity mapping formation/dilution reaction is missing from "
                f"prepared model: {mapping.formation_or_dilution_reaction_id}"
            )
        if mapping.enzyme_id not in combined_enzymes:
            raise OECapacityValidationError(
                f"capacity mapping enzyme is missing from combined data: {mapping.enzyme_id}"
            )
        multiplier = spec.dose.expression_multiplier
        if multiplier is None:
            raise OECapacityValidationError(
                "executable capacity constraints require numeric expression_multiplier."
            )
        baseline = spec.baseline_enzyme_amount
        kcat = spec.kcat
        molecular_weight = spec.molecular_weight
        if baseline is None or kcat is None or molecular_weight is None:
            raise OECapacityValidationError(
                "executable capacity constraints require baseline, kcat, and molecular weight."
            )
        expected_units = (
            ("kcat", kcat.unit, "1/h"),
            ("molecular_weight", molecular_weight.unit, "g/mol"),
            ("baseline_enzyme_amount", baseline.unit, "model_flux"),
        )
        for parameter_name, actual_unit, expected_unit in expected_units:
            if actual_unit != expected_unit:
                raise OECapacityValidationError(
                    f"{parameter_name} must use canonical unit {expected_unit}; "
                    f"received {actual_unit}."
                )
        formation_id = mapping.formation_or_dilution_reaction_id
        formation_index = reaction_index[formation_id]
        original_upper_bound = float(model.ub[formation_index])
        baseline_capacity = baseline.value_for_scenario(scenario)
        handle_key = (scenario, formation_id)
        handle_signature = (
            float(multiplier),
            float(baseline_capacity),
            float(kcat.value_for_scenario(scenario)),
            float(molecular_weight.value_for_scenario(scenario)),
        )
        existing_signature = seen_capacity_handles.get(handle_key)
        if existing_signature is not None:
            if existing_signature != handle_signature:
                raise OECapacityValidationError(
                    "conflicting capacity specs for formation handle: "
                    f"{formation_id} ({scenario.value})"
                )
            continue
        seen_capacity_handles[handle_key] = handle_signature
        metadata = (
            ("mapping_id", mapping.mapping_id),
            ("gene_id", mapping.gene_id),
            ("dose_id", spec.dose.dose_id),
        )
        changes.extend(
            (
                CapacityConstraintChange(
                    change_id=f"{mapping.mapping_id}:{scenario.value}:formation_bound",
                    scenario=scenario,
                    change_kind=ConstraintChangeKind.FORMATION_DILUTION_BOUND,
                    constraint_block="gene_capacity_formation_bound",
                    variable_id=formation_id,
                    reaction_id=mapping.reaction_id,
                    old_value=original_upper_bound,
                    new_value=float(multiplier) * baseline_capacity,
                    unit="model_flux",
                    source_ref=baseline.source_ref,
                    resource_cost_mode=spec.resource_cost_mode,
                    metadata=metadata,
                ),
                CapacityConstraintChange(
                    change_id=f"{mapping.mapping_id}:{scenario.value}:kcat",
                    scenario=scenario,
                    change_kind=ConstraintChangeKind.ENZYME_CAPACITY_COEFFICIENT,
                    constraint_block="metabolic_coupling",
                    variable_id=mapping.enzyme_variable_id,
                    reaction_id=mapping.reaction_id,
                    old_value=kcat.nominal_value,
                    new_value=kcat.value_for_scenario(scenario),
                    unit=kcat.unit,
                    source_ref=kcat.source_ref,
                    resource_cost_mode=spec.resource_cost_mode,
                    metadata=metadata,
                ),
                CapacityConstraintChange(
                    change_id=f"{mapping.mapping_id}:{scenario.value}:molecular_weight",
                    scenario=scenario,
                    change_kind=ConstraintChangeKind.PROTEIN_RESOURCE_COEFFICIENT,
                    constraint_block="protein_mass",
                    variable_id=formation_id,
                    reaction_id=mapping.reaction_id,
                    old_value=molecular_weight.nominal_value,
                    new_value=molecular_weight.value_for_scenario(scenario),
                    unit=molecular_weight.unit,
                    source_ref=molecular_weight.source_ref,
                    resource_cost_mode=spec.resource_cost_mode,
                    metadata=metadata,
                ),
            )
        )
    bundle = OECapacityConstraintBundle(
        model_fingerprint=model_fingerprints.pop(),
        plan=plan,
        changes=tuple(changes),
    )
    bundle.validate()
    return bundle


__all__ = ["build_oe_capacity_constraints"]
