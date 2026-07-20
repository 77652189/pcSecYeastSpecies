from __future__ import annotations

import json
import hashlib
from pathlib import Path
from typing import Any, Mapping

from pcsec_pichia.oe_capacity.schema import (
    ConfidenceLevel,
    CapacityAnchor,
    CapacityAnchorCatalog,
    AbsoluteCapacityAvailability,
    EvidenceSourceType,
    GeneCapacityCatalog,
    GeneCapacityParameterSet,
    GeneCapacitySpec,
    OEDoseMode,
    OEDoseSpec,
    OECapacityParameterConflictError,
    OECapacityPlan,
    OECapacityValidationError,
    OEExecutionStatus,
    OEExecutionMode,
    OECalibrationStatus,
    OEProductMode,
    OEProductState,
    ParameterPolicy,
    ResourceCostMode,
    ParameterEstimate,
    RelativeOEScenarioSpec,
)
from pcsec_pichia.oe_capacity.mapping import fingerprint_oe_capacity_model
from pcsec_pichia.screens.gene_interventions import plan_gene_overexpression


def build_oe_dose_spec(
    payload: Mapping[str, Any],
    dose_mapping: Mapping[str, Any] | None = None,
) -> OEDoseSpec:
    dose_id = str(payload.get("dose_id") or "").strip()
    promoter = str(payload.get("promoter") or "").strip()
    induction_mode = str(payload.get("induction_mode") or "").strip()
    copy_number = _optional_float(payload.get("copy_number"), "copy_number")
    explicit_multiplier = _optional_float(
        payload.get("expression_multiplier"),
        "expression_multiplier",
    )
    requested_mode = str(payload.get("dose_mode") or "").strip()
    mapped = _reviewed_dose_mapping(
        dose_mapping or {},
        dose_id=dose_id,
        promoter=promoter,
        copy_number=copy_number,
    )

    if requested_mode:
        try:
            mode = OEDoseMode(requested_mode)
        except ValueError as exc:
            raise OECapacityValidationError(
                f"unsupported dose_mode: {requested_mode}"
            ) from exc
    elif explicit_multiplier is not None:
        mode = OEDoseMode.EXPLICIT_MULTIPLIER
    elif mapped is not None:
        mode = OEDoseMode.PROMOTER_COPY_MAPPING
    else:
        mode = OEDoseMode.CATEGORICAL_ONLY

    mapping_source = ""
    multiplier = explicit_multiplier
    warnings: tuple[str, ...] = ()
    if mode is OEDoseMode.PROMOTER_COPY_MAPPING:
        if mapped is None:
            raise OECapacityValidationError(
                "promoter_copy_mapping requested but no reviewed dose mapping exists."
            )
        multiplier = _required_positive_float(
            mapped.get("expression_multiplier"),
            "mapped expression_multiplier",
        )
        mapping_source = str(mapped.get("mapping_source") or "").strip()
        if not mapping_source:
            raise OECapacityValidationError(
                "reviewed dose mapping requires mapping_source."
            )
    elif mode is OEDoseMode.CATEGORICAL_ONLY:
        if explicit_multiplier is not None:
            raise OECapacityValidationError(
                "categorical_only dose must not provide expression_multiplier."
            )
        multiplier = None
        warnings = (
            "Categorical promoter/copy input has no reviewed numeric dose mapping.",
        )

    spec = OEDoseSpec(
        dose_id=dose_id,
        dose_mode=mode,
        expression_multiplier=multiplier,
        promoter=promoter,
        copy_number=copy_number,
        induction_mode=induction_mode,
        mapping_source=mapping_source,
        warnings=warnings,
    )
    spec.validate()
    return spec


def build_gene_capacity_specs(
    gene_id: str,
    catalog: GeneCapacityCatalog,
    dose: OEDoseSpec,
    parameter_policy: ParameterPolicy,
) -> tuple[GeneCapacitySpec, ...]:
    catalog.validate()
    dose.validate()
    parameter_policy.validate()
    if dose.dose_mode is OEDoseMode.CATEGORICAL_ONLY:
        return ()
    specs: list[GeneCapacitySpec] = []
    for mapping in catalog.mappings:
        if mapping.gene_id != gene_id:
            continue
        if mapping.execution_status is not OEExecutionStatus.GENE_LEVEL_EXECUTABLE:
            continue
        selected = _select_parameter_set(
            mapping.mapping_id,
            mapping.gene_id,
            mapping.enzyme_id,
            parameter_policy,
        )
        if selected is None or selected.missing_information:
            continue
        for scenario in parameter_policy.scenarios:
            spec = GeneCapacitySpec(
                mapping=mapping,
                kcat=selected.kcat,
                molecular_weight=selected.molecular_weight,
                baseline_enzyme_amount=selected.baseline_enzyme_amount,
                complex_stoichiometry=selected.complex_stoichiometry,
                dose=dose,
                parameter_scenario=scenario,
                resource_cost_mode=ResourceCostMode.CURRENT_PROTEIN_POOL,
                capacity_anchor_binding=selected.capacity_anchor_binding,
                warnings=selected.warnings,
            )
            spec.validate()
            specs.append(spec)
    return tuple(specs)


def build_relative_oe_scenario_specs(
    gene_id: str,
    catalog: GeneCapacityCatalog,
    dose: OEDoseSpec,
    parameter_policy: ParameterPolicy,
) -> tuple[RelativeOEScenarioSpec, ...]:
    catalog.validate()
    dose.validate()
    parameter_policy.validate()
    if dose.dose_mode is OEDoseMode.CATEGORICAL_ONLY:
        return ()
    multiplier = dose.expression_multiplier
    if multiplier is None:
        return ()
    specs: list[RelativeOEScenarioSpec] = []
    for mapping in catalog.mappings:
        if (
            mapping.gene_id != gene_id
            or mapping.execution_status is not OEExecutionStatus.GENE_LEVEL_EXECUTABLE
        ):
            continue
        selected = _select_parameter_set(
            mapping.mapping_id,
            mapping.gene_id,
            mapping.enzyme_id,
            parameter_policy,
        )
        if selected is None or selected.kcat is None or selected.molecular_weight is None:
            continue
        sources = tuple(
            dict.fromkeys(
                (
                    f"mapping:{mapping.mapping_source.value}:{mapping.source_ref or mapping.mapping_id}",
                    f"kcat:{selected.kcat.source_type.value}:{selected.kcat.source_ref}",
                    "molecular_weight:"
                    f"{selected.molecular_weight.source_type.value}:"
                    f"{selected.molecular_weight.source_ref}",
                    (
                        f"dose:{dose.mapping_source}"
                        if dose.mapping_source
                        else f"dose:explicit_user_input:{dose.dose_id}"
                    ),
                )
            )
        )
        for scenario in parameter_policy.scenarios:
            spec = RelativeOEScenarioSpec(
                mapping=mapping,
                dose=dose,
                parameter_scenario=scenario,
                relative_capacity_factor=float(multiplier),
                kcat=selected.kcat,
                molecular_weight=selected.molecular_weight,
                parameter_sources=sources,
                warnings=tuple(
                    dict.fromkeys(
                        (
                            *mapping.warnings,
                            *selected.warnings,
                            *dose.warnings,
                            "Relative OE is uncalibrated and does not define an absolute model_flux capacity.",
                        )
                    )
                ),
                limitations=(
                    "cannot_interpret_as_absolute_capacity",
                    "cannot_interpret_as_true_expression_fold_change",
                    "cannot_predict_mg_per_litre_or_experimental_success",
                    "cannot_generalize_across_target_or_context",
                ),
            )
            spec.validate()
            specs.append(spec)
    return tuple(specs)


def build_current_model_parameter_policy(
    catalog: GeneCapacityCatalog,
    combined: Any,
    *,
    capacity_anchors: CapacityAnchorCatalog | None = None,
    target_id: str = "",
    context_id: str = "",
    relative_uncertainty: float = 0.0,
) -> ParameterPolicy:
    catalog.validate()
    if capacity_anchors is not None:
        capacity_anchors.validate()
    if (
        isinstance(relative_uncertainty, bool)
        or not isinstance(relative_uncertainty, (int, float))
        or relative_uncertainty < 0
        or relative_uncertainty >= 1
    ):
        raise OECapacityValidationError(
            "relative_uncertainty must be in [0, 1)."
        )
    parameter_sets: list[GeneCapacityParameterSet] = []
    for mapping in catalog.mappings:
        if mapping.execution_status is not OEExecutionStatus.GENE_LEVEL_EXECUTABLE:
            continue
        anchor = _matching_capacity_anchor(
            capacity_anchors,
            mapping=mapping,
            target_id=target_id,
            context_id=context_id,
        )
        warnings: list[str] = []
        kcat = _exact_parameter(
            combined,
            method_name="exact_enzyme_kcat",
            enzyme_id=mapping.enzyme_id,
            parameter_name="kcat",
            unit="1/h",
            source_ref=mapping.source_ref or "current combined enzyme data",
            relative_uncertainty=float(relative_uncertainty),
            source_version=catalog.model_fingerprint,
            warnings=warnings,
        )
        molecular_weight = _exact_parameter(
            combined,
            method_name="exact_enzyme_mw",
            enzyme_id=mapping.enzyme_id,
            parameter_name="molecular_weight",
            unit="g/mol",
            source_ref=mapping.source_ref or "current combined enzyme data",
            relative_uncertainty=float(relative_uncertainty),
            source_version=catalog.model_fingerprint,
            warnings=warnings,
        )
        baseline = anchor.as_parameter_estimate() if anchor is not None else None
        binding = (
            anchor.binding(
                asset_version=capacity_anchors.asset_version,
                asset_sha256=capacity_anchors.source_sha256,
            )
            if anchor is not None and capacity_anchors is not None
            else None
        )
        parameter_set = GeneCapacityParameterSet(
            parameter_set_id=f"current-{mapping.mapping_id}",
            mapping_id=mapping.mapping_id,
            gene_id=mapping.gene_id,
            enzyme_id=mapping.enzyme_id,
            kcat=kcat,
            molecular_weight=molecular_weight,
            baseline_enzyme_amount=baseline,
            complex_stoichiometry=None,
            capacity_anchor_binding=binding,
            warnings=tuple(warnings),
        )
        parameter_set.validate()
        parameter_sets.append(parameter_set)
    policy = ParameterPolicy(parameter_sets=tuple(parameter_sets))
    policy.validate()
    return policy


def load_capacity_anchor_catalog(
    source: str | Path | Mapping[str, Any],
) -> CapacityAnchorCatalog:
    source_ref = ""
    source_sha256 = ""
    if isinstance(source, Mapping):
        payload = source
        source_ref = str(payload.get("source_ref") or "")
        source_sha256 = str(payload.get("source_sha256") or "")
    else:
        path = Path(source)
        source_ref = str(path)
        try:
            raw = path.read_bytes()
            payload = json.loads(raw.decode("utf-8"))
            source_sha256 = hashlib.sha256(raw).hexdigest()
        except (OSError, json.JSONDecodeError) as exc:
            raise OECapacityValidationError(
                f"failed to load capacity anchor asset: {path}"
            ) from exc
        if not isinstance(payload, Mapping):
            raise OECapacityValidationError(
                "capacity anchor asset root must be an object."
            )
    raw_anchors = payload.get("anchors")
    if not isinstance(raw_anchors, list):
        raise OECapacityValidationError(
            "capacity anchor asset requires an anchors array."
        )
    anchors: list[CapacityAnchor] = []
    for item in raw_anchors:
        if not isinstance(item, Mapping):
            raise OECapacityValidationError("capacity anchor entries must be objects.")
        try:
            anchors.append(
                CapacityAnchor(
                    anchor_id=str(item.get("anchor_id") or ""),
                    target_id=str(item.get("target_id") or ""),
                    context_id=str(item.get("context_id") or ""),
                    gene_id=str(item.get("gene_id") or ""),
                    enzyme_id=str(item.get("enzyme_id") or ""),
                    formation_or_dilution_reaction_id=str(
                        item.get("formation_or_dilution_reaction_id") or ""
                    ),
                    model_fingerprint=str(item.get("model_fingerprint") or ""),
                    baseline_capacity=_required_positive_float(
                        item.get("baseline_capacity"), "baseline_capacity"
                    ),
                    unit=str(item.get("unit") or ""),
                    source_ref=str(item.get("source_ref") or ""),
                    source_version=str(item.get("source_version") or ""),
                    reviewed_by=str(item.get("reviewed_by") or ""),
                    reviewed_at=str(item.get("reviewed_at") or ""),
                )
            )
        except (TypeError, ValueError) as exc:
            raise OECapacityValidationError("invalid capacity anchor entry.") from exc
    try:
        schema_version = int(payload.get("schema_version", 0))
    except (TypeError, ValueError) as exc:
        raise OECapacityValidationError(
            "capacity anchor schema_version must be an integer."
        ) from exc
    catalog = CapacityAnchorCatalog(
        model_fingerprint=str(payload.get("model_fingerprint") or ""),
        anchors=tuple(anchors),
        schema_version=schema_version,
        source_ref=source_ref,
        asset_version=str(payload.get("asset_version") or ""),
        source_sha256=source_sha256,
    )
    catalog.validate()
    return catalog


def plan_gene_level_overexpression(
    model: Any,
    gene_id: str,
    target_id: str,
    context_id: str,
    dose: OEDoseSpec,
    catalog: GeneCapacityCatalog,
    parameter_policy: ParameterPolicy,
) -> OECapacityPlan:
    catalog.validate()
    dose.validate()
    parameter_policy.validate()
    model_fingerprint = fingerprint_oe_capacity_model(model)
    if catalog.model_fingerprint != model_fingerprint:
        raise OECapacityValidationError(
            "catalog model_fingerprint does not match the supplied model."
        )
    mappings = tuple(mapping for mapping in catalog.mappings if mapping.gene_id == gene_id)
    proxy_plan = plan_gene_overexpression(model, gene_id)
    proxy_reactions = tuple(proxy_plan.executable_reactions)
    explain_only = tuple(
        mapping
        for mapping in mappings
        if mapping.execution_status is not OEExecutionStatus.GENE_LEVEL_EXECUTABLE
    )
    missing = tuple(
        dict.fromkeys(
            item
            for mapping in mappings
            for item in mapping.missing_information
        )
    )
    warnings = list(proxy_plan.warnings)

    if dose.dose_mode is OEDoseMode.CATEGORICAL_ONLY:
        plan = OECapacityPlan(
            gene_id=gene_id,
            target_id=target_id,
            context_id=context_id,
            requested_dose=dose,
            execution_mode=OEExecutionMode.NOT_EXECUTABLE,
            execution_status=OEExecutionStatus.CATEGORICAL_DOSE_ONLY,
            explain_only_mappings=explain_only,
            proxy_reaction_ids=proxy_reactions,
            missing_information=tuple(dict.fromkeys((*missing, "reviewed_numeric_dose_mapping"))),
            warnings=tuple(
                dict.fromkeys(
                    (
                        *warnings,
                        "Categorical dose is explain-only until a reviewed numeric mapping exists.",
                    )
                )
            ),
            structural_mappings=mappings,
            relative_scenario_specs=(),
            product_mode=OEProductMode.NOT_EXECUTABLE,
            product_state=OEProductState.NOT_EXECUTABLE,
            absolute_capacity_availability=(
                AbsoluteCapacityAvailability.UNAVAILABLE_MISSING_REVIEWED_ANCHOR
            ),
            calibration_status=OECalibrationStatus.NOT_APPLICABLE,
            absolute_solver_allowed=False,
            model_fingerprint=model_fingerprint,
            limitations=(
                "categorical_dose_has_no_reviewed_numeric_mapping",
                "no_absolute_capacity_without_reviewed_anchor",
                "no_mg_per_litre_prediction",
            ),
        )
        plan.validate()
        return plan

    specs = build_gene_capacity_specs(gene_id, catalog, dose, parameter_policy)
    relative_specs = build_relative_oe_scenario_specs(
        gene_id,
        catalog,
        dose,
        parameter_policy,
    )
    executable_mapping_ids = {
        mapping.mapping_id
        for mapping in mappings
        if mapping.execution_status is OEExecutionStatus.GENE_LEVEL_EXECUTABLE
    }
    spec_mapping_ids = {spec.mapping.mapping_id for spec in specs}
    relative_mapping_ids = {spec.mapping.mapping_id for spec in relative_specs}
    incomplete_relative_mappings = bool(
        executable_mapping_ids - relative_mapping_ids
    )
    incomplete_expected_mappings = bool(executable_mapping_ids - spec_mapping_ids)
    if incomplete_expected_mappings:
        missing = tuple(dict.fromkeys((*missing, "capacity_parameters")))
        for mapping in mappings:
            if (
                mapping.mapping_id in spec_mapping_ids
                or mapping.mapping_id not in executable_mapping_ids
            ):
                continue
            selected = _select_parameter_set(
                mapping.mapping_id,
                mapping.gene_id,
                mapping.enzyme_id,
                parameter_policy,
            )
            if selected is None or selected.baseline_enzyme_amount is None:
                missing = tuple(
                    dict.fromkeys((*missing, "reviewed_baseline_capacity"))
                )
        # Absolute execution is an all-expected-mappings contract. A partially
        # anchored plan remains useful for relative review, but none of its
        # absolute specs may reach the solver.
        specs = ()
    if executable_mapping_ids and not specs:
        missing = tuple(
            dict.fromkeys(
                (*missing, "reviewed_baseline_capacity", "capacity_parameters")
            )
        )
    scenarios = tuple(dict.fromkeys(spec.parameter_scenario for spec in specs))
    if specs and proxy_reactions:
        mode = OEExecutionMode.COMPARISON
        status = (
            OEExecutionStatus.PARTIAL_MAPPING
            if explain_only or incomplete_expected_mappings
            else OEExecutionStatus.GENE_LEVEL_EXECUTABLE
        )
    elif specs:
        mode = OEExecutionMode.GENE_CAPACITY
        status = (
            OEExecutionStatus.PARTIAL_MAPPING
            if explain_only or incomplete_expected_mappings
            else OEExecutionStatus.GENE_LEVEL_EXECUTABLE
        )
    elif relative_specs:
        mode = OEExecutionMode.RELATIVE_GENE_CAPACITY
        status = (
            OEExecutionStatus.PARTIAL_MAPPING
            if explain_only or incomplete_relative_mappings
            else OEExecutionStatus.GENE_LEVEL_EXECUTABLE
        )
        warnings.append(
            "相对gene-capacity结果未经校准，不构成绝对model_flux锚点。"
        )
    elif proxy_reactions:
        mode = OEExecutionMode.REACTION_PROXY
        status = (
            OEExecutionStatus.PARTIAL_MAPPING
            if relative_specs and (explain_only or incomplete_expected_mappings)
            else OEExecutionStatus.PROXY_ONLY
        )
        warnings.append(
            "缺少审核过的baseline锚点，绝对gene-capacity不可用；"
            "reaction proxy是另一条独立的方向性证据路径。"
        )
    else:
        mode = OEExecutionMode.NOT_EXECUTABLE
        status = _non_executable_status(mappings)
        if not mappings:
            missing = tuple(dict.fromkeys((*missing, "gene_mapping")))

    plan = OECapacityPlan(
        gene_id=gene_id,
        target_id=target_id,
        context_id=context_id,
        requested_dose=dose,
        execution_mode=mode,
        execution_status=status,
        executable_capacity_specs=specs,
        explain_only_mappings=explain_only,
        proxy_reaction_ids=proxy_reactions,
        uncertainty_scenarios=scenarios,
        missing_information=missing,
        warnings=tuple(dict.fromkeys(warnings)),
        structural_mappings=mappings,
        available_relative_scenario_specs=relative_specs,
        relative_scenario_specs=relative_specs if not specs else (),
        product_mode=(
            OEProductMode.ABSOLUTE_CAPACITY
            if specs
            else OEProductMode.RELATIVE_UNCALIBRATED
            if relative_specs
            else OEProductMode.REACTION_PROXY
            if proxy_reactions
            else OEProductMode.NOT_EXECUTABLE
        ),
        product_state=(
            OEProductState.ABSOLUTE_AVAILABLE
            if specs
            else OEProductState.RELATIVE_UNCALIBRATED
            if relative_specs
            else OEProductState.REACTION_PROXY
            if proxy_reactions
            else OEProductState.NOT_EXECUTABLE
        ),
        absolute_capacity_availability=(
            AbsoluteCapacityAvailability.AVAILABLE_REVIEWED
            if specs
            else AbsoluteCapacityAvailability.UNAVAILABLE_MISSING_REVIEWED_ANCHOR
        ),
        calibration_status=(
            OECalibrationStatus.REVIEWED_ABSOLUTE
            if specs
            else OECalibrationStatus.RELATIVE_UNCALIBRATED
            if relative_specs
            else OECalibrationStatus.PROXY_ONLY
            if proxy_reactions
            else OECalibrationStatus.NOT_APPLICABLE
        ),
        absolute_solver_allowed=bool(specs),
        model_fingerprint=model_fingerprint,
        limitations=(
            "model_relative_only",
            "no_mg_per_litre_prediction",
            "no_experimental_success_probability",
            "no_cross_context_guarantee",
        ),
    )
    plan.validate()
    return plan


def _reviewed_dose_mapping(
    dose_mapping: Mapping[str, Any],
    *,
    dose_id: str,
    promoter: str,
    copy_number: float | None,
) -> Mapping[str, Any] | None:
    keys = [dose_id]
    if promoter:
        copy_label = "" if copy_number is None else f"{copy_number:g}"
        keys.extend((f"{promoter}|{copy_label}", promoter))
    for key in keys:
        if key and key in dose_mapping:
            item = dose_mapping[key]
            if isinstance(item, Mapping):
                return item
            return {"expression_multiplier": item, "mapping_source": f"dose_mapping:{key}"}
    return None


def _matching_capacity_anchor(
    catalog: CapacityAnchorCatalog | None,
    *,
    mapping: Any,
    target_id: str,
    context_id: str,
) -> CapacityAnchor | None:
    if catalog is None or not target_id or not context_id:
        return None
    matches = tuple(
        anchor
        for anchor in catalog.anchors
        if anchor.target_id == target_id
        and anchor.context_id == context_id
        and anchor.gene_id == mapping.gene_id
        and anchor.enzyme_id == mapping.enzyme_id
        and anchor.formation_or_dilution_reaction_id
        == mapping.formation_or_dilution_reaction_id
        and anchor.model_fingerprint == mapping.model_fingerprint
    )
    if len(matches) > 1:
        raise OECapacityValidationError(
            f"multiple capacity anchors match mapping {mapping.mapping_id}."
        )
    return matches[0] if matches else None


def _optional_float(value: Any, field_name: str) -> float | None:
    if value is None or value == "":
        return None
    return _required_positive_float(value, field_name)


def _required_positive_float(value: Any, field_name: str) -> float:
    if isinstance(value, bool):
        raise OECapacityValidationError(f"{field_name} must be a positive finite number.")
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise OECapacityValidationError(
            f"{field_name} must be a positive finite number."
        ) from exc
    if parsed <= 0 or parsed == float("inf") or parsed == float("-inf") or parsed != parsed:
        raise OECapacityValidationError(f"{field_name} must be a positive finite number.")
    return parsed


def _exact_parameter(
    combined: Any,
    *,
    method_name: str,
    enzyme_id: str,
    parameter_name: str,
    unit: str,
    source_ref: str,
    relative_uncertainty: float,
    source_version: str,
    warnings: list[str],
) -> ParameterEstimate | None:
    method = getattr(combined, method_name, None)
    if method is None:
        warnings.append(f"Combined enzyme data does not expose {method_name}.")
        return None
    try:
        value = _required_positive_float(method(enzyme_id), parameter_name)
    except (KeyError, OECapacityValidationError) as exc:
        warnings.append(str(exc))
        return None
    return ParameterEstimate(
        parameter_name=parameter_name,
        nominal_value=value,
        lower_bound=value * (1.0 - relative_uncertainty),
        upper_bound=value * (1.0 + relative_uncertainty),
        unit=unit,
        source_type=EvidenceSourceType.LOCAL_ENZYME_DATA,
        source_ref=source_ref,
        source_version=source_version,
        confidence=ConfidenceLevel.HIGH,
    )


def _select_parameter_set(
    mapping_id: str,
    gene_id: str,
    enzyme_id: str,
    policy: ParameterPolicy,
) -> GeneCapacityParameterSet | None:
    candidates = tuple(
        item
        for item in policy.parameter_sets
        if item.mapping_id == mapping_id
    )
    if not candidates:
        return None
    complete = tuple(item for item in candidates if not item.missing_information)
    ranked = sorted(
        complete or candidates,
        key=lambda item: (_parameter_priority(item), item.parameter_set_id),
    )
    best_priority = _parameter_priority(ranked[0])
    best = tuple(item for item in ranked if _parameter_priority(item) == best_priority)
    signatures = {_parameter_signature(item) for item in best}
    if len(signatures) > 1 and policy.strict_conflicts:
        raise OECapacityParameterConflictError(
            f"conflicting parameter sets at the same evidence priority for {mapping_id}."
        )
    return best[0]


def _parameter_priority(parameter_set: GeneCapacityParameterSet) -> int:
    estimates = tuple(
        item
        for item in (
            parameter_set.kcat,
            parameter_set.molecular_weight,
            parameter_set.baseline_enzyme_amount,
            parameter_set.complex_stoichiometry,
        )
        if item is not None
    )
    if not estimates:
        return 999
    return max(_SOURCE_PRIORITY[item.source_type] for item in estimates)


def _parameter_signature(parameter_set: GeneCapacityParameterSet) -> tuple[object, ...]:
    def signature(estimate: Any) -> object:
        if estimate is None:
            return None
        return (
            estimate.nominal_value,
            estimate.lower_bound,
            estimate.upper_bound,
            estimate.unit,
        )

    return (
        signature(parameter_set.kcat),
        signature(parameter_set.molecular_weight),
        signature(parameter_set.baseline_enzyme_amount),
        signature(parameter_set.complex_stoichiometry),
    )


def _non_executable_status(
    mappings: tuple[Any, ...],
) -> OEExecutionStatus:
    statuses = {mapping.execution_status for mapping in mappings}
    for status in (
        OEExecutionStatus.EXTERNAL_EVIDENCE_ONLY,
        OEExecutionStatus.COMPLEX_LIMITED,
        OEExecutionStatus.ISOENZYME_AMBIGUOUS,
        OEExecutionStatus.PARTIAL_MAPPING,
        OEExecutionStatus.UNRESOLVED,
    ):
        if status in statuses:
            return status
    return OEExecutionStatus.UNRESOLVED


_SOURCE_PRIORITY = {
    EvidenceSourceType.CURRENT_MODEL: 0,
    EvidenceSourceType.LOCAL_ENZYME_DATA: 0,
    EvidenceSourceType.REVIEWED_PICHIA_MAPPING: 1,
    EvidenceSourceType.EXTERNAL_PICHIA_MODEL: 2,
    EvidenceSourceType.PICHIA_LITERATURE: 3,
    EvidenceSourceType.HOMOLOGY_TRANSFER: 4,
    EvidenceSourceType.SMOKE_FIXTURE: 5,
}


__all__ = [
    "build_current_model_parameter_policy",
    "build_gene_capacity_specs",
    "build_oe_dose_spec",
    "load_capacity_anchor_catalog",
    "plan_gene_level_overexpression",
]
