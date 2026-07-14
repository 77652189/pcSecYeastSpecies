from __future__ import annotations

from pathlib import Path
from typing import Mapping, Sequence

from pcsec_pichia.oe_capacity.external_candidate_io import (
    _sha256_file,
    _source_artifact_matches,
)
from pcsec_pichia.oe_capacity.external_candidate_schema import (
    CapacityApplicabilityScope,
    CapacityCandidateStatus,
    CapacityConfidence,
    CapacityConversionStep,
    CapacityModelBinding,
    CapacityParameterKind,
    ExternalCapacityCandidate,
    ExternalCapacitySource,
    ExternalCapacitySourceType,
    HostCondition,
    RawCapacityMeasurement,
    _condition_context_id,
)
from pcsec_pichia.oe_capacity.schema import (
    GeneCapacityCatalog,
    OECapacityValidationError,
)


def capacity_context_id(condition: HostCondition) -> str:
    """Return the canonical current-model context identifier for a condition."""

    condition.validate()
    return _condition_context_id(condition)


def build_capacity_model_binding(
    catalog: GeneCapacityCatalog,
    *,
    target_id: str,
    context_id: str,
    gene_id: str,
    mapping_id: str = "",
    external_gene_id: str = "",
    external_protein_id: str = "",
    external_enzyme_id: str = "",
    mapping_evidence: Sequence[str] = (),
) -> CapacityModelBinding:
    catalog.validate()
    matches = tuple(
        item
        for item in catalog.mappings
        if item.gene_id == gene_id and (not mapping_id or item.mapping_id == mapping_id)
    )
    if len(matches) != 1:
        raise OECapacityValidationError(
            "capacity model binding requires exactly one current-model mapping."
        )
    mapping = matches[0]
    binding = CapacityModelBinding(
        target_id=target_id,
        context_id=context_id,
        mapping_id=mapping.mapping_id,
        model_fingerprint=catalog.model_fingerprint,
        gene_id=mapping.gene_id,
        enzyme_id=mapping.enzyme_id,
        reaction_id=mapping.reaction_id,
        formation_or_dilution_reaction_id=mapping.formation_or_dilution_reaction_id,
        mapping_evidence=tuple(
            dict.fromkeys(
                (
                    "current_model_mapping",
                    f"gpr_role:{mapping.gpr_role.value}",
                    *mapping_evidence,
                )
            )
        ),
        external_gene_id=external_gene_id,
        external_protein_id=external_protein_id,
        external_enzyme_id=external_enzyme_id,
    )
    binding.validate()
    return binding


def build_capacity_candidate(
    *,
    candidate_id: str,
    applicability_scope: CapacityApplicabilityScope,
    model_bindings: Sequence[CapacityModelBinding],
    catalogs: Mapping[str, GeneCapacityCatalog],
    sources: Sequence[ExternalCapacitySource],
    condition: HostCondition,
    abundance: RawCapacityMeasurement | None = None,
    kcat: RawCapacityMeasurement | None = None,
    direct_capacity: RawCapacityMeasurement | None = None,
    target_id: str = "",
    confidence: CapacityConfidence = CapacityConfidence.MEDIUM,
) -> ExternalCapacityCandidate:
    conflicts: list[str] = []
    missing: list[str] = []
    steps: list[CapacityConversionStep] = []
    measurements = tuple(item for item in (abundance, kcat, direct_capacity) if item is not None)
    if any(item.condition != condition for item in measurements):
        conflicts.append("candidate_measurement_condition_mismatch")
    if direct_capacity is not None and (abundance is not None or kcat is not None):
        conflicts.append("direct_and_derived_capacity_both_present")
    if direct_capacity is not None:
        if direct_capacity.parameter_kind is not CapacityParameterKind.BASELINE_CAPACITY:
            conflicts.append("direct_capacity_parameter_kind_mismatch")
        factor = _capacity_unit_factor(direct_capacity.unit, direct_capacity.condition.biomass_basis)
        if factor is None:
            missing.append("supported_capacity_unit_conversion")
            low = nominal = high = None
        else:
            low = direct_capacity.lower_bound * factor
            nominal = direct_capacity.nominal_value * factor
            high = direct_capacity.upper_bound * factor
            steps.append(
                CapacityConversionStep(
                    step_id="direct-capacity-to-model-flux",
                    input_value=direct_capacity.nominal_value,
                    input_unit=direct_capacity.unit,
                    output_value=nominal,
                    output_unit="model_flux",
                    formula="input_capacity * unit_factor",
                    factor=factor,
                    source_ref=direct_capacity.measurement_id,
                )
            )
    elif abundance is not None and kcat is not None:
        abundance_factor = _abundance_unit_factor(abundance.unit, abundance.biomass_basis or abundance.condition.biomass_basis)
        kcat_factor = _kcat_unit_factor(kcat.unit)
        if abundance_factor is None:
            missing.append("abundance_unit_or_biomass_conversion")
        if kcat_factor is None:
            missing.append("kcat_unit_conversion")
        if abundance.condition != kcat.condition:
            conflicts.append("abundance_kcat_condition_mismatch")
        if abundance_factor is None or kcat_factor is None:
            low = nominal = high = None
        else:
            low = abundance.lower_bound * abundance_factor * kcat.lower_bound * kcat_factor
            nominal = abundance.nominal_value * abundance_factor * kcat.nominal_value * kcat_factor
            high = abundance.upper_bound * abundance_factor * kcat.upper_bound * kcat_factor
            steps.extend(
                (
                    CapacityConversionStep(
                        step_id="abundance-to-mmol-per-gdw",
                        input_value=abundance.nominal_value,
                        input_unit=abundance.unit,
                        output_value=abundance.nominal_value * abundance_factor,
                        output_unit="mmol_enzyme/gDW",
                        formula="abundance * unit_factor",
                        factor=abundance_factor,
                        source_ref=abundance.measurement_id,
                    ),
                    CapacityConversionStep(
                        step_id="kcat-to-per-hour",
                        input_value=kcat.nominal_value,
                        input_unit=kcat.unit,
                        output_value=kcat.nominal_value * kcat_factor,
                        output_unit="1/h",
                        formula="kcat * time_factor",
                        factor=kcat_factor,
                        source_ref=kcat.measurement_id,
                    ),
                    CapacityConversionStep(
                        step_id="abundance-times-kcat",
                        input_value=nominal,
                        input_unit="mmol_enzyme/gDW * 1/h",
                        output_value=nominal,
                        output_unit="model_flux",
                        formula="abundance_mmol_per_gdw * kcat_per_h",
                        factor=1.0,
                        source_ref=f"{abundance.measurement_id}+{kcat.measurement_id}",
                    ),
                )
            )
    else:
        missing.append("baseline_capacity_or_abundance_and_kcat")
        low = nominal = high = None
    resolved_bindings = tuple(model_bindings)
    for binding in resolved_bindings:
        _validate_binding_against_catalog(binding, catalogs.get(binding.target_id))
        if binding.context_id != _condition_context_id(condition):
            raise OECapacityValidationError(
                "capacity binding context_id does not match the structured host condition."
            )
    source_lookup = {source.source_id: source for source in sources}
    referenced_source_ids = tuple(dict.fromkeys(item.source_id for item in measurements))
    if not referenced_source_ids or any(
        source_id not in source_lookup for source_id in referenced_source_ids
    ):
        missing.append("traceable_capacity_sources")
    referenced_sources = tuple(
        source_lookup[source_id]
        for source_id in referenced_source_ids
        if source_id in source_lookup
    )
    if any(
        source.source_type is ExternalCapacitySourceType.IDENTITY_REFERENCE
        for source in referenced_sources
    ):
        missing.append("capacity_valued_source")
    if (
        applicability_scope is CapacityApplicabilityScope.EXTERNAL_MODEL_CALIBRATED
        and not any(
            source.source_type is ExternalCapacitySourceType.EXTERNAL_ENZYME_MODEL
            for source in referenced_sources
        )
    ):
        missing.append("external_model_calibration_source")
    if any(not _source_artifact_matches(source) for source in referenced_sources):
        missing.append("source_artifact_hash_verification")
    source_review_ready = bool(referenced_source_ids) and all(
        source.terms_reviewed
        and source.license_id not in {"", "unreviewed"}
        and source.source_version not in {"", "unversioned"}
        and source.source_type is not ExternalCapacitySourceType.IDENTITY_REFERENCE
        and _source_artifact_matches(source)
        and not {
            "source_version_review_required",
            "license_review_required",
            "expected_sha256_not_predeclared",
        }.intersection(source.warnings)
        for source in referenced_sources
    )
    if applicability_scope is CapacityApplicabilityScope.HOMOLOG_TRANSFERRED:
        confidence = CapacityConfidence.LOW
        missing.append("independent_pichia_capacity_evidence")
    status = CapacityCandidateStatus.REVIEW_READY
    if conflicts or missing or not source_review_ready:
        status = CapacityCandidateStatus.REVIEW_REQUIRED
    candidate = ExternalCapacityCandidate(
        candidate_id=candidate_id,
        applicability_scope=applicability_scope,
        target_id=target_id,
        source_ids=referenced_source_ids,
        measurement_ids=tuple(item.measurement_id for item in measurements),
        model_bindings=resolved_bindings,
        condition=condition,
        nominal_capacity=nominal,
        lower_capacity=low,
        upper_capacity=high,
        unit="model_flux",
        confidence=confidence,
        status=status,
        conversion_steps=tuple(steps),
        conflicts=tuple(dict.fromkeys(conflicts)),
        missing_information=tuple(dict.fromkeys(missing)),
        warnings=(() if source_review_ready else ("source_license_version_or_terms_requires_review",)),
    )
    candidate.validate()
    return candidate


def _validate_binding_against_catalog(
    binding: CapacityModelBinding,
    catalog: GeneCapacityCatalog | None,
) -> None:
    binding.validate()
    if catalog is None:
        raise OECapacityValidationError(
            f"missing current-model catalog for target {binding.target_id}."
        )
    catalog.validate()
    if catalog.model_fingerprint != binding.model_fingerprint:
        raise OECapacityValidationError("capacity binding model_fingerprint mismatch.")
    matches = tuple(
        item for item in catalog.mappings if item.mapping_id == binding.mapping_id
    )
    if len(matches) != 1:
        raise OECapacityValidationError("capacity binding mapping_id is not in the current catalog.")
    mapping = matches[0]
    actual = (
        mapping.gene_id,
        mapping.enzyme_id,
        mapping.reaction_id,
        mapping.formation_or_dilution_reaction_id,
    )
    expected = (
        binding.gene_id,
        binding.enzyme_id,
        binding.reaction_id,
        binding.formation_or_dilution_reaction_id,
    )
    if actual != expected:
        raise OECapacityValidationError("capacity binding identity does not match the current catalog.")

def _capacity_unit_factor(unit: str, biomass_basis: str) -> float | None:
    normalized = unit.strip().lower().replace(" ", "")
    if normalized in {"model_flux", "mmol/gdw/h", "mmol/(gdw*h)"} and biomass_basis.lower() == "gdw":
        return 1.0
    return None


def _abundance_unit_factor(unit: str, biomass_basis: str) -> float | None:
    if biomass_basis.lower() != "gdw":
        return None
    normalized = unit.strip().lower().replace(" ", "")
    return {
        "mmol_enzyme/gdw": 1.0,
        "umol_enzyme/gdw": 1e-3,
        "nmol_enzyme/gdw": 1e-6,
    }.get(normalized)


def _kcat_unit_factor(unit: str) -> float | None:
    return {"1/h": 1.0, "h^-1": 1.0, "1/s": 3600.0, "s^-1": 3600.0}.get(unit.strip().lower().replace(" ", ""))
