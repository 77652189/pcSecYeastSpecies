from __future__ import annotations

from typing import Any, Callable, Mapping

import numpy as np

from pcsec_pichia.analysis.shadow_lp.constraint_spec import ConstraintBlock, ConstraintSpec, ShadowConstraintConfig
from pcsec_pichia.analysis.shadow_lp.model_adapter import ShadowTargetPreparation


LayerBuilder = Callable[[ShadowTargetPreparation, ShadowConstraintConfig], ConstraintBlock]

REFERENCE_LAYER_ORDER: tuple[str, ...] = (
    "metabolic_coupling",
    "secretory_coupling",
    "protein_mass",
    "proteasome",
    "ribosome_assembly",
    "ribosome_translation",
    "misfolding",
    "mitochondrial",
)


def build_shadow_constraint_blocks(
    prep: ShadowTargetPreparation,
    config: ShadowConstraintConfig | None = None,
) -> tuple[ConstraintBlock, ...]:
    """Build symbolic pcSec shadow constraint blocks without assembling or solving an LP."""

    resolved_config = config or ShadowConstraintConfig()
    return tuple(builder(prep, resolved_config) for builder in _BUILDERS)


def build_metabolic_coupling_block(
    prep: ShadowTargetPreparation,
    config: ShadowConstraintConfig | None = None,
) -> ConstraintBlock:
    resolved_config = config or ShadowConstraintConfig()
    model = prep.fixed_model
    reaction_index = _reaction_index(model)
    constraints: list[ConstraintSpec] = []
    mapped: set[str] = set()
    warnings: list[str] = []

    for enzyme_id, kcat in zip(prep.metabolic.enzymes, prep.metabolic.kcat):
        reaction_id = str(enzyme_id).replace("_complex", "")
        formation_id = f"{enzyme_id}_formation"
        missing = [candidate for candidate in (reaction_id, formation_id) if candidate not in reaction_index]
        if missing:
            warnings.append(
                f"metabolic_coupling missing mapping for enzyme {enzyme_id}: {', '.join(missing)}"
            )
            continue
        coefficient = float(kcat) / resolved_config.growth_rate
        constraints.append(
            ConstraintSpec(
                name=f"metabolic_coupling:{enzyme_id}",
                layer="metabolic_coupling",
                sense="eq",
                terms={reaction_id: 1.0, formation_id: -coefficient},
                rhs=0.0,
                source="pcSec metabolic enzymedata kcat divided by ShadowConstraintConfig.growth_rate",
                enabled_by_default=True,
                metadata={
                    "enzyme_id": str(enzyme_id),
                    "reaction_id": reaction_id,
                    "formation_reaction_id": formation_id,
                    "kcat": float(kcat),
                    "growth_rate": resolved_config.growth_rate,
                },
            )
        )
        mapped.update((reaction_id, formation_id))

    return _block("metabolic_coupling", constraints, mapped, warnings)


def build_secretory_coupling_block(
    prep: ShadowTargetPreparation,
    config: ShadowConstraintConfig | None = None,
) -> ConstraintBlock:
    resolved_config = config or ShadowConstraintConfig()
    model = prep.fixed_model
    reaction_index = _reaction_index(model)
    constraints: list[ConstraintSpec] = []
    mapped: set[str] = set()
    warnings: list[str] = []

    for entry in prep.secretory.unique_complex_entries():
        formation_id = _secretory_formation_reaction_id(prep.secretory, entry.complex_id)
        if formation_id not in reaction_index:
            warnings.append(f"secretory_coupling missing formation mapping for complex {entry.complex_id}: {formation_id}")
            continue

        suffix = f"_{entry.complex_id}"
        terms: dict[str, float] = {}
        for reaction_id in _reaction_ids(model):
            if not str(reaction_id).endswith(suffix):
                continue
            coefficient = prep.secretory.reaction_coefficients.get(reaction_id)
            if coefficient is None:
                continue
            terms[str(reaction_id)] = terms.get(str(reaction_id), 0.0) + float(coefficient)
            mapped.add(str(reaction_id))

        if not terms:
            warnings.append(f"secretory_coupling found no mapped reaction coefficients for complex {entry.complex_id}")
            continue

        coefficient = float(entry.kcat) / resolved_config.growth_rate
        terms[formation_id] = terms.get(formation_id, 0.0) - coefficient
        mapped.add(formation_id)
        constraints.append(
            ConstraintSpec(
                name=f"secretory_coupling:{entry.complex_id}",
                layer="secretory_coupling",
                sense="eq",
                terms=terms,
                rhs=0.0,
                source="pcSec secretory enzymedata coefficients and kcat divided by growth rate",
                enabled_by_default=True,
                metadata={
                    "complex_id": entry.complex_id,
                    "compartment": entry.compartment,
                    "formation_reaction_id": formation_id,
                    "matched_reaction_count": len(terms) - 1,
                    "kcat": float(entry.kcat),
                    "growth_rate": resolved_config.growth_rate,
                },
            )
        )

    return _block("secretory_coupling", constraints, mapped, warnings)


def build_protein_mass_block(
    prep: ShadowTargetPreparation,
    config: ShadowConstraintConfig | None = None,
) -> ConstraintBlock:
    resolved_config = config or ShadowConstraintConfig()
    model = prep.fixed_model
    reaction_index = _reaction_index(model)
    terms: dict[str, float] = {}
    mapped: set[str] = set()
    warnings: list[str] = []

    for reaction_id in _reaction_ids(model):
        if "_dilution" not in reaction_id or "dummy" in reaction_id:
            continue
        try:
            molecular_weight = prep.combined.molecular_weight_for_dilution_reaction(reaction_id)
        except KeyError as exc:
            warnings.append(f"protein_mass missing molecular weight mapping for {reaction_id}: {exc}")
            continue
        terms[reaction_id] = float(molecular_weight) / 1000.0
        mapped.add(reaction_id)

    for reaction_id in ("dilute_dummy", "dilute_dummyER"):
        if reaction_id not in reaction_index:
            warnings.append(f"protein_mass missing dummy protein reaction: {reaction_id}")
            continue
        terms[reaction_id] = terms.get(reaction_id, 0.0) + 40.0
        mapped.add(reaction_id)

    constraints: list[ConstraintSpec] = []
    modeled_fraction, modeled_fraction_warning = _modeled_protein_fraction(model)
    if modeled_fraction_warning:
        warnings.append(modeled_fraction_warning)
    if terms:
        total_modeled_protein = resolved_config.total_protein_content * modeled_fraction
        constraints.append(
            ConstraintSpec(
                name="protein_mass:total_modeled_protein",
                layer="protein_mass",
                sense="eq",
                terms=terms,
                rhs=resolved_config.growth_rate * total_modeled_protein,
                source="pcSec combined enzymedata molecular weights and modeled biomass protein fraction",
                enabled_by_default=True,
                metadata={
                    "growth_rate": resolved_config.growth_rate,
                    "total_protein_content": resolved_config.total_protein_content,
                    "modeled_protein_fraction": modeled_fraction,
                    "term_count": len(terms),
                },
            )
        )
    else:
        warnings.append("protein_mass found no dilution or dummy protein terms")

    if "dilute_dummyER" in reaction_index:
        constraints.append(
            ConstraintSpec(
                name="protein_mass:unmodeled_er",
                layer="protein_mass",
                sense="eq",
                terms={"dilute_dummyER": 40.0},
                rhs=resolved_config.growth_rate * resolved_config.unmodeled_er_protein_fraction,
                source="pcSec unmodeled ER protein fraction and fixed dummy molecular weight",
                enabled_by_default=True,
                metadata={
                    "growth_rate": resolved_config.growth_rate,
                    "unmodeled_er_protein_fraction": resolved_config.unmodeled_er_protein_fraction,
                    "dummy_molecular_weight": 40.0,
                },
            )
        )

    return _block("protein_mass", constraints, mapped, warnings)


def build_proteasome_block(
    prep: ShadowTargetPreparation,
    config: ShadowConstraintConfig | None = None,
) -> ConstraintBlock:
    resolved_config = config or ShadowConstraintConfig()
    model = prep.fixed_model
    reaction_index = _reaction_index(model)
    formation_id = "Mach_proteasome_complex_formation"
    mapped: set[str] = set()
    warnings: list[str] = []
    if formation_id not in reaction_index:
        return _block("proteasome", [], mapped, [f"proteasome missing formation reaction: {formation_id}"])

    terms: dict[str, float] = {}
    for reaction_id in _reaction_ids(model):
        if reaction_id.endswith("_subunit_degradation"):
            protein_id = reaction_id.replace("_subunit_degradation", "").replace("r_", "")
            try:
                coefficient = prep.combined.exact_protein_length(protein_id) / 467.0
            except KeyError as exc:
                warnings.append(f"proteasome missing protein length mapping for {protein_id}: {exc}")
                continue
        elif reaction_id.endswith("_sp_degradation"):
            protein_id = None
            coefficient = 25.0 / 467.0
        else:
            continue
        terms[reaction_id] = float(coefficient)
        mapped.add(reaction_id)

    if not terms:
        return _block("proteasome", [], mapped, [*warnings, "proteasome found no degradation terms"])

    try:
        formation_coefficient = prep.combined.exact_enzyme_kcat("Mach_proteasome_complex") / resolved_config.growth_rate
    except KeyError as exc:
        return _block("proteasome", [], mapped, [*warnings, f"proteasome missing kcat mapping: {exc}"])

    terms[formation_id] = terms.get(formation_id, 0.0) - formation_coefficient
    mapped.add(formation_id)
    constraint = ConstraintSpec(
        name="proteasome:degradation_capacity",
        layer="proteasome",
        sense="eq",
        terms=terms,
        rhs=0.0,
        source="pcSec protein degradation lengths and proteasome kcat divided by growth rate",
        enabled_by_default=True,
        metadata={
            "formation_reaction_id": formation_id,
            "proteasome_kcat_over_growth": formation_coefficient,
            "growth_rate": resolved_config.growth_rate,
            "degradation_term_count": len(terms) - 1,
            "protein_length_denominator": 467.0,
            "signal_peptide_length": 25.0,
        },
    )
    return _block("proteasome", [constraint], mapped, warnings)


def build_ribosome_assembly_block(
    prep: ShadowTargetPreparation,
    config: ShadowConstraintConfig | None = None,
) -> ConstraintBlock:
    resolved_config = config or ShadowConstraintConfig()
    model = prep.fixed_model
    reaction_index = _reaction_index(model)
    ribosome_id = "Mach_Ribosome_complex_formation"
    assembly_id = "Mach_Ribosome_Assembly_Factors_complex_formation"
    missing = [reaction_id for reaction_id in (ribosome_id, assembly_id) if reaction_id not in reaction_index]
    if missing:
        return _block("ribosome_assembly", [], set(), [f"ribosome_assembly missing reaction mapping: {', '.join(missing)}"])

    try:
        coefficient = prep.combined.exact_enzyme_kcat("Mach_Ribosome_Assembly_Factors_complex") / resolved_config.growth_rate
    except KeyError as exc:
        return _block("ribosome_assembly", [], set(), [f"ribosome_assembly missing kcat mapping: {exc}"])

    constraint = ConstraintSpec(
        name="ribosome_assembly:assembly_factor_capacity",
        layer="ribosome_assembly",
        sense="eq",
        terms={ribosome_id: 1.0, assembly_id: -coefficient},
        rhs=0.0,
        source="pcSec ribosome assembly factor kcat divided by growth rate",
        enabled_by_default=True,
        metadata={
            "ribosome_formation_reaction_id": ribosome_id,
            "assembly_factor_formation_reaction_id": assembly_id,
            "assembly_factor_kcat_over_growth": coefficient,
            "growth_rate": resolved_config.growth_rate,
        },
    )
    return _block("ribosome_assembly", [constraint], {ribosome_id, assembly_id}, [])


def build_ribosome_translation_block(
    prep: ShadowTargetPreparation,
    config: ShadowConstraintConfig | None = None,
) -> ConstraintBlock:
    resolved_config = config or ShadowConstraintConfig()
    reason = "ribosome_translation disabled to match current pcSec reference defaults"
    if resolved_config.enable_ribosome_translation:
        reason = "ribosome_translation builders are intentionally not implemented in Phase 2"
    return ConstraintBlock(
        layer_id="ribosome_translation",
        constraints=(),
        counts={"ribosome_translation": 0},
        warnings=(reason,),
        metadata={"enabled_by_config": resolved_config.enable_ribosome_translation, "target_id": prep.target_id},
    )


def build_misfolding_block(
    prep: ShadowTargetPreparation,
    config: ShadowConstraintConfig | None = None,
) -> ConstraintBlock:
    resolved_config = config or ShadowConstraintConfig()
    reason = "misfolding disabled to match current pcSec reference defaults"
    if resolved_config.enable_misfolding:
        reason = "misfolding builders are intentionally not implemented in Phase 2"
    return ConstraintBlock(
        layer_id="misfolding",
        constraints=(),
        counts={"misfolding": 0},
        warnings=(reason,),
        metadata={"enabled_by_config": resolved_config.enable_misfolding, "target_id": prep.target_id},
    )


def build_mitochondrial_block(
    prep: ShadowTargetPreparation,
    config: ShadowConstraintConfig | None = None,
) -> ConstraintBlock:
    resolved_config = config or ShadowConstraintConfig()
    model = prep.fixed_model
    reaction_index = _reaction_index(model)
    terms: dict[str, float] = {}
    mapped: set[str] = set()
    warnings: list[str] = []

    for reaction_id in _collect_compartment_reactions(model, ("m", "mm")):
        dilution_id = f"{reaction_id}_complex_dilution"
        if dilution_id not in reaction_index:
            warnings.append(f"mitochondrial missing dilution reaction for {reaction_id}: {dilution_id}")
            continue
        enzyme_id = dilution_id.replace("_dilution", "")
        try:
            molecular_weight = prep.combined.exact_enzyme_mw(enzyme_id)
        except KeyError as exc:
            warnings.append(f"mitochondrial missing enzyme molecular weight for {enzyme_id}: {exc}")
            continue
        terms[dilution_id] = float(molecular_weight) / 1000.0
        mapped.add(dilution_id)

    if not terms:
        return _block("mitochondrial", [], mapped, [*warnings, "mitochondrial found no dilution terms"])

    rhs = resolved_config.growth_rate * resolved_config.mitochondrial_protein_fraction
    constraint = ConstraintSpec(
        name="mitochondrial:protein_mass_upper_bound",
        layer="mitochondrial",
        sense="le",
        terms=terms,
        rhs=rhs,
        source="pcSec mitochondrial compartment reactions and enzyme molecular weights",
        enabled_by_default=True,
        metadata={
            "growth_rate": resolved_config.growth_rate,
            "mitochondrial_protein_fraction": resolved_config.mitochondrial_protein_fraction,
            "term_count": len(terms),
        },
    )
    return _block("mitochondrial", [constraint], mapped, warnings)


def _block(
    layer_id: str,
    constraints: list[ConstraintSpec],
    mapped: set[str],
    warnings: list[str],
) -> ConstraintBlock:
    return ConstraintBlock(
        layer_id=layer_id,
        constraints=tuple(constraints),
        counts={layer_id: len(constraints)},
        mapped_reaction_count=len(mapped),
        missing_mapping_count=len(warnings),
        warnings=tuple(warnings),
        metadata={"constraint_count": len(constraints)},
    )


def _reaction_index(model: Any) -> Mapping[str, int]:
    reaction_index = getattr(model, "reaction_index")
    if callable(reaction_index):
        reaction_index = reaction_index()
    return reaction_index


def _reaction_ids(model: Any) -> tuple[str, ...]:
    return tuple(str(reaction_id) for reaction_id in getattr(model, "rxns"))


def _secretory_formation_reaction_id(secretory: Any, complex_id: str) -> str:
    method = getattr(secretory, "formation_reaction_id_for_complex", None)
    if callable(method):
        return str(method(complex_id))
    return f"{complex_id}_formation"


def _modeled_protein_fraction(model: Any) -> tuple[float, str | None]:
    method = getattr(model, "modeled_protein_fraction", None)
    if callable(method):
        try:
            return float(method("BIOMASS", "PROTEIN[c]")), None
        except TypeError:
            return float(method()), None
        except KeyError as exc:
            return 1.0, f"protein_mass could not resolve modeled protein fraction: {exc}"

    try:
        metabolite_index = list(getattr(model, "mets")).index("PROTEIN[c]")
        reaction_index = _reaction_index(model)["BIOMASS"]
    except (KeyError, ValueError) as exc:
        return 1.0, f"protein_mass could not resolve modeled protein fraction: {exc}"
    return 1.0 + float(model.s_matrix[metabolite_index, reaction_index]), None


def _collect_compartment_reactions(model: Any, compartment_ids: tuple[str, ...]) -> list[str]:
    tokens = tuple(f"[{compartment_id}]" for compartment_id in compartment_ids)
    metabolite_indices = [
        index
        for index, metabolite_id in enumerate(getattr(model, "mets"))
        if any(token in str(metabolite_id) for token in tokens)
    ]
    if not metabolite_indices:
        return []

    submatrix = model.s_matrix[metabolite_indices, :].tocsc()
    has_stoichiometry = np.diff(submatrix.indptr) > 0
    rules = tuple(getattr(model, "rules", ()))
    return [
        reaction_id
        for index, reaction_id in enumerate(_reaction_ids(model))
        if has_stoichiometry[index] and index < len(rules) and str(rules[index]).strip() not in {"", "[]"}
    ]


_BUILDERS: tuple[LayerBuilder, ...] = (
    build_metabolic_coupling_block,
    build_secretory_coupling_block,
    build_protein_mass_block,
    build_proteasome_block,
    build_ribosome_assembly_block,
    build_ribosome_translation_block,
    build_misfolding_block,
    build_mitochondrial_block,
)
