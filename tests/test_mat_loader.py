from __future__ import annotations

import json

import pytest

from pcsec_pichia.adapters.mat_loader import MatStructLoader
from pcsec_pichia.core.pichia_model import collect_compartment_reactions
from pcsec_pichia.core.paths import ProjectPaths


def test_mat_loader_reads_pcsec_pichia_model_summary() -> None:
    loader = MatStructLoader(ProjectPaths.discover())

    summary = loader.load_pcsec_pichia_summary()

    assert summary.reaction_count == 29026
    assert summary.metabolite_count == 20195
    assert summary.gene_count == 1025
    assert summary.stoichiometric_shape == (20195, 29026)
    assert summary.rxn_gene_shape == (29026, 1025)
    assert all(summary.key_reactions.values())


def test_pichia_model_can_change_reaction_bounds_without_matlab() -> None:
    model = MatStructLoader(ProjectPaths.discover()).load_pcsec_pichia_model()
    index = model.reaction_index["Ex_glc_D"]

    changed = model.change_rxn_bounds("Ex_glc_D", lower=-5.0, upper=2.0)

    assert model.lb[index] != -5.0 or model.ub[index] != 2.0
    assert changed.lb[index] == -5.0
    assert changed.ub[index] == 2.0


def test_mat_loader_exports_model_and_enzymedata_summary(tmp_path) -> None:
    output = tmp_path / "pichia_model_summary.json"

    MatStructLoader(ProjectPaths.discover()).write_model_summary_json(output)

    payload = json.loads(output.read_text(encoding="utf-8"))
    assert payload["model"]["reaction_count"] == 29026
    assert {item["variable_name"] for item in payload["enzymedata"]} == {"enzymedataSEC", "enzymedataMachine"}
    assert all(item["enzyme_count"] > 0 for item in payload["enzymedata"])


def test_mat_loader_reads_metabolic_enzymedata_for_coupling_constraints() -> None:
    loader = MatStructLoader(ProjectPaths.discover())
    model = loader.load_pcsec_pichia_model()
    enzymedata = loader.load_metabolic_enzymedata()

    assert enzymedata.enzyme_count == 2732
    assert len(enzymedata.kcat) == 2732
    assert all(
        enzymedata.reaction_id_for_enzyme(enzyme_id) in model.reaction_index
        and enzymedata.formation_reaction_id_for_enzyme(enzyme_id) in model.reaction_index
        for enzyme_id in enzymedata.enzymes
    )


def test_mat_loader_reads_secretory_enzymedata_for_coupling_constraints() -> None:
    loader = MatStructLoader(ProjectPaths.discover())
    model = loader.load_pcsec_pichia_model()
    enzymedata = loader.load_secretory_enzymedata()

    entries = enzymedata.unique_complex_entries()
    assert len(entries) == 58
    assert len(enzymedata.reaction_coefficients) == 8224
    assert all(enzymedata.formation_reaction_id_for_complex(entry.complex_id) in model.reaction_index for entry in entries)
    assert all(
        any(reaction_id.endswith(f"_{entry.complex_id}") for reaction_id in enzymedata.reaction_coefficients)
        for entry in entries
    )


def test_mat_loader_combines_enzymedata_for_protein_mass_constraints() -> None:
    loader = MatStructLoader(ProjectPaths.discover())
    model = loader.load_pcsec_pichia_model()
    combined = loader.load_combined_enzymedata()
    dilution_rxns = [reaction_id for reaction_id in model.rxns if "_dilution" in reaction_id and "dummy" not in reaction_id]

    assert combined.enzyme_count == 2793
    assert combined.protein_count == 1417
    assert len(combined.kcat) == 2793
    assert len(combined.protein_length) == 1417
    assert len(dilution_rxns) == 4210
    assert all(combined.molecular_weight_for_dilution_reaction(reaction_id) > 0 for reaction_id in dilution_rxns)
    assert combined.exact_enzyme_kcat("Mach_proteasome_complex") == pytest.approx(60.0)
    assert model.modeled_protein_fraction("BIOMASS", "PROTEIN[c]") == pytest.approx(0.879787)


def test_collect_compartment_reactions_matches_mitochondrial_dilution_set() -> None:
    loader = MatStructLoader(ProjectPaths.discover())
    model = loader.load_pcsec_pichia_model()
    combined = loader.load_combined_enzymedata()
    mito_reactions = collect_compartment_reactions(model, ("m", "mm"))
    mito_dilution_reactions = [f"{reaction_id}_complex_dilution" for reaction_id in mito_reactions]

    assert len(mito_reactions) == 522
    assert all(reaction_id in model.reaction_index for reaction_id in mito_dilution_reactions)
    assert all(
        combined.exact_enzyme_mw(reaction_id.replace("_dilution", "")) > 0
        for reaction_id in mito_dilution_reactions
    )
    assert combined.exact_enzyme_kcat("Mach_Ribosome_Assembly_Factors_complex") == pytest.approx(120000.0)


def test_combined_enzymedata_maps_degradation_reactions_to_protein_lengths() -> None:
    loader = MatStructLoader(ProjectPaths.discover())
    model = loader.load_pcsec_pichia_model()
    combined = loader.load_combined_enzymedata()
    degradation_reactions = [
        reaction_id
        for reaction_id in model.rxns
        if reaction_id.endswith("_subunit_degradation") or reaction_id.endswith("_sp_degradation")
    ]
    subunit_degradation_reactions = [
        reaction_id for reaction_id in degradation_reactions if reaction_id.endswith("_subunit_degradation")
    ]

    assert len(degradation_reactions) == 1417
    assert subunit_degradation_reactions
    assert all(
        combined.exact_protein_length(
            reaction_id.replace("_subunit_degradation", "").replace("r_", "")
        )
        > 0
        for reaction_id in subunit_degradation_reactions
    )
