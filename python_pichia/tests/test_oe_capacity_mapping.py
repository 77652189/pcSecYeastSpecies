from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np

from pcsec_pichia.core.pichia_enzymes import CombinedEnzymeData, MetabolicEnzymeData
from pcsec_pichia.oe_capacity import (
    OEExecutionStatus,
    build_gene_enzyme_reaction_catalog,
    summarize_gene_capacity_catalog,
)


def test_current_model_single_gene_maps_to_exact_enzyme_formation_handle() -> None:
    model = SimpleNamespace(
        rxns=["R1", "R1_complex_formation"],
        rules=["x(1)", ""],
        gr_rules=["G1", ""],
        genes=["G1"],
        gene_index={"G1": 0},
    )
    metabolic = MetabolicEnzymeData(
        source_file=Path("Enzymedata/metabolic.mat"),
        enzymes=["R1_complex"],
        kcat=np.array([120.0]),
    )
    combined = CombinedEnzymeData(
        source_files=(Path("Enzymedata/combined.mat"),),
        enzymes=["R1_complex"],
        kcat=np.array([120.0]),
        enzyme_mw=np.array([60.0]),
        proteins=[],
        protein_length=np.array([]),
        protein_mw=np.array([]),
    )

    catalog = build_gene_enzyme_reaction_catalog(model, metabolic, combined)
    catalog.validate()

    assert len(catalog.mappings) == 1
    mapping = catalog.mappings[0]
    assert mapping.gene_id == "G1"
    assert mapping.reaction_id == "R1"
    assert mapping.enzyme_id == "R1_complex"
    assert mapping.formation_or_dilution_reaction_id == "R1_complex_formation"
    assert mapping.execution_status is OEExecutionStatus.GENE_LEVEL_EXECUTABLE
    coverage = summarize_gene_capacity_catalog(catalog)
    assert coverage.total_mappings == 1
    assert coverage.gene_count == 1
    assert dict(coverage.by_status) == {"gene_level_executable": 1}


def test_current_model_gene_without_gpr_is_preserved_as_unresolved() -> None:
    model = SimpleNamespace(
        rxns=["R1"],
        rules=[""],
        gr_rules=[""],
        genes=["G1"],
        gene_index={"G1": 0},
    )
    metabolic = MetabolicEnzymeData(
        source_file=Path("Enzymedata/metabolic.mat"),
        enzymes=[],
        kcat=np.array([]),
    )
    combined = CombinedEnzymeData(
        source_files=(Path("Enzymedata/combined.mat"),),
        enzymes=[],
        kcat=np.array([]),
        enzyme_mw=np.array([]),
        proteins=[],
        protein_length=np.array([]),
        protein_mw=np.array([]),
    )

    catalog = build_gene_enzyme_reaction_catalog(model, metabolic, combined)

    assert len(catalog.mappings) == 1
    mapping = catalog.mappings[0]
    assert mapping.execution_status is OEExecutionStatus.UNRESOLVED
    assert mapping.reaction_id == ""
    assert "model_gpr_rule" in mapping.missing_information
