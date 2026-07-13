from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np

from pcsec_pichia.core.pichia_enzymes import CombinedEnzymeData, MetabolicEnzymeData
from pcsec_pichia.oe_capacity import (
    OEExecutionStatus,
    build_gene_enzyme_reaction_catalog,
)


def _model_for_rule(rule: str, gr_rule: str, genes: list[str]) -> SimpleNamespace:
    return SimpleNamespace(
        rxns=["R1", "R1_complex_formation"],
        rules=[rule, ""],
        gr_rules=[gr_rule, ""],
        genes=genes,
        gene_index={gene_id: index for index, gene_id in enumerate(genes)},
    )


def _enzyme_data() -> tuple[MetabolicEnzymeData, CombinedEnzymeData]:
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
    return metabolic, combined


def _catalog_for_rule(rule: str, gr_rule: str, genes: list[str]):
    metabolic, combined = _enzyme_data()
    return build_gene_enzyme_reaction_catalog(
        _model_for_rule(rule, gr_rule, genes),
        metabolic,
        combined,
    )


def test_isoenzyme_or_rule_is_not_gene_level_executable() -> None:
    catalog = _catalog_for_rule("x(1) | x(2)", "G1 or G2", ["G1", "G2"])

    assert catalog.mappings
    assert {
        mapping.execution_status for mapping in catalog.mappings
    } == {OEExecutionStatus.ISOENZYME_AMBIGUOUS}
    assert all(
        mapping.execution_status is not OEExecutionStatus.GENE_LEVEL_EXECUTABLE
        for mapping in catalog.mappings
    )


def test_complex_and_rule_is_not_gene_level_executable() -> None:
    catalog = _catalog_for_rule("x(1) & x(2)", "G1 and G2", ["G1", "G2"])

    assert catalog.mappings
    assert {
        mapping.execution_status for mapping in catalog.mappings
    } == {OEExecutionStatus.COMPLEX_LIMITED}
    assert all(
        mapping.execution_status is not OEExecutionStatus.GENE_LEVEL_EXECUTABLE
        for mapping in catalog.mappings
    )


def test_mixed_gpr_rule_remains_partial_mapping() -> None:
    catalog = _catalog_for_rule(
        "x(1) | (x(2) & x(3))",
        "G1 or (G2 and G3)",
        ["G1", "G2", "G3"],
    )

    assert catalog.mappings
    assert {
        mapping.execution_status for mapping in catalog.mappings
    } == {OEExecutionStatus.PARTIAL_MAPPING}


def test_external_only_gene_absent_from_current_model_stays_external_only() -> None:
    metabolic, combined = _enzyme_data()
    catalog = build_gene_enzyme_reaction_catalog(
        _model_for_rule("x(1)", "G1", ["G1"]),
        metabolic,
        combined,
        external_evidence=(
            {
                "gene_id": "G_EXTERNAL",
                "reaction_id": "R_EXTERNAL",
                "enzyme_id": "E_EXTERNAL",
                "gpr_role": "single_gene",
                "source_type": "external_pichia_model",
                "source_ref": "external-model/v1",
            },
        ),
    )

    external = next(
        mapping for mapping in catalog.mappings if mapping.gene_id == "G_EXTERNAL"
    )
    assert external.execution_status is OEExecutionStatus.EXTERNAL_EVIDENCE_ONLY
    assert external.model_fingerprint == ""
    assert "current_model_gene_enzyme_reaction_mapping" in external.missing_information


def test_matching_external_evidence_adds_traceability_without_upgrading_status() -> None:
    metabolic, combined = _enzyme_data()
    model = _model_for_rule("x(1) | x(2)", "G1 or G2", ["G1", "G2"])
    baseline = build_gene_enzyme_reaction_catalog(model, metabolic, combined)
    enriched = build_gene_enzyme_reaction_catalog(
        model,
        metabolic,
        combined,
        external_evidence=(
            {
                "gene_id": "G1",
                "reaction_id": "R1",
                "source_type": "pichia_literature",
                "source_ref": "literature/example-1",
            },
        ),
    )

    baseline_g1 = next(mapping for mapping in baseline.mappings if mapping.gene_id == "G1")
    enriched_g1 = next(mapping for mapping in enriched.mappings if mapping.gene_id == "G1")
    assert enriched_g1.execution_status is baseline_g1.execution_status
    assert enriched_g1.mapping_id == baseline_g1.mapping_id
    assert enriched_g1.warnings == (
        "External traceability only: literature/example-1",
    )


def test_model_fingerprint_and_mapping_ids_are_stable_for_identical_inputs() -> None:
    metabolic, combined = _enzyme_data()
    model = _model_for_rule("x(1)", "G1", ["G1"])

    first = build_gene_enzyme_reaction_catalog(model, metabolic, combined)
    second = build_gene_enzyme_reaction_catalog(model, metabolic, combined)

    assert first.model_fingerprint == second.model_fingerprint
    assert first.model_fingerprint
    assert tuple(mapping.mapping_id for mapping in first.mappings) == tuple(
        mapping.mapping_id for mapping in second.mappings
    )
