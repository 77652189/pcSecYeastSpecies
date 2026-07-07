from __future__ import annotations

from types import SimpleNamespace
from typing import Mapping

import pytest

from pcsec_pichia.analysis.shadow_lp import (
    REFERENCE_LAYER_ORDER,
    build_metabolic_coupling_block,
    build_shadow_constraint_blocks,
    prepare_builtin_shadow_target,
)


EXPECTED_LAYER_COUNTS = {
    "metabolic_coupling": 2732,
    "secretory_coupling": 58,
    "protein_mass": 2,
    "proteasome": 1,
    "ribosome_assembly": 1,
    "ribosome_translation": 0,
    "misfolding": 0,
    "mitochondrial": 1,
}


@pytest.fixture(scope="module")
def prepared_targets() -> Mapping[str, object]:
    return {
        "hLF": prepare_builtin_shadow_target("hLF"),
        "OPN_ALPHA_FULL_PROJECT": prepare_builtin_shadow_target("OPN_ALPHA_FULL_PROJECT"),
    }


@pytest.fixture(scope="module")
def target_blocks(prepared_targets: Mapping[str, object]) -> Mapping[str, dict[str, object]]:
    return {
        target_id: {block.layer_id: block for block in build_shadow_constraint_blocks(prep)}
        for target_id, prep in prepared_targets.items()
    }


def test_builders_emit_reference_layer_counts_for_hlf_and_opn(
    target_blocks: Mapping[str, dict[str, object]],
) -> None:
    for blocks in target_blocks.values():
        assert tuple(blocks) == REFERENCE_LAYER_ORDER
        assert {layer_id: block.counts[layer_id] for layer_id, block in blocks.items()} == EXPECTED_LAYER_COUNTS

        for block in blocks.values():
            for constraint in block.constraints:
                assert constraint.source
                assert constraint.metadata
                assert all(isinstance(reaction_id, str) for reaction_id in constraint.terms)
            if block.layer_id not in {"ribosome_translation", "misfolding"}:
                assert block.missing_mapping_count == 0
                assert block.warnings == ()

    hlf_counts = {layer_id: block.counts[layer_id] for layer_id, block in target_blocks["hLF"].items()}
    opn_counts = {
        layer_id: block.counts[layer_id]
        for layer_id, block in target_blocks["OPN_ALPHA_FULL_PROJECT"].items()
    }
    assert hlf_counts == opn_counts


def test_representative_coefficients_and_rhs_match_full_ladder_golden_payload(
    target_blocks: Mapping[str, dict[str, object]],
) -> None:
    blocks = target_blocks["hLF"]

    metabolic = blocks["metabolic_coupling"].constraints[0]
    assert metabolic.terms["FACOAE140_no_1_fwd"] == pytest.approx(1.0)
    assert metabolic.terms["FACOAE140_no_1_fwd_complex_formation"] == pytest.approx(-134866.8)
    assert metabolic.rhs == pytest.approx(0.0)

    secretory = blocks["secretory_coupling"].constraints[0]
    secretory_formation_id = (
        "sec_Apl6p_Aps3p_Apm3p_Apl5p_Vam3p_Clc1p_Chc1p_Arf1p_Swa2p_Vps1p_complex_formation"
    )
    assert secretory.terms[secretory_formation_id] == pytest.approx(-2151.1026453250893)
    assert secretory.rhs == pytest.approx(0.0)

    protein_mass = blocks["protein_mass"].constraints[0]
    assert protein_mass.rhs == pytest.approx(0.032552119, rel=1e-8)
    assert protein_mass.terms["dilute_dummy"] == pytest.approx(40.0)

    er_mass = blocks["protein_mass"].constraints[1]
    assert er_mass.rhs == pytest.approx(0.004)
    assert er_mass.terms == {"dilute_dummyER": 40.0}

    proteasome = blocks["proteasome"].constraints[0]
    assert proteasome.terms["Mach_proteasome_complex_formation"] == pytest.approx(-600.0)

    ribosome = blocks["ribosome_assembly"].constraints[0]
    assert ribosome.terms["Mach_Ribosome_complex_formation"] == pytest.approx(1.0)
    assert ribosome.terms["Mach_Ribosome_Assembly_Factors_complex_formation"] == pytest.approx(-1_200_000.0)

    mitochondrial = blocks["mitochondrial"].constraints[0]
    assert mitochondrial.sense == "le"
    assert mitochondrial.rhs == pytest.approx(0.005)


def test_disabled_reference_layers_emit_zero_count_blocks(
    target_blocks: Mapping[str, dict[str, object]],
) -> None:
    blocks = target_blocks["hLF"]

    assert blocks["ribosome_translation"].constraints == ()
    assert blocks["ribosome_translation"].counts == {"ribosome_translation": 0}
    assert blocks["ribosome_translation"].warnings
    assert blocks["misfolding"].constraints == ()
    assert blocks["misfolding"].counts == {"misfolding": 0}
    assert blocks["misfolding"].warnings


def test_missing_mapping_is_reported_as_warning_not_silently_dropped() -> None:
    prep = SimpleNamespace(
        target_id="fake",
        fixed_model=SimpleNamespace(rxns=("R1",), reaction_index={"R1": 0}),
        metabolic=SimpleNamespace(enzymes=("missing_complex",), kcat=(10.0,)),
    )

    block = build_metabolic_coupling_block(prep)

    assert block.constraints == ()
    assert block.counts == {"metabolic_coupling": 0}
    assert block.missing_mapping_count == 1
    assert block.warnings == (
        "metabolic_coupling missing mapping for enzyme missing_complex: missing, missing_complex_formation",
    )
