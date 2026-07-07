from __future__ import annotations

from typing import Mapping

import pytest

from pcsec_pichia.analysis.shadow_lp import (
    FORMAL_SHADOW_LADDER_ORDER,
    attach_reference_validation,
    run_shadow_hardcode_audit,
    run_shadow_ladder,
)


@pytest.fixture(scope="module")
def validated_ladders() -> Mapping[str, object]:
    payload = {}
    for target_id in ("hLF", "OPN_ALPHA_FULL_PROJECT"):
        ladder = run_shadow_ladder(target_id)
        payload[target_id] = attach_reference_validation(ladder)
    return payload


def test_shadow_ladder_layer_order_and_skipped_layers(validated_ladders: Mapping[str, object]) -> None:
    for ladder in validated_ladders.values():
        assert tuple(layer.layer_id for layer in ladder.layers) == FORMAL_SHADOW_LADDER_ORDER
        assert ladder.final_layer_id == "mitochondrial"
        assert ladder.final_layer.layer_id == "mitochondrial"
        assert set(ladder.layers[0].key_fluxes) == {
            "BIOMASS",
            "Ex_glc_D",
            "Ex_o2",
            ladder.exchange_reaction_id,
        }

        skipped = {layer.layer_id: layer for layer in ladder.layers if layer.status == "skipped"}
        assert set(skipped) == {"ribosome_translation", "misfolding"}
        assert skipped["ribosome_translation"].skipped_layers["ribosome_translation"]
        assert skipped["misfolding"].skipped_layers["misfolding"]


def test_shadow_ladder_reference_validation_matches_pcsec_reference(
    validated_ladders: Mapping[str, object],
) -> None:
    for ladder in validated_ladders.values():
        validation = ladder.reference_validation
        assert validation is not None
        assert validation["final_alignment_status"] == "aligned"
        assert validation["objective_rel_diff"] <= 1e-4
        assert validation["constraint_count_diff"] == 0
        assert ladder.final_layer.success is True
        assert ladder.final_layer.objective == pytest.approx(validation["reference_objective"], rel=1e-4)


def test_shadow_ladder_forbidden_solver_audit_passes() -> None:
    audit = run_shadow_hardcode_audit()

    assert audit.passed
    assert audit.production_shadow_path_has_no_forbidden_reference_solver_call
    assert audit.reference_values_only_used_after_solve_for_comparison
