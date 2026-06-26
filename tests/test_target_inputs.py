from __future__ import annotations

from app.core.opn_inputs import BuiltinOpnInputProvider
from pcsec_pichia.core.target_inputs import (
    LeaderCandidateInput,
    TargetInputSet,
    TargetProteinInput,
    TargetRegistry,
    build_pcsec_target_row,
    draft_hlf_input_set,
    target_input_set_from_opn_input_set,
)


def test_opn_inputs_convert_to_generic_target_schema() -> None:
    opn_input_set = BuiltinOpnInputProvider().load_input_set()
    generic = target_input_set_from_opn_input_set(opn_input_set)

    assert generic.ready_for_model is True
    assert generic.target.target_id == "OPN"
    assert len(generic.leaders) == len(opn_input_set.candidates)

    row = build_pcsec_target_row(generic.target, generic.leaders[0])
    assert row["Protein name"] == "OPN_ALPHA_FULL_PROJECT"
    assert str(row["sequence"]).endswith(generic.target.mature_sequence)
    assert int(row["Length"]) == len(str(row["sequence"]))


def test_hlf_draft_is_readable_but_not_model_ready() -> None:
    hlf = draft_hlf_input_set()

    assert hlf.target.target_id == "hLF"
    assert hlf.ready_for_model is False
    assert "参数待确认" in hlf.target.readiness_message()


def test_target_registry_accepts_future_targets_without_engine_changes() -> None:
    registry = TargetRegistry()
    registry.register(draft_hlf_input_set())
    registry.register(
        TargetInputSet(
            target=TargetProteinInput(
                target_id="TEST_PROTEIN",
                protein_name="test protein",
                abbreviation="TEST",
                mature_sequence="MSTNPKPQR",
                parameter_status="ready_for_model",
            ),
            leaders=[
                LeaderCandidateInput(
                    candidate_id="TEST_LEADER",
                    leader_sequence="MKFAISTLLIILQAAAVFAA",
                    signal_peptide_sequence="MKFAISTLLIILQAAAVFAA",
                )
            ],
        )
    )

    targets = {target.target_id for target in registry.list_targets()}
    assert {"hLF", "TEST_PROTEIN"}.issubset(targets)
    assert registry.get("TEST_PROTEIN").ready_for_model is True
