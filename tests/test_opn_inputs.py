from __future__ import annotations

from app.core.opn_inputs import CsvOpnCandidateInputProvider


def test_csv_opn_candidate_input_provider_loads_candidates() -> None:
    content = (
        "candidate_id,leader_sequence,signal_peptide_sequence,category,processing_route,"
        "source_note,rationale,caution\n"
        "SIGSCOUT_REP_001,mkfaistlliilqaaavfaa,MKFAISTLLIILQAAAVFAA,"
        "pichia_native_signal,signal peptidase only,SigScout representative,"
        "candidate selected by representative workflow,needs wet-lab confirmation\n"
    )

    result = CsvOpnCandidateInputProvider(content=content, source_name="sigscout csv").load_input_set()

    assert result.errors == []
    assert result.source_name == "sigscout csv"
    assert result.target.protein_name == "OPN"
    assert result.candidates[0].candidate_id == "SIGSCOUT_REP_001"
    assert result.candidates[0].leader_sequence == "MKFAISTLLIILQAAAVFAA"


def test_csv_opn_candidate_input_provider_reports_missing_columns() -> None:
    result = CsvOpnCandidateInputProvider(content="candidate_id,leader_sequence\nBROKEN,ABC\n").load_input_set()

    assert result.candidates == []
    assert "Missing required OPN candidate columns" in result.errors[0]
