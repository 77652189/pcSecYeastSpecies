from __future__ import annotations

from app.services.signal_peptide_rules import score_signal_peptide


def test_rule_screening_accepts_typical_secretory_signal_peptide() -> None:
    result = score_signal_peptide("MKALLLALLALAAASAGA")

    assert result.passed is True
    assert result.score >= 65
    assert result.h_region_max_hydrophobicity >= 1.8
    assert result.c_region_small_neutral_rule is True


def test_rule_screening_rejects_low_complexity_non_hydrophobic_sequence() -> None:
    result = score_signal_peptide("MNNNNNNNNNNNNNNNNN")

    assert result.passed is False
    assert "疏水核心偏弱" in "；".join(result.risks)
    assert "较长重复残基" in "；".join(result.risks)
