from __future__ import annotations

from pcsec_pichia.homology.review_rules import (
    ALIAS_CONFIRMED_BY_RBH,
    LOW_IDENTITY_REVIEW_REQUIRED,
    MODEL_READY_RBH_HIGH_CONFIDENCE,
    NO_RECIPROCAL_HIT,
    PARALOG_RISK_REVIEW_REQUIRED,
    RBH_NOT_IN_MODEL,
    SEQUENCE_NAME_CONFLICT,
    classify_homology_review_status,
    classify_name_consistency,
)


def test_high_confidence_requires_rbh_thresholds_and_model_gene() -> None:
    status, warnings = classify_homology_review_status(
        is_rbh=True,
        identity_pct=70,
        query_coverage=95,
        subject_coverage=90,
        in_model_gene_index=True,
    )

    assert status == MODEL_READY_RBH_HIGH_CONFIDENCE
    assert warnings == ()


def test_rbh_not_in_model_is_homology_evidence_only() -> None:
    status, warnings = classify_homology_review_status(
        is_rbh=True,
        identity_pct=70,
        query_coverage=95,
        subject_coverage=90,
        in_model_gene_index=False,
    )

    assert status == RBH_NOT_IN_MODEL
    assert warnings == ("RBH candidate is not present in current Pichia GEM gene_index",)


def test_low_identity_and_non_rbh_cannot_be_high_confidence() -> None:
    low_identity, _ = classify_homology_review_status(
        is_rbh=True,
        identity_pct=24,
        query_coverage=95,
        subject_coverage=90,
        in_model_gene_index=True,
    )
    non_rbh, _ = classify_homology_review_status(
        is_rbh=False,
        identity_pct=90,
        query_coverage=95,
        subject_coverage=90,
        in_model_gene_index=True,
    )

    assert low_identity == LOW_IDENTITY_REVIEW_REQUIRED
    assert non_rbh == NO_RECIPROCAL_HIT


def test_paralog_risk_has_review_status() -> None:
    status, warnings = classify_homology_review_status(
        is_rbh=True,
        identity_pct=70,
        query_coverage=95,
        subject_coverage=90,
        in_model_gene_index=True,
        paralog_count=2,
    )

    assert status == PARALOG_RISK_REVIEW_REQUIRED
    assert "paralog risk" in warnings[0]


def test_name_consistency_flags_aliases_and_conflicts() -> None:
    assert (
        classify_name_consistency(
            internal_common_name="DOA10",
            external_gene_name="SSM4",
            external_aliases=("DOA10",),
            is_rbh=True,
        )
        == ALIAS_CONFIRMED_BY_RBH
    )
    assert (
        classify_name_consistency(
            internal_common_name="SSA1",
            external_gene_name="SSA2",
            external_aliases=(),
            is_rbh=True,
        )
        == SEQUENCE_NAME_CONFLICT
    )
