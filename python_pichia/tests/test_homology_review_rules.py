from __future__ import annotations

from pcsec_pichia.homology.review_rules import (
    ALIAS_CONFIRMED_BY_RBH,
    EXTERNAL_ALIAS_CONFIRMED,
    EXTERNAL_CONFLICT,
    EXTERNAL_LOCUS_CONFIRMED,
    EXTERNAL_MATCH_CONFIRMED,
    EXTERNAL_REFERENCE_INCOMPLETE,
    LOW_IDENTITY_REVIEW_REQUIRED,
    MODEL_READY_RBH_HIGH_CONFIDENCE,
    NO_RECIPROCAL_HIT,
    PARALOG_RISK_REVIEW_REQUIRED,
    RBH_NOT_IN_MODEL,
    RULE_TRANSFER_LOW_CONFIDENCE,
    RULE_TRANSFER_NOT_SUPPORTED,
    RULE_TRANSFER_READY,
    RULE_TRANSFER_SUPPORTED_NOT_MODEL_OPERABLE,
    RULE_TRANSFER_UNRESOLVED,
    SEQUENCE_NAME_CONFLICT,
    classify_external_name_crosscheck,
    classify_homology_review_status,
    classify_name_consistency,
    classify_rule_transfer_status,
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


def test_rule_transfer_statuses_preserve_model_and_confidence_boundaries() -> None:
    ready, ready_warnings = classify_rule_transfer_status(
        homology_review_status=MODEL_READY_RBH_HIGH_CONFIDENCE,
        is_rbh=True,
        in_model_gene_index=True,
    )
    not_model, _ = classify_rule_transfer_status(
        homology_review_status=RBH_NOT_IN_MODEL,
        is_rbh=True,
        in_model_gene_index=False,
    )
    low_confidence, _ = classify_rule_transfer_status(
        homology_review_status=LOW_IDENTITY_REVIEW_REQUIRED,
        is_rbh=True,
        in_model_gene_index=True,
    )
    unresolved, _ = classify_rule_transfer_status(
        homology_review_status="unresolved_query_symbol",
        is_rbh=False,
        in_model_gene_index=False,
    )
    not_supported, _ = classify_rule_transfer_status(
        homology_review_status=NO_RECIPROCAL_HIT,
        is_rbh=False,
        in_model_gene_index=True,
    )

    assert ready == RULE_TRANSFER_READY
    assert ready_warnings == ()
    assert not_model == RULE_TRANSFER_SUPPORTED_NOT_MODEL_OPERABLE
    assert low_confidence == RULE_TRANSFER_LOW_CONFIDENCE
    assert unresolved == RULE_TRANSFER_UNRESOLVED
    assert not_supported == RULE_TRANSFER_NOT_SUPPORTED


def test_external_name_crosscheck_classifies_name_alias_locus_and_conflict() -> None:
    name_status, name_warnings = classify_external_name_crosscheck(
        internal_common_name="KAR2 / BiP",
        external_gene_name="KAR2",
        external_locus_tag="PAS_chr2-1_0140",
        external_aliases=(),
        reference_gene_name="KAR2",
        reference_locus_tag="PAS_chr2-1_0140",
        reference_aliases=(),
    )
    alias_status, _ = classify_external_name_crosscheck(
        internal_common_name="DOA10",
        external_gene_name="DOA10",
        external_locus_tag="PAS_chr3_0123",
        external_aliases=(),
        reference_gene_name="SSM4",
        reference_locus_tag="PAS_chr3_0123",
        reference_aliases=("DOA10",),
    )
    locus_status, _ = classify_external_name_crosscheck(
        internal_common_name="VPS1",
        external_gene_name="",
        external_locus_tag="PAS_chr1_0440",
        external_aliases=(),
        reference_gene_name="",
        reference_locus_tag="PAS_chr1_0440",
        reference_aliases=(),
    )
    conflict_status, conflict_warnings = classify_external_name_crosscheck(
        internal_common_name="HRD1",
        external_gene_name="HRD1",
        external_locus_tag="PAS_chr4_0156",
        external_aliases=(),
        reference_gene_name="PEP4",
        reference_locus_tag="PAS_chr4_0156",
        reference_aliases=(),
    )
    incomplete_status, incomplete_warnings = classify_external_name_crosscheck(
        internal_common_name="SEC12",
        external_gene_name="",
        external_locus_tag="",
        external_aliases=(),
        reference_gene_name="",
        reference_locus_tag="",
        reference_aliases=(),
    )

    assert name_status == EXTERNAL_MATCH_CONFIRMED
    assert name_warnings == ()
    assert alias_status == EXTERNAL_ALIAS_CONFIRMED
    assert locus_status == EXTERNAL_LOCUS_CONFIRMED
    assert conflict_status == EXTERNAL_CONFLICT
    assert "does not match" in conflict_warnings[0]
    assert incomplete_status == EXTERNAL_REFERENCE_INCOMPLETE
    assert "lacks comparable" in incomplete_warnings[0]
