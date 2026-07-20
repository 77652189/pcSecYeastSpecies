from __future__ import annotations

from app.ui.views.simulation import _prefill_field_values


def test_gene_candidate_routes_to_gene_inputs() -> None:
    values = _prefill_field_values("PAS_chr1-4_0047", "KO", "gene")

    assert values["pichia_draft_ko_genes"] == "PAS_chr1-4_0047"
    assert values["pichia_draft_ko_reactions"] == ""
    assert values["pichia_draft_oe_genes"] == ""
    assert values["pichia_draft_oe_reactions"] == ""


def test_catalog_reaction_candidate_routes_to_reaction_inputs() -> None:
    values = _prefill_field_values("sec_PDI1_ERV2_Ero1p_complex_formation", "OE", "catalog_reaction")

    assert values["pichia_draft_oe_reactions"] == "sec_PDI1_ERV2_Ero1p_complex_formation"
    assert values["pichia_draft_oe_genes"] == ""
    assert values["pichia_draft_ko_genes"] == ""
    assert values["pichia_draft_ko_reactions"] == ""


def test_complex_oe_hypothesis_candidate_routes_to_reaction_inputs_not_gene_inputs() -> None:
    """Regression test: candidate_kind != "catalog_reaction" used to be read as "is a gene",
    which would have sent this reaction id into the gene-ID box (and let it silently fail
    GPR resolution) instead of the reaction-ID box where it actually belongs.
    """
    values = _prefill_field_values("ATPS3m_no_1_fwd", "OE", "complex_oe_hypothesis")

    assert values["pichia_draft_oe_reactions"] == "ATPS3m_no_1_fwd"
    assert values["pichia_draft_oe_genes"] == ""


def test_unknown_future_candidate_kind_fails_safe_to_reaction_routing() -> None:
    values = _prefill_field_values("some_id", "KO", "some_new_kind_nobody_added_yet")

    assert values["pichia_draft_ko_reactions"] == "some_id"
    assert values["pichia_draft_ko_genes"] == ""
