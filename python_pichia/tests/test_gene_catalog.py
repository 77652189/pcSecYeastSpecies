from __future__ import annotations

from pcsec_pichia.services.gene_catalog import (
    SECRETION_GENE_CATALOG,
    get_ko_reactions_for_selection,
    get_oe_reactions_for_selection,
)


def test_every_entry_is_executable_in_the_direction_its_intervention_declares() -> None:
    """Regression test for a curation bug found in review: several catalog entries declared
    intervention="KO" but stored their reaction id in oe_reaction_id (or vice versa for
    intervention="OE"), so get_ko_reactions_for_selection/get_oe_reactions_for_selection -
    which route purely by which field is populated, not by the intervention label - silently
    returned nothing for the declared direction and surfaced the entry under the wrong one
    instead. gene_id (used for GPR-based KO) counts as KO-capable for "both" entries that
    rely on gene-level resolution rather than a direct reaction id.
    """
    mismatched_ko: list[str] = []
    mismatched_oe: list[str] = []
    for entry in SECRETION_GENE_CATALOG:
        if entry.intervention in ("KO", "both") and not (entry.gene_id or entry.ko_reaction_id):
            mismatched_ko.append(entry.common_name)
        if entry.intervention in ("OE", "both") and not (entry.gene_id or entry.oe_reaction_id):
            mismatched_oe.append(entry.common_name)

    assert mismatched_ko == []
    assert mismatched_oe == []


def test_vacuolar_sorting_ko_entry_routes_to_ko_selection_not_oe() -> None:
    values = get_ko_reactions_for_selection(["AP-3 衔接蛋白复合体"])

    assert values == ["sec_Apl6p_Aps3p_Apm3p_Apl5p_Vam3p_Clc1p_Chc1p_Arf1p_Swa2p_Vps1p_complex_formation"]
    assert get_oe_reactions_for_selection(["AP-3 衔接蛋白复合体"]) == []


def test_ubc6_ubc7_and_doa10_route_to_oe_selection() -> None:
    assert get_oe_reactions_for_selection(["UBC6/UBC7"]) == [
        "sec_Ubc6p_Ubc7p_Yos9p_Hrd1p_Hrd3p_Der1p_Usa1p_complex_formation"
    ]
    assert get_oe_reactions_for_selection(["DOA10"]) == ["sec_Ubc6p_Ubc7p_Doa10p_complex_formation"]
