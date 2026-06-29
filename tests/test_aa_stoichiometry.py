from __future__ import annotations

import pytest

from pcsec_pichia.adapters.aa_stoichiometry import AminoAcidStoichiometry
from pcsec_pichia.core.paths import ProjectPaths


def test_amino_acid_stoichiometry_loads_translation_and_degradation_tables() -> None:
    tables = AminoAcidStoichiometry.from_workbook(ProjectPaths.discover().pichia_aa_id_xlsx)

    assert len(tables.amino_acids) == 20
    assert tables.translation_substrates["M"] == "mettrna[c]"
    assert tables.translation_products["M"] == "trnamet[c]"
    assert tables.degradation_products["M"] == "met_L[c]"
    assert tables.energy_substrates == ("h2o[c]", "atp[c]", "gtp[c]")
    assert tables.energy_products == ("h[c]", "adp[c]", "gdp[c]", "pi[c]")


def test_translation_stoichiometry_matches_matlab_countaa_energy_rules() -> None:
    tables = AminoAcidStoichiometry.from_workbook(ProjectPaths.discover().pichia_aa_id_xlsx)

    stoichiometry = tables.translation_stoichiometry("TEST", "MST")

    assert stoichiometry["mettrna[c]"] == pytest.approx(-2.0)
    assert stoichiometry["trnamet[c]"] == pytest.approx(2.0)
    assert stoichiometry["sertrna[c]"] == pytest.approx(-1.0)
    assert stoichiometry["trnaser[c]"] == pytest.approx(1.0)
    assert stoichiometry["h2o[c]"] == pytest.approx(-18.0)
    assert stoichiometry["atp[c]"] == pytest.approx(-9.0)
    assert stoichiometry["gtp[c]"] == pytest.approx(-9.0)
    assert stoichiometry["h[c]"] == pytest.approx(18.0)
    assert stoichiometry["adp[c]"] == pytest.approx(9.0)
    assert stoichiometry["gdp[c]"] == pytest.approx(9.0)
    assert stoichiometry["pi[c]"] == pytest.approx(18.0)
    assert stoichiometry["TEST_peptide[c]"] == pytest.approx(1.0)


def test_degradation_stoichiometry_matches_matlab_countaa_deg_rules() -> None:
    tables = AminoAcidStoichiometry.from_workbook(ProjectPaths.discover().pichia_aa_id_xlsx)

    signal_peptide = tables.signal_peptide_degradation_stoichiometry("TEST", "MST")
    subunit = tables.subunit_degradation_stoichiometry("TEST", "MST")

    assert signal_peptide["met_L[c]"] == pytest.approx(1.0)
    assert signal_peptide["h2o[c]"] == pytest.approx(-3.0)
    assert signal_peptide["atp[c]"] == pytest.approx(-3.0)
    assert signal_peptide["h[c]"] == pytest.approx(3.0)
    assert signal_peptide["adp[c]"] == pytest.approx(3.0)
    assert signal_peptide["pi[c]"] == pytest.approx(3.0)
    assert signal_peptide["TEST_sp[c]"] == pytest.approx(-1.0)

    assert subunit["met_L[c]"] == pytest.approx(2.0)
    assert subunit["TEST_subunit[c]"] == pytest.approx(-1.0)
