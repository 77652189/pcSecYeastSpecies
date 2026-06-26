from __future__ import annotations

import pytest

from pcsec_pichia.adapters.mat_loader import MatStructLoader
from pcsec_pichia.core.pichia_media import apply_glucose_reference_conditions, set_fixed_growth_rate, set_media_pp
from pcsec_pichia.core.paths import ProjectPaths


@pytest.fixture(scope="module")
def pichia_model():
    return MatStructLoader(ProjectPaths.discover()).load_pcsec_pichia_model()


def _bounds(model, reaction_id: str) -> tuple[float, float]:
    index = model.reaction_index[reaction_id]
    return float(model.lb[index]), float(model.ub[index])


def test_set_media_pp_type_4_matches_matlab_supplement_bounds(pichia_model) -> None:
    configured = set_media_pp(pichia_model, media_type=4)
    model = configured.model

    assert configured.missing_reactions == ()
    assert _bounds(model, "Ex_nh4") == (0.0, 1000.0)
    assert _bounds(model, "Ex_o2") == (0.0, 1000.0)
    assert _bounds(model, "Ex_h2o") == (0.0, 1000.0)
    assert _bounds(model, "Ex_btn") == (-2.0, 1000.0)
    assert _bounds(model, "Ex_arg_L") == (-0.08, 1000.0)
    assert _bounds(model, "Ex_ura") == (-0.08, 1000.0)
    assert _bounds(model, "Ex_ala_L") == (0.0, 1000.0)


def test_set_media_pp_can_open_minimal_exchanges_for_basic_fba_smokes(pichia_model) -> None:
    configured = set_media_pp(pichia_model, media_type=4, open_minimal_exchanges=True)
    model = configured.model

    assert _bounds(model, "Ex_nh4") == (-1000.0, 1000.0)
    assert _bounds(model, "Ex_o2") == (-1000.0, 1000.0)
    assert _bounds(model, "Ex_h2o") == (-1000.0, 1000.0)


def test_glucose_reference_conditions_block_glycerol_and_methanol_growth(pichia_model) -> None:
    configured = apply_glucose_reference_conditions(pichia_model, media_type=2)
    model = configured.model

    assert configured.missing_reactions == ()
    assert _bounds(model, "Ex_nh4") == (0.0, 1000.0)
    assert _bounds(model, "Ex_h2o") == (0.0, 1000.0)
    assert _bounds(model, "Ex_glc_D") == (-1000.0, 1000.0)
    assert _bounds(model, "Ex_o2") == (-1000.0, 1000.0)
    assert _bounds(model, "Ex_glyc") == (0.0, 1000.0)
    assert _bounds(model, "BIOMASS_glyc") == (0.0, 0.0)
    assert _bounds(model, "Ex_meoh") == (0.0, 1000.0)
    assert _bounds(model, "BIOMASS_meoh") == (0.0, 0.0)


def test_set_fixed_growth_rate_sets_biomass_equality(pichia_model) -> None:
    fixed = set_fixed_growth_rate(pichia_model, 0.10)

    assert _bounds(fixed, "BIOMASS") == (0.10, 0.10)
