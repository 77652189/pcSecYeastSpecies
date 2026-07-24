from __future__ import annotations

from pcsec_pichia.strain_modifications import StrainModifications, apply_strain_modifications


class _FakeModel:
    def __init__(self, reaction_index: dict[str, int], ub: list[float]) -> None:
        self.reaction_index = reaction_index
        self.ub = ub
        self.bound_changes: dict[str, tuple[float | None, float | None]] | None = None

    def with_bounds(self, changes: dict[str, tuple[float | None, float | None]]) -> "_FakeModel":
        new = _FakeModel(self.reaction_index, self.ub)
        new.bound_changes = dict(changes)
        return new


class _FakeSecretory:
    def __init__(self, calls: tuple[tuple[str, float], ...] = ()) -> None:
        self.calls = calls

    def with_complex_kcat_multiplier(self, complex_id: str, factor: float) -> "_FakeSecretory":
        return _FakeSecretory((*self.calls, (complex_id, factor)))


class _FakeCombined:
    def __init__(self, calls: tuple[tuple[str, float], ...] = ()) -> None:
        self.calls = calls

    def with_enzyme_kcat_multiplier(self, enzyme_id: str, factor: float) -> "_FakeCombined":
        return _FakeCombined((*self.calls, (enzyme_id, factor)))


def _model() -> _FakeModel:
    return _FakeModel(
        reaction_index={
            "sec_PDI_complex_formation": 0,
            "Mach_translocation_formation": 1,
            "R_plain": 2,
            "R_ko": 3,
        },
        ub=[10.0, 10.0, 4.0, 7.0],
    )


def test_empty_spec_is_noop_returns_inputs_unchanged() -> None:
    model, sec, comb = _model(), _FakeSecretory(), _FakeCombined()
    out_model, out_sec, out_comb, applied, warnings = apply_strain_modifications(
        model, sec, comb, StrainModifications()
    )
    assert out_model is model and out_sec is sec and out_comb is comb
    assert applied == () and warnings == ()


def test_ko_reaction_pins_bounds_to_zero() -> None:
    out_model, _, _, applied, warnings = apply_strain_modifications(
        _model(), _FakeSecretory(), _FakeCombined(), StrainModifications(ko_reaction_ids=("R_ko",))
    )
    assert out_model.bound_changes == {"R_ko": (0.0, 0.0)}
    assert applied == ({"reaction_id": "R_ko", "kind": "KO", "capacity_basis": "reaction_bounds_zero"},)
    assert warnings == ()


def test_oe_sec_complex_multiplies_secretory_kcat_by_factor() -> None:
    _, out_sec, _, applied, _ = apply_strain_modifications(
        _model(), _FakeSecretory(), _FakeCombined(),
        StrainModifications(oe_reaction_ids=("sec_PDI_complex_formation",), oe_factor=3.0),
    )
    # complex_id is the reaction id minus the _formation suffix (matches run_pcsec_oe_screen).
    assert out_sec.calls == (("sec_PDI_complex", 3.0),)
    assert applied[0]["capacity_basis"] == "secretory_complex_kcat_multiplier"
    assert applied[0]["factor"] == 3.0


def test_oe_machine_complex_multiplies_combined_kcat() -> None:
    _, _, out_comb, applied, _ = apply_strain_modifications(
        _model(), _FakeSecretory(), _FakeCombined(),
        StrainModifications(oe_reaction_ids=("Mach_translocation_formation",), oe_factor=2.0),
    )
    assert out_comb.calls == (("Mach_translocation", 2.0),)
    assert applied[0]["capacity_basis"] == "machine_complex_kcat_multiplier"


def test_oe_plain_reaction_scales_upper_bound() -> None:
    out_model, out_sec, out_comb, applied, _ = apply_strain_modifications(
        _model(), _FakeSecretory(), _FakeCombined(),
        StrainModifications(oe_reaction_ids=("R_plain",), oe_factor=2.5),
    )
    assert out_model.bound_changes == {"R_plain": (None, 10.0)}  # ub 4.0 * 2.5
    assert out_sec.calls == () and out_comb.calls == ()
    assert applied[0]["capacity_basis"] == "reaction_upper_bound"


def test_unknown_reaction_ids_are_warned_not_faked() -> None:
    _, _, _, applied, warnings = apply_strain_modifications(
        _model(), _FakeSecretory(), _FakeCombined(),
        StrainModifications(ko_reaction_ids=("nope_ko",), oe_reaction_ids=("nope_oe",)),
    )
    assert applied == ()
    assert any("nope_ko" in w for w in warnings)
    assert any("nope_oe" in w for w in warnings)


def test_ko_wins_when_same_reaction_requested_for_ko_and_oe() -> None:
    out_model, out_sec, _, applied, warnings = apply_strain_modifications(
        _model(), _FakeSecretory(), _FakeCombined(),
        StrainModifications(ko_reaction_ids=("sec_PDI_complex_formation",), oe_reaction_ids=("sec_PDI_complex_formation",)),
    )
    assert out_model.bound_changes == {"sec_PDI_complex_formation": (0.0, 0.0)}
    assert out_sec.calls == ()  # OE dropped
    assert [entry["kind"] for entry in applied] == ["KO"]
    assert any("both KO and OE" in w for w in warnings)


def test_neutral_oe_factor_skips_oe_with_warning_but_keeps_ko() -> None:
    out_model, out_sec, _, applied, warnings = apply_strain_modifications(
        _model(), _FakeSecretory(), _FakeCombined(),
        StrainModifications(ko_reaction_ids=("R_ko",), oe_reaction_ids=("sec_PDI_complex_formation",), oe_factor=1.0),
    )
    assert out_model.bound_changes == {"R_ko": (0.0, 0.0)}
    assert out_sec.calls == ()
    assert [entry["kind"] for entry in applied] == ["KO"]
    assert any("no capacity effect" in w for w in warnings)
