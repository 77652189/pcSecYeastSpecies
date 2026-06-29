from __future__ import annotations

from pathlib import Path

import pytest

from pcsec_pichia.adapters.lp_alignment import compare_lp_files, compare_lp_math, write_lp_math_alignment_json
from pcsec_pichia.core.paths import ProjectPaths


def _write_minimal_lp(
    path: Path,
    objective: str = "X3",
    c1_rhs: str = "0",
    cm1_coefficient: str = "10.000000000000000",
    x3_upper: str = "1.000000",
) -> None:
    path.write_text(
        "\n".join(
            [
                "Maximize",
                f"obj: {objective}",
                "Subject To",
                f"C1: 1.000000000000000 X1 -1.000000000000000 X2 = {c1_rhs}",
                f"CM1: X2 - {cm1_coefficient} X3 = 0",
                "Bounds",
                "0.000000 <= X1 <= +infinity",
                "0.000000 <= X2 <= +infinity",
                f"0.000000 <= X3 <= {x3_upper}",
                "End",
            ]
        ),
        encoding="utf-8",
    )


def test_compare_lp_files_reports_matching_structure(tmp_path: Path) -> None:
    reference = tmp_path / "reference.lp"
    candidate = tmp_path / "candidate.lp"
    _write_minimal_lp(reference)
    _write_minimal_lp(candidate)

    comparison = compare_lp_files(reference, candidate)

    assert comparison.is_match is True
    assert comparison.differing_fields == {}
    assert "objective_variable" in comparison.matching_fields
    assert comparison.to_dict()["is_match"] is True


def test_compare_lp_files_reports_objective_mismatch(tmp_path: Path) -> None:
    reference = tmp_path / "reference.lp"
    candidate = tmp_path / "candidate.lp"
    _write_minimal_lp(reference, objective="X3")
    _write_minimal_lp(candidate, objective="X2")

    comparison = compare_lp_files(reference, candidate)

    assert comparison.is_match is False
    assert comparison.differing_fields["objective_variable"].reference == "X3"
    assert comparison.differing_fields["objective_variable"].candidate == "X2"


def test_compare_lp_math_treats_equivalent_rhs_formatting_as_match(tmp_path: Path) -> None:
    reference = tmp_path / "reference.lp"
    candidate = tmp_path / "candidate.lp"
    _write_minimal_lp(reference, c1_rhs="0")
    _write_minimal_lp(candidate, c1_rhs="0.000000000000000")

    comparison = compare_lp_math(reference, candidate)

    assert comparison.is_match is True
    assert comparison.constraint_difference_count == 0
    assert comparison.bound_difference_count == 0
    assert comparison.to_dict()["is_match"] is True


def test_compare_lp_math_reports_rhs_coefficient_and_bound_differences(tmp_path: Path) -> None:
    reference = tmp_path / "reference.lp"
    candidate = tmp_path / "candidate.lp"
    _write_minimal_lp(reference)
    _write_minimal_lp(
        candidate,
        c1_rhs="0.1",
        cm1_coefficient="11.000000000000000",
        x3_upper="2.000000",
    )

    comparison = compare_lp_math(reference, candidate)

    assert comparison.is_match is False
    assert comparison.constraint_difference_count == 2
    assert comparison.bound_difference_count == 1
    differences_by_name = {difference.name: difference for difference in comparison.constraint_differences}
    assert differences_by_name["C1"].kind == "structure_or_rhs"
    assert differences_by_name["C1"].candidate_rhs == pytest.approx(0.1)
    assert differences_by_name["CM1"].kind == "coefficient"
    assert differences_by_name["CM1"].variable == "X3"
    assert differences_by_name["CM1"].absolute_difference == pytest.approx(1.0)
    assert comparison.bound_differences[0].variable == "X3"
    assert comparison.bound_differences[0].candidate_upper == pytest.approx(2.0)


def test_write_lp_math_alignment_json(tmp_path: Path) -> None:
    reference = tmp_path / "reference.lp"
    candidate = tmp_path / "candidate.lp"
    _write_minimal_lp(reference)
    _write_minimal_lp(candidate)
    comparison = compare_lp_math(reference, candidate)

    output = write_lp_math_alignment_json(comparison, tmp_path / "math_alignment.json")

    assert output.exists()
    assert '"is_match": true' in output.read_text(encoding="utf-8")


def test_compare_existing_opn_target_lp_math_when_available() -> None:
    paths = ProjectPaths.discover()
    matlab_lp = (
        paths.opn_run_dir
        / "Simulation_dilutionOPN_ALPHA_FULL_PROJECT_mu0p1_media4_ratio1em08_misfolddefault_noMisfoldEq_noRiboEq_PP.lp"
    )
    python_lp = (
        paths.local_runs_dir
        / "pichia_python"
        / "python_OPN_ALPHA_FULL_PROJECT_mu0p10_media4_ratio1em08_noMisfoldEq_noRiboEq.lp"
    )
    if not matlab_lp.exists() or not python_lp.exists():
        pytest.skip("Local MATLAB/Python OPN target LP alignment artifacts have not been generated on this machine.")

    comparison = compare_lp_math(matlab_lp, python_lp)

    assert comparison.is_match is True
    assert comparison.constraint_difference_count == 0
    assert comparison.bound_difference_count == 0
