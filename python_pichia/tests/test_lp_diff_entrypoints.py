from __future__ import annotations

from pathlib import Path

from scipy import sparse

from pcsec_pichia.alignment.lp_diff import diff_lp_files, parse_lp, write_lp_diff_outputs, write_matrix_lp_from_rows


class TinyModel:
    reaction_index = {"obj_rxn": 1}
    lb = [0.0, 0.0]
    ub = [1000.0, 10.0]


def test_parse_lp_reads_objective_constraints_and_bounds(tmp_path: Path) -> None:
    lp_path = tmp_path / "tiny.lp"
    lp_path.write_text(
        "\n".join(
            [
                "Maximize",
                "obj: X2",
                "Subject To",
                "C1: 1 X1 -2 X2 = 0",
                "C2: 3 X2 <= 4",
                "Bounds",
                "0.000000 <= X1 <= +infinity",
                "-1.000000 <= X2 <= 5.000000",
                "End",
            ]
        ),
        encoding="utf-8",
    )

    parsed = parse_lp(lp_path)

    assert parsed.optimization_sense == "Maximize"
    assert parsed.objective.coefficients == {2: 1.0}
    assert parsed.constraints[0].label == "C1"
    assert parsed.constraints[0].coefficients == {1: 1.0, 2: -2.0}
    assert parsed.constraints[1].sense == "<="
    assert parsed.constraints[1].rhs == 4.0
    assert parsed.bounds[1].upper == float("inf")
    assert parsed.bounds[2].lower == -1.0


def test_lp_diff_classifies_labels_coefficients_rhs_and_bounds(tmp_path: Path) -> None:
    matlab_lp = tmp_path / "matlab.lp"
    python_lp = tmp_path / "python.lp"
    matlab_lp.write_text(
        "\n".join(
            [
                "Maximize",
                "obj: X1",
                "Subject To",
                "C1: 1 X1 = 0",
                "CM1: 1 X1 -2 X2 = 0",
                "C3: 1 X3 <= 1",
                "Bounds",
                "0.000000 <= X1 <= +infinity",
                "0.000000 <= X2 <= 10.000000",
                "0.000000 <= X3 <= 3.000000",
                "End",
            ]
        ),
        encoding="utf-8",
    )
    python_lp.write_text(
        "\n".join(
            [
                "Maximize",
                "obj: X2",
                "Subject To",
                "C1: 1 X1 = 0",
                "C2: 1 X1 -3 X2 = 0",
                "C3: 1 X4 <= 2",
                "Bounds",
                "0.000000 <= X1 <= +infinity",
                "0.000000 <= X2 <= 8.000000",
                "0.000000 <= X3 <= 3.000000",
                "End",
            ]
        ),
        encoding="utf-8",
    )

    diff = diff_lp_files(matlab_lp, python_lp)

    assert diff["objective_diff"]["coefficient_difference_count"] == 2
    assert diff["row_diff_summary"]["ordered_label_difference_count"] == 1
    assert diff["row_diff_summary"]["coefficient_difference_count"] == 2
    assert diff["row_diff_summary"]["rhs_difference_count"] == 1
    assert diff["row_diff_summary"]["sparsity_difference_count"] == 1
    assert diff["bound_diff_summary"]["difference_count"] == 1

    outputs = write_lp_diff_outputs(diff, tmp_path / "out")
    assert outputs["summary"].exists()
    assert outputs["top_rows"].exists()
    assert outputs["report"].exists()


def test_write_matrix_lp_from_rows_preserves_explicit_labels_and_order(tmp_path: Path) -> None:
    lp_path = tmp_path / "ordered.lp"

    write_matrix_lp_from_rows(
        model=TinyModel(),
        rows=[
            ("CM1", sparse.csr_matrix([[1.0, -2.0]]), "=", 0.0),
            ("Cmito", sparse.csr_matrix([[0.0, 3.0]]), "<=", 0.5),
        ],
        objective_reaction="obj_rxn",
        path=lp_path,
    )
    parsed = parse_lp(lp_path)

    assert parsed.objective.coefficients == {2: 1.0}
    assert [row.label for row in parsed.constraints] == ["CM1", "Cmito"]
    assert parsed.constraints[1].sense == "<="
    assert parsed.constraints[1].rhs == 0.5
