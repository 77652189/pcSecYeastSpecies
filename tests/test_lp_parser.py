from __future__ import annotations

from pathlib import Path

import pytest

from pcsec_pichia.adapters.lp_parser import parse_lp_file
from pcsec_pichia.core.paths import ProjectPaths


def test_lp_parser_handles_wrapped_constraints(tmp_path: Path) -> None:
    lp = tmp_path / "wrapped.lp"
    lp.write_text(
        "\n".join(
            [
                "Maximize",
                "obj: X3",
                "Subject To",
                "C1: 1.0 X1 -2.0 X2",
                " + 3.0 X3 = 0",
                "CM2: X3 - 4.0 X4 = 0",
                "Cmito: X4 <= 1.0",
                "Bounds",
                "0 <= X1 <= +infinity",
                "-10 <= X2 <= 0",
                "0 <= X3 <= 1000",
                "0 <= X4 <= 1",
                "End",
            ]
        ),
        encoding="utf-8",
    )

    summary = parse_lp_file(lp)

    assert summary.optimization_sense == "Maximize"
    assert summary.objective_name == "obj"
    assert summary.objective_variable == "X3"
    assert summary.constraint_count == 3
    assert summary.bounds_count == 4
    assert summary.distinct_variable_count == 4
    assert summary.constraint_prefix_counts == {"C": 1, "CM": 1, "Cmito": 1}
    assert summary.constraint_sense_counts == {"=": 2, "<=": 1}


def test_lp_parser_summarizes_existing_matlab_reference_lp_when_available() -> None:
    paths = ProjectPaths.discover()
    lp = (
        paths.local_runs_dir
        / "PPA_GLC_ref_smoke"
        / "Simulation_dilutionref_mu0p1_media4_misfolddefault_openMisfoldDilution_noMisfoldEq_noRiboEq_PP.lp"
    )
    if not lp.exists():
        pytest.skip("Local MATLAB reference LP has not been generated on this machine.")

    summary = parse_lp_file(lp)

    assert summary.optimization_sense == "Maximize"
    assert summary.objective_variable == "X545"
    assert summary.constraint_count == 22990
    assert summary.bounds_count == 29026
    assert summary.distinct_variable_count == 29026
    assert summary.max_variable_index == 29026
    assert summary.constraint_prefix_counts == {"C": 20195, "CM": 2794, "Cmito": 1}
    assert summary.constraint_sense_counts == {"=": 22989, "<=": 1}
