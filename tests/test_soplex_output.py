from __future__ import annotations

from pathlib import Path

import pytest

from pcsec_pichia.adapters.soplex_output import compare_soplex_objectives, parse_soplex_output, write_soplex_summary_json
from pcsec_pichia.core.paths import ProjectPaths


def test_parse_soplex_output_reads_status_objective_time_and_iterations(tmp_path: Path) -> None:
    output = tmp_path / "solve.lp.out"
    output.write_text(
        "\n".join(
            [
                "SoPlex status       : problem is solved [optimal]",
                "Solving time (sec)  : 27.04",
                "Iterations          : 22698",
                "Objective value     : -1.07282263e+00",
                "SoPlex status       : problem is solved [optimal]",
                "  Objective value   : -1.07282263e+00",
            ]
        ),
        encoding="utf-8",
    )

    summary = parse_soplex_output(output)

    assert summary.is_optimal is True
    assert summary.status == "problem is solved [optimal]"
    assert summary.objective_value == pytest.approx(-1.07282263)
    assert summary.solving_time_seconds == pytest.approx(27.04)
    assert summary.iterations == 22698
    assert summary.diagnostic == "optimal"


def test_compare_soplex_objectives_reports_tolerance(tmp_path: Path) -> None:
    reference = tmp_path / "reference.out"
    candidate = tmp_path / "candidate.out"
    for path, value in [(reference, "-1.00000000e+00"), (candidate, "-1.00000001e+00")]:
        path.write_text(
            "\n".join(
                [
                    "SoPlex status       : problem is solved [optimal]",
                    f"Objective value     : {value}",
                ]
            ),
            encoding="utf-8",
        )

    comparison = compare_soplex_objectives(reference, candidate, tolerance=1e-7)

    assert comparison.both_optimal is True
    assert comparison.absolute_difference == pytest.approx(1e-8)
    assert comparison.within_tolerance is True

    output = write_soplex_summary_json(comparison, tmp_path / "comparison.json")
    assert output.exists()


def test_parse_soplex_output_marks_rational_numerical_difficulty(tmp_path: Path) -> None:
    output = tmp_path / "exact.lp.out"
    output.write_text(
        "\n".join(
            [
                "SoPlex status       : error [unspecified]",
                "Objective value     : 0.00000000e+00",
                " --- termination despite violations (numerical difficulties, bound range = 3e+10)",
                "Max. bound violation = 1/4",
                "Max. row violation = 3/8",
                "Solution (rational) : ",
                "  Objective value   : 0",
                "Iterations          : 27198",
                "Numerics            :",
                "  Condition Number  : 0.00000000e+00",
            ]
        ),
        encoding="utf-8",
    )

    summary = parse_soplex_output(output)

    assert summary.is_optimal is False
    assert summary.status == "error [unspecified]"
    assert summary.solution_type == "rational"
    assert summary.objective_value == pytest.approx(0.0)
    assert summary.max_bound_violation == pytest.approx(0.25)
    assert summary.max_row_violation == pytest.approx(0.375)
    assert summary.termination_despite_violations is True
    assert summary.diagnostic == "rational_numerical_difficulty"


def test_parse_soplex_output_marks_optimal_with_unscaled_violations(tmp_path: Path) -> None:
    output = tmp_path / "unscaled.lp.out"
    output.write_text(
        "\n".join(
            [
                "SoPlex status       : problem is solved [optimal with unscaled violations]",
                "Solving time (sec)  : 30.41",
                "Iterations          : 26354",
                "Solution (real)     : ",
                "  Objective value   : 1.37026651e-02",
                "Numerics            :",
                "  Condition Number  : 6.27512119e+16",
            ]
        ),
        encoding="utf-8",
    )

    summary = parse_soplex_output(output)

    assert summary.is_optimal is False
    assert summary.status == "problem is solved [optimal with unscaled violations]"
    assert summary.solution_type == "real"
    assert summary.objective_value == pytest.approx(0.0137026651)
    assert summary.condition_number == pytest.approx(6.27512119e16)
    assert summary.diagnostic == "optimal_with_unscaled_violations"


def test_parse_existing_matlab_opn_soplex_output_when_available() -> None:
    output = (
        ProjectPaths.discover().opn_run_dir
        / "Simulation_dilutionOPN_ALPHA_FULL_PROJECT_mu0p1_media4_ratio1em08_misfolddefault_noMisfoldEq_noRiboEq_PP.lp.float.out"
    )
    if not output.exists():
        pytest.skip("Local MATLAB OPN SoPlex output has not been generated on this machine.")

    summary = parse_soplex_output(output)

    assert summary.is_optimal is True
    assert summary.objective_value == pytest.approx(-1.07282263)
    assert summary.iterations == 22698


def test_parse_existing_python_exact_soplex_output_when_available() -> None:
    output = (
        ProjectPaths.discover().local_runs_dir
        / "pichia_python"
        / "python_OPN_ALPHA_FULL_PROJECT_mu0p10_media4_ratio1em08_noMisfoldEq_noRiboEq.lp.exact.out"
    )
    if not output.exists():
        pytest.skip("Local Python exact SoPlex output has not been generated on this machine.")

    summary = parse_soplex_output(output)

    assert summary.is_optimal is False
    assert summary.status == "error [unspecified]"
    assert summary.solution_type == "rational"
    assert summary.diagnostic == "rational_numerical_difficulty"
    assert summary.termination_despite_violations is True


def test_compare_existing_dsb_secretory_matlab_and_python_soplex_outputs_when_available() -> None:
    paths = ProjectPaths.discover()
    matlab_output = (
        paths.local_runs_dir
        / "pichia_python"
        / "dsb_secretory_baseline"
        / "Simulation_dilutionDSB80_TEST_LEADER_mu0p10_media4_dsb_secretory_noMisfoldEq_noRiboEq_PP.lp.float.out"
    )
    python_output = (
        paths.local_runs_dir
        / "pichia_python"
        / "dsb_secretory_smoke"
        / "python_DSB80_TEST_LEADER_mu0p10_dsb_secretory_target_reference_constraints.lp.float.out"
    )
    if not matlab_output.exists() or not python_output.exists():
        pytest.skip("Local DSB-only MATLAB/Python SoPlex outputs have not been generated on this machine.")

    matlab_summary = parse_soplex_output(matlab_output)
    python_summary = parse_soplex_output(python_output)
    comparison = compare_soplex_objectives(matlab_output, python_output, tolerance=1e-12)

    assert matlab_summary.status == "problem is solved [optimal with unscaled violations]"
    assert python_summary.status == matlab_summary.status
    assert matlab_summary.diagnostic == "optimal_with_unscaled_violations"
    assert python_summary.diagnostic == matlab_summary.diagnostic
    assert comparison.within_tolerance is True
    assert comparison.absolute_difference == pytest.approx(0.0)
    assert comparison.both_optimal is False


def test_parse_existing_dsb_og_secretory_python_soplex_output_when_available() -> None:
    output = (
        ProjectPaths.discover().local_runs_dir
        / "pichia_python"
        / "dsb_og_secretory_smoke"
        / "python_DSB_OG80_TEST_LEADER_mu0p10_dsb_og_secretory_target_reference_constraints.lp.float.out"
    )
    if not output.exists():
        pytest.skip("Local DSB+OG Python SoPlex output has not been generated on this machine.")

    summary = parse_soplex_output(output)

    assert summary.status == "problem is solved [optimal]"
    assert summary.is_optimal is True
    assert summary.solution_type == "real"
    assert summary.diagnostic == "optimal"
    assert summary.objective_value == pytest.approx(0.0132246623)


def test_compare_existing_dsb_og_secretory_matlab_and_python_soplex_outputs_when_available() -> None:
    paths = ProjectPaths.discover()
    matlab_output = (
        paths.local_runs_dir
        / "pichia_python"
        / "dsb_og_secretory_baseline"
        / "Simulation_dilutionDSB_OG80_TEST_LEADER_mu0p10_media4_dsb_og_secretory_noMisfoldEq_noRiboEq_PP.lp.float.out"
    )
    python_output = (
        paths.local_runs_dir
        / "pichia_python"
        / "dsb_og_secretory_smoke"
        / "python_DSB_OG80_TEST_LEADER_mu0p10_dsb_og_secretory_target_reference_constraints.lp.float.out"
    )
    if not matlab_output.exists() or not python_output.exists():
        pytest.skip("Local DSB+OG MATLAB/Python SoPlex outputs have not been generated on this machine.")

    matlab_summary = parse_soplex_output(matlab_output)
    python_summary = parse_soplex_output(python_output)
    comparison = compare_soplex_objectives(matlab_output, python_output, tolerance=1e-12)

    assert matlab_summary.status == "problem is solved [optimal]"
    assert python_summary.status == matlab_summary.status
    assert matlab_summary.diagnostic == "optimal"
    assert python_summary.diagnostic == matlab_summary.diagnostic
    assert comparison.within_tolerance is True
    assert comparison.absolute_difference == pytest.approx(0.0)
    assert comparison.both_optimal is True


def test_parse_existing_dsb_ng_secretory_python_soplex_output_when_available() -> None:
    output = (
        ProjectPaths.discover().local_runs_dir
        / "pichia_python"
        / "dsb_ng_secretory_smoke"
        / "python_DSB_NG80_TEST_LEADER_mu0p10_dsb_ng_secretory_target_reference_constraints.lp.float.out"
    )
    if not output.exists():
        pytest.skip("Local DSB+NG Python SoPlex output has not been generated on this machine.")

    summary = parse_soplex_output(output)

    assert summary.status == "problem is solved [optimal]"
    assert summary.is_optimal is True
    assert summary.solution_type == "real"
    assert summary.diagnostic == "optimal"
    assert summary.objective_value == pytest.approx(0.0121802802)


def test_compare_existing_dsb_ng_secretory_matlab_and_python_soplex_outputs_when_available() -> None:
    paths = ProjectPaths.discover()
    matlab_output = (
        paths.local_runs_dir
        / "pichia_python"
        / "dsb_ng_secretory_baseline"
        / "Simulation_dilutionDSB_NG80_TEST_LEADER_mu0p10_media4_dsb_ng_secretory_noMisfoldEq_noRiboEq_PP.lp.float.out"
    )
    python_output = (
        paths.local_runs_dir
        / "pichia_python"
        / "dsb_ng_secretory_smoke"
        / "python_DSB_NG80_TEST_LEADER_mu0p10_dsb_ng_secretory_target_reference_constraints.lp.float.out"
    )
    if not matlab_output.exists() or not python_output.exists():
        pytest.skip("Local DSB+NG MATLAB/Python SoPlex outputs have not been generated on this machine.")

    matlab_summary = parse_soplex_output(matlab_output)
    python_summary = parse_soplex_output(python_output)
    comparison = compare_soplex_objectives(matlab_output, python_output, tolerance=1e-12)

    assert matlab_summary.status == "problem is solved [optimal]"
    assert python_summary.status == matlab_summary.status
    assert matlab_summary.diagnostic == "optimal"
    assert python_summary.diagnostic == matlab_summary.diagnostic
    assert comparison.both_optimal is True
    assert comparison.within_tolerance is True
    assert comparison.absolute_difference == pytest.approx(0.0)


def test_parse_existing_dsb_ng_og_secretory_python_soplex_output_when_available() -> None:
    output = (
        ProjectPaths.discover().local_runs_dir
        / "pichia_python"
        / "dsb_ng_og_secretory_smoke"
        / "python_DSB_NG_OG80_TEST_LEADER_mu0p10_dsb_ng_og_secretory_target_reference_constraints.lp.float.out"
    )
    if not output.exists():
        pytest.skip("Local DSB+NG+OG Python SoPlex output has not been generated on this machine.")

    summary = parse_soplex_output(output)

    assert summary.status == "problem is solved [optimal]"
    assert summary.is_optimal is True
    assert summary.solution_type == "real"
    assert summary.diagnostic == "optimal"
    assert summary.objective_value is not None
    assert summary.objective_value > 0


def test_compare_existing_dsb_ng_og_secretory_matlab_and_python_soplex_outputs_when_available() -> None:
    paths = ProjectPaths.discover()
    matlab_output = (
        paths.local_runs_dir
        / "pichia_python"
        / "dsb_ng_og_secretory_baseline"
        / "Simulation_dilutionDSB_NG_OG80_TEST_LEADER_mu0p10_media4_dsb_ng_og_secretory_noMisfoldEq_noRiboEq_PP.lp.float.out"
    )
    python_output = (
        paths.local_runs_dir
        / "pichia_python"
        / "dsb_ng_og_secretory_smoke"
        / "python_DSB_NG_OG80_TEST_LEADER_mu0p10_dsb_ng_og_secretory_target_reference_constraints.lp.float.out"
    )
    if not matlab_output.exists() or not python_output.exists():
        pytest.skip("Local DSB+NG+OG MATLAB/Python SoPlex outputs have not been generated on this machine.")

    matlab_summary = parse_soplex_output(matlab_output)
    python_summary = parse_soplex_output(python_output)
    comparison = compare_soplex_objectives(matlab_output, python_output, tolerance=1e-12)

    assert matlab_summary.status == "problem is solved [optimal]"
    assert python_summary.status == matlab_summary.status
    assert matlab_summary.diagnostic == "optimal"
    assert python_summary.diagnostic == matlab_summary.diagnostic
    assert comparison.both_optimal is True
    assert comparison.within_tolerance is True
    assert comparison.absolute_difference == pytest.approx(0.0)


def test_parse_existing_ng_og_secretory_python_soplex_output_when_available() -> None:
    output = (
        ProjectPaths.discover().local_runs_dir
        / "pichia_python"
        / "ng_og_secretory_smoke"
        / "python_NG_OG80_TEST_LEADER_mu0p10_ng_og_secretory_target_reference_constraints.lp.float.out"
    )
    if not output.exists():
        pytest.skip("Local NG+OG Python SoPlex output has not been generated on this machine.")

    summary = parse_soplex_output(output)

    assert summary.status == "problem is solved [optimal]"
    assert summary.is_optimal is True
    assert summary.solution_type == "real"
    assert summary.diagnostic == "optimal"
    assert summary.objective_value == pytest.approx(0.0121292303)


def test_compare_existing_ng_og_secretory_matlab_and_python_soplex_outputs_when_available() -> None:
    paths = ProjectPaths.discover()
    matlab_output = (
        paths.local_runs_dir
        / "pichia_python"
        / "ng_og_secretory_baseline"
        / "Simulation_dilutionNG_OG80_TEST_LEADER_mu0p10_media4_ng_og_secretory_noMisfoldEq_noRiboEq_PP.lp.float.out"
    )
    python_output = (
        paths.local_runs_dir
        / "pichia_python"
        / "ng_og_secretory_smoke"
        / "python_NG_OG80_TEST_LEADER_mu0p10_ng_og_secretory_target_reference_constraints.lp.float.out"
    )
    if not matlab_output.exists() or not python_output.exists():
        pytest.skip("Local NG+OG MATLAB/Python SoPlex outputs have not been generated on this machine.")

    matlab_summary = parse_soplex_output(matlab_output)
    python_summary = parse_soplex_output(python_output)
    comparison = compare_soplex_objectives(matlab_output, python_output, tolerance=1e-12)

    assert matlab_summary.status == "problem is solved [optimal]"
    assert python_summary.status == matlab_summary.status
    assert matlab_summary.diagnostic == "optimal"
    assert python_summary.diagnostic == matlab_summary.diagnostic
    assert comparison.both_optimal is True
    assert comparison.within_tolerance is True
    assert comparison.absolute_difference == pytest.approx(0.0)


def test_parse_existing_ng_secretory_python_soplex_output_when_available() -> None:
    output = (
        ProjectPaths.discover().local_runs_dir
        / "pichia_python"
        / "ng_secretory_smoke"
        / "python_NG80_TEST_LEADER_mu0p10_ng_secretory_target_reference_constraints.lp.float.out"
    )
    if not output.exists():
        pytest.skip("Local NG-only Python SoPlex output has not been generated on this machine.")

    summary = parse_soplex_output(output)

    assert summary.status == "problem is solved [optimal]"
    assert summary.is_optimal is True
    assert summary.solution_type == "real"
    assert summary.diagnostic == "optimal"
    assert summary.objective_value == pytest.approx(0.00706758633)


def test_compare_existing_ng_secretory_matlab_and_python_soplex_outputs_when_available() -> None:
    paths = ProjectPaths.discover()
    matlab_output = (
        paths.local_runs_dir
        / "pichia_python"
        / "ng_secretory_baseline"
        / "Simulation_dilutionNG80_TEST_LEADER_mu0p10_media4_ng_secretory_noMisfoldEq_noRiboEq_PP.lp.float.out"
    )
    python_output = (
        paths.local_runs_dir
        / "pichia_python"
        / "ng_secretory_smoke"
        / "python_NG80_TEST_LEADER_mu0p10_ng_secretory_target_reference_constraints.lp.float.out"
    )
    if not matlab_output.exists() or not python_output.exists():
        pytest.skip("Local NG-only MATLAB/Python SoPlex outputs have not been generated on this machine.")

    matlab_summary = parse_soplex_output(matlab_output)
    python_summary = parse_soplex_output(python_output)
    comparison = compare_soplex_objectives(matlab_output, python_output, tolerance=1e-12)

    assert matlab_summary.status == "problem is solved [optimal]"
    assert python_summary.status == matlab_summary.status
    assert matlab_summary.diagnostic == "optimal"
    assert python_summary.diagnostic == matlab_summary.diagnostic
    assert comparison.both_optimal is True
    assert comparison.within_tolerance is True
    assert comparison.absolute_difference == pytest.approx(0.0)


def test_parse_existing_og_secretory_python_soplex_output_when_available() -> None:
    output = (
        ProjectPaths.discover().local_runs_dir
        / "pichia_python"
        / "og_secretory_smoke"
        / "python_OG80_TEST_LEADER_mu0p10_og_secretory_target_reference_constraints.lp.float.out"
    )
    if not output.exists():
        pytest.skip("Local OG-only Python SoPlex output has not been generated on this machine.")

    summary = parse_soplex_output(output)

    assert summary.status == "problem is solved [optimal]"
    assert summary.is_optimal is True
    assert summary.solution_type == "real"
    assert summary.diagnostic == "optimal"
    assert summary.objective_value == pytest.approx(0.0131338977)


def test_compare_existing_og_secretory_matlab_and_python_soplex_outputs_when_available() -> None:
    paths = ProjectPaths.discover()
    matlab_output = (
        paths.local_runs_dir
        / "pichia_python"
        / "og_secretory_baseline"
        / "Simulation_dilutionOG80_TEST_LEADER_mu0p10_media4_og_secretory_noMisfoldEq_noRiboEq_PP.lp.float.out"
    )
    python_output = (
        paths.local_runs_dir
        / "pichia_python"
        / "og_secretory_smoke"
        / "python_OG80_TEST_LEADER_mu0p10_og_secretory_target_reference_constraints.lp.float.out"
    )
    if not matlab_output.exists() or not python_output.exists():
        pytest.skip("Local OG-only MATLAB/Python SoPlex outputs have not been generated on this machine.")

    matlab_summary = parse_soplex_output(matlab_output)
    python_summary = parse_soplex_output(python_output)
    comparison = compare_soplex_objectives(matlab_output, python_output, tolerance=1e-12)

    assert matlab_summary.status == "problem is solved [optimal]"
    assert python_summary.status == matlab_summary.status
    assert matlab_summary.diagnostic == "optimal"
    assert python_summary.diagnostic == matlab_summary.diagnostic
    assert comparison.both_optimal is True
    assert comparison.within_tolerance is True
    assert comparison.absolute_difference == pytest.approx(0.0)


def test_parse_existing_gpi_secretory_python_soplex_output_when_available() -> None:
    output = (
        ProjectPaths.discover().local_runs_dir
        / "pichia_python"
        / "gpi_secretory_smoke"
        / "python_GPI80_TEST_LEADER_mu0p10_gpi_secretory_target_reference_constraints.lp.float.out"
    )
    if not output.exists():
        pytest.skip("Local GPI-only Python SoPlex output has not been generated on this machine.")

    summary = parse_soplex_output(output)

    assert summary.status == "problem is solved [optimal]"
    assert summary.is_optimal is True
    assert summary.solution_type == "real"
    assert summary.diagnostic == "optimal"
    assert summary.objective_value == pytest.approx(0.00177893632)


def test_compare_existing_gpi_secretory_matlab_and_python_soplex_outputs_when_available() -> None:
    paths = ProjectPaths.discover()
    matlab_output = (
        paths.local_runs_dir
        / "pichia_python"
        / "gpi_secretory_baseline"
        / "Simulation_dilutionGPI80_TEST_LEADER_mu0p10_media4_gpi_secretory_noMisfoldEq_noRiboEq_PP.lp.float.out"
    )
    python_output = (
        paths.local_runs_dir
        / "pichia_python"
        / "gpi_secretory_smoke"
        / "python_GPI80_TEST_LEADER_mu0p10_gpi_secretory_target_reference_constraints.lp.float.out"
    )
    if not matlab_output.exists() or not python_output.exists():
        pytest.skip("Local GPI-only MATLAB/Python SoPlex outputs have not been generated on this machine.")

    matlab_summary = parse_soplex_output(matlab_output)
    python_summary = parse_soplex_output(python_output)
    comparison = compare_soplex_objectives(matlab_output, python_output, tolerance=1e-12)

    assert matlab_summary.status == "problem is solved [optimal]"
    assert python_summary.status == matlab_summary.status
    assert matlab_summary.diagnostic == "optimal"
    assert python_summary.diagnostic == matlab_summary.diagnostic
    assert comparison.both_optimal is True
    assert comparison.within_tolerance is True
    assert comparison.absolute_difference == pytest.approx(0.0)


def test_parse_existing_gpi_og_secretory_python_soplex_output_when_available() -> None:
    output = (
        ProjectPaths.discover().local_runs_dir
        / "pichia_python"
        / "gpi_og_secretory_smoke"
        / "python_GPI_OG80_TEST_LEADER_mu0p10_gpi_og_secretory_target_reference_constraints.lp.float.out"
    )
    if not output.exists():
        pytest.skip("Local GPI+OG Python SoPlex output has not been generated on this machine.")

    summary = parse_soplex_output(output)

    assert summary.status == "problem is solved [optimal]"
    assert summary.is_optimal is True
    assert summary.solution_type == "real"
    assert summary.diagnostic == "optimal"
    assert summary.objective_value == pytest.approx(0.00296832742)


def test_compare_existing_gpi_og_secretory_matlab_and_python_soplex_outputs_when_available() -> None:
    paths = ProjectPaths.discover()
    matlab_output = (
        paths.local_runs_dir
        / "pichia_python"
        / "gpi_og_secretory_baseline"
        / "Simulation_dilutionGPI_OG80_TEST_LEADER_mu0p10_media4_gpi_og_secretory_noMisfoldEq_noRiboEq_PP.lp.float.out"
    )
    python_output = (
        paths.local_runs_dir
        / "pichia_python"
        / "gpi_og_secretory_smoke"
        / "python_GPI_OG80_TEST_LEADER_mu0p10_gpi_og_secretory_target_reference_constraints.lp.float.out"
    )
    if not matlab_output.exists() or not python_output.exists():
        pytest.skip("Local GPI+OG MATLAB/Python SoPlex outputs have not been generated on this machine.")

    matlab_summary = parse_soplex_output(matlab_output)
    python_summary = parse_soplex_output(python_output)
    comparison = compare_soplex_objectives(matlab_output, python_output, tolerance=1e-12)

    assert matlab_summary.status == "problem is solved [optimal]"
    assert python_summary.status == matlab_summary.status
    assert matlab_summary.diagnostic == "optimal"
    assert python_summary.diagnostic == matlab_summary.diagnostic
    assert comparison.both_optimal is True
    assert comparison.within_tolerance is True
    assert comparison.absolute_difference == pytest.approx(0.0)


def test_parse_existing_gpi_ng_secretory_python_soplex_output_when_available() -> None:
    output = (
        ProjectPaths.discover().local_runs_dir
        / "pichia_python"
        / "gpi_ng_secretory_smoke"
        / "python_GPI_NG80_TEST_LEADER_mu0p10_gpi_ng_secretory_target_reference_constraints.lp.float.out"
    )
    if not output.exists():
        pytest.skip("Local GPI+NG Python SoPlex output has not been generated on this machine.")

    summary = parse_soplex_output(output)

    assert summary.status == "problem is solved [optimal]"
    assert summary.is_optimal is True
    assert summary.solution_type == "real"
    assert summary.diagnostic == "optimal"
    assert summary.objective_value == pytest.approx(0.00691270370)


def test_compare_existing_gpi_ng_secretory_matlab_and_python_soplex_outputs_when_available() -> None:
    paths = ProjectPaths.discover()
    matlab_output = (
        paths.local_runs_dir
        / "pichia_python"
        / "gpi_ng_secretory_baseline"
        / "Simulation_dilutionGPI_NG80_TEST_LEADER_mu0p10_media4_gpi_ng_secretory_noMisfoldEq_noRiboEq_PP.lp.float.out"
    )
    python_output = (
        paths.local_runs_dir
        / "pichia_python"
        / "gpi_ng_secretory_smoke"
        / "python_GPI_NG80_TEST_LEADER_mu0p10_gpi_ng_secretory_target_reference_constraints.lp.float.out"
    )
    if not matlab_output.exists() or not python_output.exists():
        pytest.skip("Local GPI+NG MATLAB/Python SoPlex outputs have not been generated on this machine.")

    matlab_summary = parse_soplex_output(matlab_output)
    python_summary = parse_soplex_output(python_output)
    comparison = compare_soplex_objectives(matlab_output, python_output, tolerance=1e-12)

    assert matlab_summary.status == "problem is solved [optimal]"
    assert python_summary.status == matlab_summary.status
    assert matlab_summary.diagnostic == "optimal"
    assert python_summary.diagnostic == matlab_summary.diagnostic
    assert comparison.both_optimal is True
    assert comparison.within_tolerance is True
    assert comparison.absolute_difference == pytest.approx(0.0)


def test_parse_existing_gpi_ng_og_secretory_python_soplex_output_when_available() -> None:
    output = (
        ProjectPaths.discover().local_runs_dir
        / "pichia_python"
        / "gpi_ng_og_secretory_smoke"
        / "python_GPI_NG_OG80_TEST_LEADER_mu0p10_gpi_ng_og_secretory_target_reference_constraints.lp.float.out"
    )
    if not output.exists():
        pytest.skip("Local GPI+NG+OG Python SoPlex output has not been generated on this machine.")

    summary = parse_soplex_output(output)

    assert summary.status == "problem is solved [optimal]"
    assert summary.is_optimal is True
    assert summary.solution_type == "real"
    assert summary.diagnostic == "optimal"
    assert summary.objective_value == pytest.approx(0.00680967942)


def test_compare_existing_gpi_ng_og_secretory_matlab_and_python_soplex_outputs_when_available() -> None:
    paths = ProjectPaths.discover()
    matlab_output = (
        paths.local_runs_dir
        / "pichia_python"
        / "gpi_ng_og_secretory_baseline"
        / "Simulation_dilutionGPI_NG_OG80_TEST_LEADER_mu0p10_media4_gpi_ng_og_secretory_noMisfoldEq_noRiboEq_PP.lp.float.out"
    )
    python_output = (
        paths.local_runs_dir
        / "pichia_python"
        / "gpi_ng_og_secretory_smoke"
        / "python_GPI_NG_OG80_TEST_LEADER_mu0p10_gpi_ng_og_secretory_target_reference_constraints.lp.float.out"
    )
    if not matlab_output.exists() or not python_output.exists():
        pytest.skip("Local GPI+NG+OG MATLAB/Python SoPlex outputs have not been generated on this machine.")

    matlab_summary = parse_soplex_output(matlab_output)
    python_summary = parse_soplex_output(python_output)
    comparison = compare_soplex_objectives(matlab_output, python_output, tolerance=1e-12)

    assert matlab_summary.status == "problem is solved [optimal]"
    assert python_summary.status == matlab_summary.status
    assert matlab_summary.diagnostic == "optimal"
    assert python_summary.diagnostic == matlab_summary.diagnostic
    assert comparison.both_optimal is True
    assert comparison.within_tolerance is True
    assert comparison.absolute_difference == pytest.approx(0.0)


def test_parse_existing_transmembrane_erm_python_soplex_output_when_available() -> None:
    output = (
        ProjectPaths.discover().local_runs_dir
        / "pichia_python"
        / "transmembrane_erm_smoke"
        / "python_TM80_TEST_LEADER_mu0p10_transmembrane_erm_target_reference_constraints.lp.float.out"
    )
    if not output.exists():
        pytest.skip("Local transmembrane ERM Python SoPlex output has not been generated on this machine.")

    summary = parse_soplex_output(output)

    assert summary.status == "problem is solved [optimal]"
    assert summary.is_optimal is True
    assert summary.solution_type == "real"
    assert summary.diagnostic == "optimal"
    assert summary.objective_value == pytest.approx(0.0132063581)


def test_compare_existing_transmembrane_erm_matlab_and_python_soplex_outputs_when_available() -> None:
    paths = ProjectPaths.discover()
    matlab_output = (
        paths.local_runs_dir
        / "pichia_python"
        / "transmembrane_erm_baseline"
        / "Simulation_dilutionTM80_TEST_LEADER_mu0p10_media4_transmembrane_erm_noMisfoldEq_noRiboEq_PP.lp.float.out"
    )
    python_output = (
        paths.local_runs_dir
        / "pichia_python"
        / "transmembrane_erm_smoke"
        / "python_TM80_TEST_LEADER_mu0p10_transmembrane_erm_target_reference_constraints.lp.float.out"
    )
    if not matlab_output.exists() or not python_output.exists():
        pytest.skip("Local transmembrane ERM MATLAB/Python SoPlex outputs have not been generated on this machine.")

    matlab_summary = parse_soplex_output(matlab_output)
    python_summary = parse_soplex_output(python_output)
    comparison = compare_soplex_objectives(matlab_output, python_output, tolerance=1e-12)

    assert matlab_summary.status == "problem is solved [optimal]"
    assert python_summary.status == matlab_summary.status
    assert matlab_summary.diagnostic == "optimal"
    assert python_summary.diagnostic == matlab_summary.diagnostic
    assert comparison.both_optimal is True
    assert comparison.within_tolerance is True
    assert comparison.absolute_difference == pytest.approx(0.0)


def test_parse_existing_transmembrane_er_python_soplex_output_when_available() -> None:
    output = (
        ProjectPaths.discover().local_runs_dir
        / "pichia_python"
        / "transmembrane_er_smoke"
        / "python_TM80ER_TEST_LEADER_mu0p10_transmembrane_er_target_reference_constraints.lp.float.out"
    )
    if not output.exists():
        pytest.skip("Local transmembrane ER Python SoPlex output has not been generated on this machine.")

    summary = parse_soplex_output(output)

    assert summary.status == "problem is solved [optimal]"
    assert summary.is_optimal is True
    assert summary.solution_type == "real"
    assert summary.diagnostic == "optimal"
    assert summary.objective_value == pytest.approx(0.0132063581)


def test_compare_existing_transmembrane_er_matlab_and_python_soplex_outputs_when_available() -> None:
    paths = ProjectPaths.discover()
    matlab_output = (
        paths.local_runs_dir
        / "pichia_python"
        / "transmembrane_er_baseline"
        / "Simulation_dilutionTM80ER_TEST_LEADER_mu0p10_media4_transmembrane_er_noMisfoldEq_noRiboEq_PP.lp.float.out"
    )
    python_output = (
        paths.local_runs_dir
        / "pichia_python"
        / "transmembrane_er_smoke"
        / "python_TM80ER_TEST_LEADER_mu0p10_transmembrane_er_target_reference_constraints.lp.float.out"
    )
    if not matlab_output.exists() or not python_output.exists():
        pytest.skip("Local transmembrane ER MATLAB/Python SoPlex outputs have not been generated on this machine.")

    matlab_summary = parse_soplex_output(matlab_output)
    python_summary = parse_soplex_output(python_output)
    comparison = compare_soplex_objectives(matlab_output, python_output, tolerance=1e-12)

    assert matlab_summary.status == "problem is solved [optimal]"
    assert python_summary.status == matlab_summary.status
    assert matlab_summary.diagnostic == "optimal"
    assert python_summary.diagnostic == matlab_summary.diagnostic
    assert comparison.both_optimal is True
    assert comparison.within_tolerance is True
    assert comparison.absolute_difference == pytest.approx(0.0)


def test_parse_existing_transmembrane_c_python_soplex_output_when_available() -> None:
    output = (
        ProjectPaths.discover().local_runs_dir
        / "pichia_python"
        / "transmembrane_c_smoke"
        / "python_TM80C_TEST_LEADER_mu0p10_transmembrane_c_target_reference_constraints.lp.float.out"
    )
    if not output.exists():
        pytest.skip("Local transmembrane cytosolic Python SoPlex output has not been generated on this machine.")

    summary = parse_soplex_output(output)

    assert summary.status == "problem is solved [optimal]"
    assert summary.is_optimal is True
    assert summary.solution_type == "real"
    assert summary.diagnostic == "optimal"
    assert summary.objective_value == pytest.approx(0.0132452337)


def test_compare_existing_transmembrane_c_matlab_and_python_soplex_outputs_when_available() -> None:
    paths = ProjectPaths.discover()
    matlab_output = (
        paths.local_runs_dir
        / "pichia_python"
        / "transmembrane_c_baseline"
        / "Simulation_dilutionTM80C_TEST_LEADER_mu0p10_media4_transmembrane_c_noMisfoldEq_noRiboEq_PP.lp.float.out"
    )
    python_output = (
        paths.local_runs_dir
        / "pichia_python"
        / "transmembrane_c_smoke"
        / "python_TM80C_TEST_LEADER_mu0p10_transmembrane_c_target_reference_constraints.lp.float.out"
    )
    if not matlab_output.exists() or not python_output.exists():
        pytest.skip("Local transmembrane cytosolic MATLAB/Python SoPlex outputs have not been generated on this machine.")

    matlab_summary = parse_soplex_output(matlab_output)
    python_summary = parse_soplex_output(python_output)
    comparison = compare_soplex_objectives(matlab_output, python_output, tolerance=1e-12)

    assert matlab_summary.status == "problem is solved [optimal]"
    assert python_summary.status == matlab_summary.status
    assert matlab_summary.diagnostic == "optimal"
    assert python_summary.diagnostic == matlab_summary.diagnostic
    assert comparison.both_optimal is True
    assert comparison.within_tolerance is True
    assert comparison.absolute_difference == pytest.approx(0.0)


def test_parse_existing_transmembrane_e_python_soplex_output_when_available() -> None:
    output = (
        ProjectPaths.discover().local_runs_dir
        / "pichia_python"
        / "transmembrane_e_smoke"
        / "python_TME80_TEST_LEADER_mu0p10_transmembrane_e_target_reference_constraints.lp.float.out"
    )
    if not output.exists():
        pytest.skip("Local transmembrane extracellular Python SoPlex output has not been generated on this machine.")

    summary = parse_soplex_output(output)

    assert summary.status == "problem is solved [optimal with unscaled violations]"
    assert summary.is_optimal is False
    assert summary.solution_type == "real"
    assert summary.diagnostic == "optimal_with_unscaled_violations"
    assert summary.objective_value == pytest.approx(0.0131904717)


def test_compare_existing_transmembrane_e_matlab_and_python_soplex_outputs_when_available() -> None:
    paths = ProjectPaths.discover()
    matlab_output = (
        paths.local_runs_dir
        / "pichia_python"
        / "transmembrane_e_baseline"
        / "Simulation_dilutionTME80_TEST_LEADER_mu0p10_media4_transmembrane_e_noMisfoldEq_noRiboEq_PP.lp.float.out"
    )
    python_output = (
        paths.local_runs_dir
        / "pichia_python"
        / "transmembrane_e_smoke"
        / "python_TME80_TEST_LEADER_mu0p10_transmembrane_e_target_reference_constraints.lp.float.out"
    )
    if not matlab_output.exists() or not python_output.exists():
        pytest.skip("Local transmembrane extracellular MATLAB/Python SoPlex outputs have not been generated on this machine.")

    matlab_summary = parse_soplex_output(matlab_output)
    python_summary = parse_soplex_output(python_output)
    comparison = compare_soplex_objectives(matlab_output, python_output, tolerance=1e-12)

    assert matlab_summary.status == "problem is solved [optimal with unscaled violations]"
    assert python_summary.status == matlab_summary.status
    assert matlab_summary.diagnostic == "optimal_with_unscaled_violations"
    assert python_summary.diagnostic == matlab_summary.diagnostic
    assert comparison.both_optimal is False
    assert comparison.within_tolerance is True
    assert comparison.absolute_difference == pytest.approx(0.0)


def test_parse_existing_co_translation_secretory_python_soplex_output_when_available() -> None:
    output = (
        ProjectPaths.discover().local_runs_dir
        / "pichia_python"
        / "co_translation_smoke"
        / "python_CO80_TEST_LEADER_mu0p10_co_translation_secretory_target_reference_constraints.lp.float.out"
    )
    if not output.exists():
        pytest.skip("Local co-translation Python SoPlex output has not been generated on this machine.")

    summary = parse_soplex_output(output)

    assert summary.status == "problem is solved [optimal]"
    assert summary.is_optimal is True
    assert summary.solution_type == "real"
    assert summary.diagnostic == "optimal"
    assert summary.objective_value == pytest.approx(0.013179831)


def test_compare_existing_co_translation_secretory_matlab_and_python_soplex_outputs_when_available() -> None:
    paths = ProjectPaths.discover()
    matlab_output = (
        paths.local_runs_dir
        / "pichia_python"
        / "co_translation_baseline"
        / "Simulation_dilutionCO80_TEST_LEADER_mu0p10_media4_co_translation_noMisfoldEq_noRiboEq_PP.lp.float.out"
    )
    python_output = (
        paths.local_runs_dir
        / "pichia_python"
        / "co_translation_smoke"
        / "python_CO80_TEST_LEADER_mu0p10_co_translation_secretory_target_reference_constraints.lp.float.out"
    )
    if not matlab_output.exists() or not python_output.exists():
        pytest.skip("Local co-translation MATLAB/Python SoPlex outputs have not been generated on this machine.")

    matlab_summary = parse_soplex_output(matlab_output)
    python_summary = parse_soplex_output(python_output)
    comparison = compare_soplex_objectives(matlab_output, python_output, tolerance=1e-12)

    assert matlab_summary.status == "problem is solved [optimal]"
    assert python_summary.status == matlab_summary.status
    assert matlab_summary.diagnostic == "optimal"
    assert python_summary.diagnostic == matlab_summary.diagnostic
    assert comparison.both_optimal is True
    assert comparison.within_tolerance is True
    assert comparison.absolute_difference == pytest.approx(0.0)


def test_parse_existing_soluble_ce_secretory_python_soplex_output_when_available() -> None:
    output = (
        ProjectPaths.discover().local_runs_dir
        / "pichia_python"
        / "ce_secretory_smoke"
        / "python_CE80_TEST_LEADER_mu0p10_soluble_ce_secretory_target_reference_constraints.lp.float.out"
    )
    if not output.exists():
        pytest.skip("Local soluble CE Python SoPlex output has not been generated on this machine.")

    summary = parse_soplex_output(output)

    assert summary.status == "problem is solved [optimal]"
    assert summary.is_optimal is True
    assert summary.solution_type == "real"
    assert summary.diagnostic == "optimal"
    assert summary.objective_value == pytest.approx(0.0131918537)


def test_compare_existing_soluble_ce_secretory_matlab_and_python_soplex_outputs_when_available() -> None:
    paths = ProjectPaths.discover()
    matlab_output = (
        paths.local_runs_dir
        / "pichia_python"
        / "ce_secretory_baseline"
        / "Simulation_dilutionCE80_TEST_LEADER_mu0p10_media4_ce_secretory_noMisfoldEq_noRiboEq_PP.lp.float.out"
    )
    python_output = (
        paths.local_runs_dir
        / "pichia_python"
        / "ce_secretory_smoke"
        / "python_CE80_TEST_LEADER_mu0p10_soluble_ce_secretory_target_reference_constraints.lp.float.out"
    )
    if not matlab_output.exists() or not python_output.exists():
        pytest.skip("Local soluble CE MATLAB/Python SoPlex outputs have not been generated on this machine.")

    matlab_summary = parse_soplex_output(matlab_output)
    python_summary = parse_soplex_output(python_output)
    comparison = compare_soplex_objectives(matlab_output, python_output, tolerance=1e-12)

    assert matlab_summary.status == "problem is solved [optimal]"
    assert python_summary.status == matlab_summary.status
    assert matlab_summary.diagnostic == "optimal"
    assert python_summary.diagnostic == matlab_summary.diagnostic
    assert comparison.both_optimal is True
    assert comparison.within_tolerance is True
    assert comparison.absolute_difference == pytest.approx(0.0)


def test_parse_existing_soluble_vm_secretory_python_soplex_output_when_available() -> None:
    output = (
        ProjectPaths.discover().local_runs_dir
        / "pichia_python"
        / "vm_secretory_smoke"
        / "python_VM80_TEST_LEADER_mu0p10_soluble_vm_secretory_target_reference_constraints.lp.float.out"
    )
    if not output.exists():
        pytest.skip("Local soluble VM Python SoPlex output has not been generated on this machine.")

    summary = parse_soplex_output(output)

    assert summary.status == "problem is solved [optimal]"
    assert summary.is_optimal is True
    assert summary.solution_type == "real"
    assert summary.diagnostic == "optimal"
    assert summary.objective_value == pytest.approx(0.0132110308)


def test_compare_existing_soluble_vm_secretory_matlab_and_python_soplex_outputs_when_available() -> None:
    paths = ProjectPaths.discover()
    matlab_output = (
        paths.local_runs_dir
        / "pichia_python"
        / "vm_secretory_baseline"
        / "Simulation_dilutionVM80_TEST_LEADER_mu0p10_media4_vm_secretory_noMisfoldEq_noRiboEq_PP.lp.float.out"
    )
    python_output = (
        paths.local_runs_dir
        / "pichia_python"
        / "vm_secretory_smoke"
        / "python_VM80_TEST_LEADER_mu0p10_soluble_vm_secretory_target_reference_constraints.lp.float.out"
    )
    if not matlab_output.exists() or not python_output.exists():
        pytest.skip("Local soluble VM MATLAB/Python SoPlex outputs have not been generated on this machine.")

    matlab_summary = parse_soplex_output(matlab_output)
    python_summary = parse_soplex_output(python_output)
    comparison = compare_soplex_objectives(matlab_output, python_output, tolerance=1e-12)

    assert matlab_summary.status == "problem is solved [optimal]"
    assert python_summary.status == matlab_summary.status
    assert matlab_summary.diagnostic == "optimal"
    assert python_summary.diagnostic == matlab_summary.diagnostic
    assert comparison.both_optimal is True
    assert comparison.within_tolerance is True
    assert comparison.absolute_difference == pytest.approx(0.0)


def test_parse_existing_soluble_v_secretory_python_soplex_output_when_available() -> None:
    output = (
        ProjectPaths.discover().local_runs_dir
        / "pichia_python"
        / "v_secretory_smoke"
        / "python_V80_TEST_LEADER_mu0p10_soluble_v_secretory_target_reference_constraints.lp.float.out"
    )
    if not output.exists():
        pytest.skip("Local soluble V Python SoPlex output has not been generated on this machine.")

    summary = parse_soplex_output(output)

    assert summary.status == "problem is solved [optimal with unscaled violations]"
    assert summary.is_optimal is False
    assert summary.solution_type == "real"
    assert summary.diagnostic == "optimal_with_unscaled_violations"
    assert summary.objective_value == pytest.approx(0.0131529186)


def test_compare_existing_soluble_v_secretory_matlab_and_python_soplex_outputs_when_available() -> None:
    paths = ProjectPaths.discover()
    matlab_output = (
        paths.local_runs_dir
        / "pichia_python"
        / "v_secretory_baseline"
        / "Simulation_dilutionV80_TEST_LEADER_mu0p10_media4_v_secretory_noMisfoldEq_noRiboEq_PP.lp.float.out"
    )
    python_output = (
        paths.local_runs_dir
        / "pichia_python"
        / "v_secretory_smoke"
        / "python_V80_TEST_LEADER_mu0p10_soluble_v_secretory_target_reference_constraints.lp.float.out"
    )
    if not matlab_output.exists() or not python_output.exists():
        pytest.skip("Local soluble V MATLAB/Python SoPlex outputs have not been generated on this machine.")

    matlab_summary = parse_soplex_output(matlab_output)
    python_summary = parse_soplex_output(python_output)
    comparison = compare_soplex_objectives(matlab_output, python_output, tolerance=1e-12)

    assert matlab_summary.status == "problem is solved [optimal with unscaled violations]"
    assert python_summary.status == matlab_summary.status
    assert matlab_summary.diagnostic == "optimal_with_unscaled_violations"
    assert python_summary.diagnostic == matlab_summary.diagnostic
    assert comparison.both_optimal is False
    assert comparison.within_tolerance is True
    assert comparison.absolute_difference == pytest.approx(0.0)


def test_parse_existing_soluble_other_secretory_python_soplex_output_when_available() -> None:
    output = (
        ProjectPaths.discover().local_runs_dir
        / "pichia_python"
        / "other_secretory_smoke"
        / "python_M80_TEST_LEADER_mu0p10_soluble_other_secretory_target_reference_constraints.lp.float.out"
    )
    if not output.exists():
        pytest.skip("Local soluble Other Python SoPlex output has not been generated on this machine.")

    summary = parse_soplex_output(output)

    assert summary.status == "problem is solved [optimal]"
    assert summary.is_optimal is True
    assert summary.solution_type == "real"
    assert summary.diagnostic == "optimal"
    assert summary.objective_value == pytest.approx(0.0132299639)


def test_compare_existing_soluble_other_secretory_matlab_and_python_soplex_outputs_when_available() -> None:
    paths = ProjectPaths.discover()
    matlab_output = (
        paths.local_runs_dir
        / "pichia_python"
        / "other_secretory_baseline"
        / "Simulation_dilutionM80_TEST_LEADER_mu0p10_media4_other_secretory_noMisfoldEq_noRiboEq_PP.lp.float.out"
    )
    python_output = (
        paths.local_runs_dir
        / "pichia_python"
        / "other_secretory_smoke"
        / "python_M80_TEST_LEADER_mu0p10_soluble_other_secretory_target_reference_constraints.lp.float.out"
    )
    if not matlab_output.exists() or not python_output.exists():
        pytest.skip("Local soluble Other MATLAB/Python SoPlex outputs have not been generated on this machine.")

    matlab_summary = parse_soplex_output(matlab_output)
    python_summary = parse_soplex_output(python_output)
    comparison = compare_soplex_objectives(matlab_output, python_output, tolerance=1e-12)

    assert matlab_summary.status == "problem is solved [optimal]"
    assert python_summary.status == matlab_summary.status
    assert matlab_summary.diagnostic == "optimal"
    assert python_summary.diagnostic == matlab_summary.diagnostic
    assert comparison.both_optimal is True
    assert comparison.within_tolerance is True
    assert comparison.absolute_difference == pytest.approx(0.0)
