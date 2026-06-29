from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import pytest

from pcsec_pichia.adapters.process_runner import CommandResult
from pcsec_pichia.adapters.soplex_solver import DockerSoplexSolver, write_soplex_solve_result_json


@dataclass
class FakeRunner:
    returncode: int = 0
    last_args: list[str] | None = None
    last_cwd: Path | None = None
    last_timeout_seconds: int | None = None

    def run(self, args: list[str], cwd: Path, timeout_seconds=None) -> CommandResult:
        self.last_args = args
        self.last_cwd = cwd
        self.last_timeout_seconds = timeout_seconds
        volume = args[args.index("-v") + 1]
        run_dir = Path(volume.split(":/work", 1)[0])
        shell_command = args[-1]
        output_match = re.search(r">\s+(.+)$", shell_command)
        output_name = output_match.group(1).strip("'\"") if output_match else "solve.lp.float.out"
        (run_dir / output_name).write_text(
            "SoPlex status       : problem is solved [optimal]\n"
            "Solution (real)     : \n"
            "  Objective value   : -1.07282263e+00\n"
            "Iterations          : 22698\n",
            encoding="utf-8",
        )
        return CommandResult(args=args, returncode=self.returncode, stdout="docker ok", stderr="")


def test_docker_soplex_solver_builds_float_command_with_expected_flags(tmp_path: Path) -> None:
    lp = tmp_path / "case.lp"
    out = tmp_path / "case.lp.float.out"
    lp.write_text("Maximize\nobj: X1\nSubject To\nC1: X1 = 0\nBounds\n0 <= X1 <= 1\nEnd\n", encoding="utf-8")

    args = DockerSoplexSolver(tmp_path).float_docker_args(lp, out, timeout_seconds=120)
    shell_command = args[-1]

    assert args[:3] == ["docker", "run", "--rm"]
    assert args[args.index("-w") + 1] == "/work"
    assert "soplex -s0 -g5 -t120 -q" in shell_command
    assert "--readmode=0 --solvemode=0" in shell_command
    assert "--real:fpfeastol=1e-3 --real:fpopttol=1e-3" in shell_command
    assert "case.lp > case.lp.float.out" in shell_command


def test_docker_soplex_solver_runs_and_parses_summary(tmp_path: Path) -> None:
    runner = FakeRunner()
    lp = tmp_path / "case.lp"
    lp.write_text("Maximize\nobj: X1\nSubject To\nC1: X1 = 0\nBounds\n0 <= X1 <= 1\nEnd\n", encoding="utf-8")
    output = tmp_path / "case.lp.float.out"

    result = DockerSoplexSolver(tmp_path, runner=runner).solve_float(lp, output, timeout_seconds=120)

    assert result.success is True
    assert result.output_path == output.resolve()
    assert result.command_result.returncode == 0
    assert result.summary is not None
    assert result.summary.is_optimal is True
    assert result.summary.objective_value == pytest.approx(-1.07282263)
    assert result.summary.objective_text == "-1.07282263e+00"
    assert result.summary.solution_type == "real"
    assert runner.last_cwd == tmp_path
    assert runner.last_timeout_seconds == 180

    saved = write_soplex_solve_result_json(result, tmp_path / "solve_summary.json")
    assert saved.exists()
    assert '"success": true' in saved.read_text(encoding="utf-8")


def test_docker_soplex_solver_requires_same_lp_and_output_directory(tmp_path: Path) -> None:
    lp = tmp_path / "case.lp"
    lp.write_text("", encoding="utf-8")
    other_dir = tmp_path / "other"
    other_dir.mkdir()

    with pytest.raises(ValueError, match="same directory"):
        DockerSoplexSolver(tmp_path).solve_float(lp, other_dir / "case.lp.float.out")
