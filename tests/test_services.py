from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from app.adapters.process_runner import CommandResult
from app.adapters.soplex_parser import parse_soplex_output
from app.core.paths import ProjectPaths
from app.services.health import HealthService
from app.services.simulation import parse_objective_from_text


def test_soplex_parser_reads_optimal_objective() -> None:
    text = """
    SoPlex status       : problem is solved [optimal]
    Objective value     : -1.05755453e+00
    """

    summary = parse_soplex_output(text)

    assert summary.optimal is True
    assert summary.objective_value == "-1.05755453e+00"
    assert summary.status_line == "problem is solved [optimal]"


def test_objective_parser_uses_last_objective() -> None:
    assert parse_objective_from_text("Objective value: 1\nObjective value   : 2") == "2"


@dataclass
class FakePowerShell:
    def run_script(self, script: Path, cwd: Path, args=None, timeout_seconds=None) -> CommandResult:
        return CommandResult(args=["powershell"], returncode=0, stdout="preflight ok", stderr="")


def test_health_service_includes_preflight_result() -> None:
    paths = ProjectPaths.discover()
    report = HealthService(paths, FakePowerShell()).check()

    assert any(item.name == "本地依赖预检脚本" and item.status == "ok" for item in report.items)
    assert report.preflight_output == "preflight ok"
