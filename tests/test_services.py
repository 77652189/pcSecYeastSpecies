from __future__ import annotations

import ast
from dataclasses import dataclass
from pathlib import Path

from pcsec_pichia.adapters.process_runner import CommandResult
from app.adapters.soplex_parser import parse_soplex_output
from pcsec_pichia.core.paths import ProjectPaths
from app.services.health import HealthService
from app.services.simulation import parse_objective_from_text


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_app_no_longer_imports_deleted_process_runner_adapter() -> None:
    deleted_adapter_imports: list[str] = []
    for path in (REPO_ROOT / "app").rglob("*.py"):
        module_ast = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(module_ast):
            imported: list[str] = []
            if isinstance(node, ast.ImportFrom) and node.module:
                imported.append(node.module)
            elif isinstance(node, ast.Import):
                imported.extend(alias.name for alias in node.names)
            for module_name in imported:
                if module_name == "app.adapters.process_runner":
                    deleted_adapter_imports.append(str(path.relative_to(REPO_ROOT)))

    assert deleted_adapter_imports == []


def test_soplex_parser_reads_optimal_objective() -> None:
    text = """
    SoPlex status       : problem is solved [optimal]
    Objective value     : -1.05755453e+00
    """

    summary = parse_soplex_output(text)

    assert summary.optimal is True
    assert summary.objective_value == "-1.05755453e+00"
    assert summary.status_line == "problem is solved [optimal]"
    assert summary.diagnostic == "optimal"


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
