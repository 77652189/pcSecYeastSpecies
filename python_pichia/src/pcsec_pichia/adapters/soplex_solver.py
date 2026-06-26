from __future__ import annotations

import json
import shlex
from dataclasses import dataclass, field
from pathlib import Path

from pcsec_pichia.adapters.process_runner import CommandResult, ProcessRunner
from pcsec_pichia.adapters.soplex_output import SoplexOutputSummary, parse_soplex_output


@dataclass(frozen=True)
class SoplexSolveResult:
    lp_path: Path
    output_path: Path
    command_result: CommandResult
    summary: SoplexOutputSummary | None

    @property
    def success(self) -> bool:
        return self.command_result.returncode == 0 and bool(self.summary and self.summary.is_optimal)

    def to_dict(self) -> dict[str, object]:
        return {
            "lp_path": str(self.lp_path),
            "output_path": str(self.output_path),
            "returncode": self.command_result.returncode,
            "success": self.success,
            "summary": self.summary.to_dict() if self.summary else None,
            "stdout": self.command_result.stdout,
            "stderr": self.command_result.stderr,
        }


@dataclass(frozen=True)
class DockerSoplexSolver:
    repo_root: Path
    runner: ProcessRunner = field(default_factory=ProcessRunner)
    image: str = "pcsec-soplex:24.04"

    def solve_float(
        self,
        lp_path: Path,
        output_path: Path | None = None,
        timeout_seconds: int = 300,
    ) -> SoplexSolveResult:
        lp_path = lp_path.resolve()
        output_path = (output_path or lp_path.with_suffix(lp_path.suffix + ".float.out")).resolve()
        if lp_path.parent != output_path.parent:
            raise ValueError("LP file and SoPlex output must be in the same directory for the Docker workdir mount.")
        command_result = self.runner.run(
            self.float_docker_args(lp_path, output_path, timeout_seconds=timeout_seconds),
            cwd=self.repo_root,
            timeout_seconds=timeout_seconds + 60,
        )
        summary = parse_soplex_output(output_path) if output_path.exists() else None
        return SoplexSolveResult(
            lp_path=lp_path,
            output_path=output_path,
            command_result=command_result,
            summary=summary,
        )

    def float_docker_args(self, lp_path: Path, output_path: Path, timeout_seconds: int = 300) -> list[str]:
        lp_path = lp_path.resolve()
        output_path = output_path.resolve()
        if lp_path.parent != output_path.parent:
            raise ValueError("LP file and SoPlex output must be in the same directory for the Docker workdir mount.")
        shell_command = self.float_shell_command(
            lp_name=lp_path.name,
            output_name=output_path.name,
            timeout_seconds=timeout_seconds,
        )
        return [
            "docker",
            "run",
            "--rm",
            "-v",
            f"{lp_path.parent}:/work",
            "-w",
            "/work",
            self.image,
            "sh",
            "-lc",
            shell_command,
        ]

    @staticmethod
    def float_shell_command(lp_name: str, output_name: str, timeout_seconds: int = 300) -> str:
        return (
            f"timeout {int(timeout_seconds)} soplex -s0 -g5 -t{int(timeout_seconds)} -q "
            "--readmode=0 --solvemode=0 --real:fpfeastol=1e-3 --real:fpopttol=1e-3 "
            f"{shlex.quote(lp_name)} > {shlex.quote(output_name)}"
        )


def write_soplex_solve_result_json(result: SoplexSolveResult, output_path: Path) -> Path:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(result.to_dict(), ensure_ascii=False, indent=2), encoding="utf-8")
    return output_path
