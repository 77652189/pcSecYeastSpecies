from __future__ import annotations

from pathlib import Path

from app.adapters.process_runner import CommandResult, ProcessRunner


class PowerShellAdapter:
    def __init__(self, runner: ProcessRunner | None = None) -> None:
        self.runner = runner or ProcessRunner()

    def run_script(
        self,
        script: Path,
        cwd: Path,
        args: list[str] | None = None,
        timeout_seconds: int | None = None,
    ) -> CommandResult:
        command = [
            "powershell",
            "-ExecutionPolicy",
            "Bypass",
            "-File",
            str(script),
            *(args or []),
        ]
        return self.runner.run(command, cwd=cwd, timeout_seconds=timeout_seconds)
