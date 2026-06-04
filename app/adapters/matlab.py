from __future__ import annotations

import os
from pathlib import Path

from app.adapters.process_runner import CommandResult, ProcessRunner


class MatlabAdapter:
    def __init__(self, runner: ProcessRunner | None = None) -> None:
        self.runner = runner or ProcessRunner()

    def find_executable(self) -> Path | None:
        explicit = os.environ.get("MATLAB_EXE")
        if explicit and Path(explicit).exists():
            return Path(explicit)
        for root in [Path(r"C:\Program Files\MATLAB"), Path(r"C:\Program Files (x86)\MATLAB")]:
            if not root.exists():
                continue
            candidates = sorted(root.glob(r"*\bin\matlab.exe"), reverse=True)
            if candidates:
                return candidates[0]
        return None

    def run_batch(self, repo_root: Path, matlab_command: str, timeout_seconds: int | None = None) -> CommandResult:
        matlab = self.find_executable()
        if matlab is None:
            raise FileNotFoundError("MATLAB executable not found. Set MATLAB_EXE or install MATLAB R2020b+.")
        return self.runner.run([str(matlab), "-batch", matlab_command], cwd=repo_root, timeout_seconds=timeout_seconds)
