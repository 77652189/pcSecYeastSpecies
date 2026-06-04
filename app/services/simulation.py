from __future__ import annotations

import re
import threading
from dataclasses import dataclass, field
from pathlib import Path

from app.adapters.matlab import MatlabAdapter
from app.adapters.powershell import PowerShellAdapter
from app.adapters.soplex_parser import parse_soplex_file
from app.core.models import SimulationResult
from app.core.paths import ProjectPaths


MU_MIN = 0.01
MU_MAX = 0.44


@dataclass
class SimulationService:
    paths: ProjectPaths
    matlab: MatlabAdapter
    powershell: PowerShellAdapter
    _lock: threading.Lock = field(default_factory=threading.Lock)

    def run_sce_glucose_smoke(self, mu: float = 0.10, timeout_seconds: int = 300) -> SimulationResult:
        if not MU_MIN <= mu <= MU_MAX:
            raise ValueError("生长速率 mu 必须在 0.01 到 0.44 h^-1 之间。")
        if not self._lock.acquire(blocking=False):
            return SimulationResult(success=False, mu=mu, message="已有仿真任务正在运行，请稍后再试。")
        try:
            mu_command = f"{mu:.2f}"
            matlab_command = (
                f"cd('{self.paths.repo_root.as_posix()}'); "
                "startup_pcsec_local('yeast'); "
                f"local_smoke_sce_glc({mu_command});"
            )
            matlab_result = self.matlab.run_batch(self.paths.repo_root, matlab_command, timeout_seconds=180)
            if matlab_result.returncode != 0:
                return SimulationResult(
                    success=False,
                    mu=mu,
                    message="MATLAB 生成线性规划 LP 文件失败。",
                    command_output=matlab_result.combined_output,
                )

            soplex_result = self.powershell.run_script(
                self.paths.run_soplex_script,
                self.paths.repo_root,
                args=["-TimeoutSeconds", str(timeout_seconds)],
                timeout_seconds=timeout_seconds + 60,
            )
            output = "\n".join([matlab_result.combined_output, soplex_result.combined_output])
            out_file = self._latest_soplex_output()
            summary = parse_soplex_file(out_file) if out_file else None
            success = soplex_result.returncode == 0 and bool(summary and summary.optimal)
            return SimulationResult(
                success=success,
                mu=mu,
                message="仿真完成，求解器 SoPlex 已返回可用结果。" if success else "SoPlex 求解未通过，请查看命令输出。",
                lp_file=self._latest_lp_file(),
                output_file=out_file,
                objective_value=summary.objective_value if summary else None,
                command_output=output,
            )
        finally:
            self._lock.release()

    def _latest_lp_file(self) -> Path | None:
        files = list(self.paths.smoke_run_dir.glob("*.lp"))
        return max(files, key=lambda p: p.stat().st_mtime) if files else None

    def _latest_soplex_output(self) -> Path | None:
        files = list(self.paths.smoke_run_dir.glob("*.lp.out"))
        return max(files, key=lambda p: p.stat().st_mtime) if files else None


def parse_objective_from_text(text: str) -> str | None:
    matches = re.findall(r"Objective value\s*:?\s*(.+)", text)
    return matches[-1].strip() if matches else None
