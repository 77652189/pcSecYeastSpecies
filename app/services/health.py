from __future__ import annotations

import shutil
from dataclasses import dataclass

from app.adapters.powershell import PowerShellAdapter
from app.core.models import HealthItem, HealthReport
from pcsec_pichia.core.paths import ProjectPaths


@dataclass
class HealthService:
    paths: ProjectPaths
    powershell: PowerShellAdapter

    def check(self) -> HealthReport:
        items = [
            self._file_item("项目说明 README", self.paths.repo_root / "README.md"),
            self._file_item("酿酒酵母模型 pcSecYeast.mat", self.paths.models_dir / "pcSecYeast.mat"),
            self._file_item("结果目录 Results", self.paths.results_dir),
            self._file_item("本地运行目录 local_runs", self.paths.local_runs_dir),
            self._command_item("Docker", "docker"),
            self._command_item("PowerShell", "powershell"),
        ]
        preflight_output = ""
        if self.paths.local_preflight_script.exists():
            result = self.powershell.run_script(self.paths.local_preflight_script, self.paths.repo_root, timeout_seconds=60)
            preflight_output = result.combined_output
            items.append(
                HealthItem(
                    name="本地依赖预检脚本",
                    status="ok" if result.returncode == 0 else "error",
                    detail="预检脚本已执行" if result.returncode == 0 else "预检脚本返回错误，请查看输出",
                )
            )
        return HealthReport(items=items, preflight_output=preflight_output)

    def _file_item(self, name: str, path) -> HealthItem:
        return HealthItem(
            name=name,
            status="ok" if path.exists() else "missing",
            detail=str(path),
        )

    def _command_item(self, name: str, command: str) -> HealthItem:
        found = shutil.which(command)
        return HealthItem(name=name, status="ok" if found else "missing", detail=found or "未在 PATH 中找到")
