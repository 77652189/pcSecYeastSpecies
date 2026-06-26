from __future__ import annotations

from dataclasses import dataclass

from app.adapters.soplex_parser import parse_soplex_file
from pcsec_pichia.core.paths import ProjectPaths


@dataclass
class RunLogService:
    paths: ProjectPaths

    def recent_files(self, limit: int = 20) -> list[dict[str, object]]:
        if not self.paths.local_runs_dir.exists():
            return []
        files = sorted(
            [path for path in self.paths.local_runs_dir.rglob("*") if path.is_file()],
            key=lambda p: p.stat().st_mtime,
            reverse=True,
        )
        return [
            {
                "文件": str(path.relative_to(self.paths.repo_root)),
                "大小KB": round(path.stat().st_size / 1024, 1),
                "修改时间": path.stat().st_mtime,
            }
            for path in files[:limit]
        ]

    def latest_soplex_summary(self):
        outputs = list(self.paths.local_runs_dir.rglob("*.lp.out")) if self.paths.local_runs_dir.exists() else []
        if not outputs:
            return None, None
        latest = max(outputs, key=lambda p: p.stat().st_mtime)
        return latest, parse_soplex_file(latest)
