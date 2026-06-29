"""Legacy MATLAB-backed OPN smoke engine.

This facade preserves the historical app service path that shells out to
MATLAB and SoPlex. It is intentionally not part of the Python corrected draft
pipeline used by ``app.services.pichia_secretion_service``.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

from app.adapters.matlab import MatlabAdapter
from pcsec_pichia.adapters.process_runner import ProcessRunner
from pcsec_pichia.adapters.soplex_solver import DockerSoplexSolver, SoplexSolveResult
from pcsec_pichia.core.paths import ProjectPaths
from pcsec_pichia.engines.base import PichiaSimulationRequest, PichiaSimulationRunResult


def _matlab_string(value: str) -> str:
    return value.replace("'", "''")


@dataclass
class MatlabPichiaEngine:
    paths: ProjectPaths
    matlab: MatlabAdapter
    runner: ProcessRunner = field(default_factory=ProcessRunner)
    soplex_solver: DockerSoplexSolver | None = None
    engine_name: str = "matlab_pichia"

    def __post_init__(self) -> None:
        if self.soplex_solver is None:
            self.soplex_solver = DockerSoplexSolver(self.paths.repo_root, runner=self.runner)

    def run_target_smoke(self, request: PichiaSimulationRequest) -> PichiaSimulationRunResult:
        if request.target_id.upper() != "OPN":
            return PichiaSimulationRunResult(
                success=False,
                target_id=request.target_id,
                candidate_id=request.candidate_id,
                mu=request.mu,
                production_ratio=request.production_ratio,
                media_type=request.media_type,
                message=(
                    "当前 MATLAB engine 只封装了已有 OPN smoke 工作流；"
                    f"{request.target_id} 需要等 target 添加阶段完成后再运行。"
                ),
            )

        self.paths.opn_run_dir.mkdir(parents=True, exist_ok=True)
        matlab_command = self._build_opn_matlab_command(request)
        matlab_result = self.matlab.run_batch(self.paths.repo_root, matlab_command, timeout_seconds=220)
        if matlab_result.returncode != 0:
            return PichiaSimulationRunResult(
                success=False,
                target_id=request.target_id,
                candidate_id=request.candidate_id,
                mu=request.mu,
                production_ratio=request.production_ratio,
                media_type=request.media_type,
                message="MATLAB 生成 OPN 候选 LP 文件失败。",
                command_output=matlab_result.combined_output,
            )

        lp_file = self._latest_candidate_lp(request.candidate_id)
        if lp_file is None:
            return PichiaSimulationRunResult(
                success=False,
                target_id=request.target_id,
                candidate_id=request.candidate_id,
                mu=request.mu,
                production_ratio=request.production_ratio,
                media_type=request.media_type,
                message="MATLAB 已结束，但没有找到新生成的 OPN LP 文件。",
                command_output=matlab_result.combined_output,
            )

        output_file = lp_file.with_suffix(lp_file.suffix + ".float.out")
        soplex_result = self._run_soplex_float(lp_file, output_file, request.timeout_seconds)
        summary = soplex_result.summary
        success = soplex_result.success
        output = "\n".join([matlab_result.combined_output, soplex_result.command_result.combined_output])
        return PichiaSimulationRunResult(
            success=success,
            target_id=request.target_id,
            candidate_id=request.candidate_id,
            mu=request.mu,
            production_ratio=request.production_ratio,
            media_type=request.media_type,
            message="OPN 候选验证完成，SoPlex 返回 optimal。" if success else "OPN 候选求解未通过，请查看命令输出。",
            lp_file=lp_file,
            output_file=output_file,
            objective_value=summary.objective_text if summary else None,
            command_output=output,
        )

    def _build_opn_matlab_command(self, request: PichiaSimulationRequest) -> str:
        repo = self.paths.repo_root.as_posix()
        candidate = _matlab_string(request.candidate_id)
        return (
            f"cd('{repo}'); "
            "opts=struct('mediaType',4,'writeMisfoldingConstraints',false,'writeRibosomeConstraint',false); "
            f"local_opn_pichia_glc({request.mu:.4g},{request.production_ratio:.4g},[],opts,'{candidate}');"
        )

    def _latest_candidate_lp(self, candidate_id: str) -> Path | None:
        safe_candidate = re.sub(r"[^A-Za-z0-9]+", "_", candidate_id).strip("_")
        files = list(self.paths.opn_run_dir.glob(f"*{safe_candidate}*.lp"))
        return max(files, key=lambda path: path.stat().st_mtime) if files else None

    def _run_soplex_float(self, lp_file: Path, output_file: Path, timeout_seconds: int) -> SoplexSolveResult:
        assert self.soplex_solver is not None
        return self.soplex_solver.solve_float(lp_file, output_file, timeout_seconds=timeout_seconds)
