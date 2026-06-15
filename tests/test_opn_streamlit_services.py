from __future__ import annotations

import re
import shutil
from dataclasses import dataclass
from pathlib import Path

from app.adapters.process_runner import CommandResult
from app.core.paths import ProjectPaths
from app.services.opn import DEFAULT_OPN_CANDIDATE, OpnCandidateCatalog, OpnSimulationService


def test_opn_candidate_catalog_loads_demo_candidates() -> None:
    paths = ProjectPaths.discover()
    candidates = OpnCandidateCatalog(paths).list_candidates()
    by_id = {candidate.candidate_id: candidate for candidate in candidates}

    assert DEFAULT_OPN_CANDIDATE in by_id
    assert by_id[DEFAULT_OPN_CANDIDATE].category == "pichia_native_signal"
    assert by_id[DEFAULT_OPN_CANDIDATE].construct_length == 316


def test_opn_candidate_ranking_has_actionable_shortlist() -> None:
    paths = ProjectPaths.discover()
    rankings = OpnCandidateCatalog(paths).rank_candidates()
    by_id = {ranking.candidate_id: ranking for ranking in rankings}

    assert by_id["OPN_PPA_PASCHR3_0030"].rank == 1
    assert by_id["OPN_PPA_PASCHR3_0030"].recommendation == "首选实验候选"
    assert by_id["OPN_PPA_DDDK18"].rank == 2
    assert by_id["OPN_ALPHA_FULL_PROJECT"].experimental_role == "必须保留的对照"


def test_opn_construct_design_exports_full_sequences() -> None:
    paths = ProjectPaths.discover()
    designs = OpnCandidateCatalog(paths).construct_designs()
    by_id = {design.candidate_id: design for design in designs}

    assert by_id["OPN_PPA_PASCHR3_0030"].experimental_role == "优先构建 1"
    assert by_id["OPN_PPA_PASCHR3_0030"].mature_opn_length == 298
    assert by_id["OPN_PPA_PASCHR3_0030"].full_protein_sequence.startswith("MKFAISTLLIILQAAAVFAA")
    assert by_id["OPN_PPA_PASCHR3_0030"].contains_alpha_pro_region is False
    assert by_id["OPN_ALPHA_FULL_PROJECT"].contains_alpha_pro_region is True
    assert "Kex2" in by_id["OPN_ALPHA_FULL_PROJECT"].kex2_risk


@dataclass
class FakeMatlab:
    paths: ProjectPaths

    def run_batch(self, repo_root: Path, matlab_command: str, timeout_seconds=None) -> CommandResult:
        self.paths.opn_run_dir.mkdir(parents=True, exist_ok=True)
        lp = self.paths.opn_run_dir / (
            "Simulation_dilutionOPN_PPA_DDDK18_mu0p1_media4_"
            "ratio1em08_misfolddefault_noMisfoldEq_noRiboEq_PP.lp"
        )
        lp.write_text("Maximize\nobj: X1\nSubject To\nC1: 1 X1 = 0\nBounds\n0 <= X1 <= +infinity\nEnd\n", encoding="utf-8")
        return CommandResult(args=["matlab"], returncode=0, stdout=matlab_command, stderr="")


class FakeRunner:
    def run(self, args: list[str], cwd: Path, timeout_seconds=None) -> CommandResult:
        volume = args[args.index("-v") + 1]
        run_dir = Path(volume.split(":/work", 1)[0])
        shell_command = args[-1]
        match = re.search(r">\s+(.+)$", shell_command)
        output_name = match.group(1).strip("'\"") if match else "fake.lp.float.out"
        (run_dir / output_name).write_text(
            "SoPlex status       : problem is solved [optimal]\n"
            "Objective value     : -1.07457773e+00\n",
            encoding="utf-8",
        )
        return CommandResult(args=args, returncode=0, stdout="", stderr="")


def test_opn_simulation_service_returns_soplex_summary(tmp_path: Path) -> None:
    repo_paths = ProjectPaths.discover()
    paths = ProjectPaths(tmp_path)
    paths.opn_candidate_meta_csv.parent.mkdir(parents=True)
    shutil.copyfile(repo_paths.opn_candidate_meta_csv, paths.opn_candidate_meta_csv)
    service = OpnSimulationService(paths, FakeMatlab(paths), FakeRunner())

    result = service.run_candidate_smoke(DEFAULT_OPN_CANDIDATE, mu=0.10, production_ratio=1e-8, timeout_seconds=120)

    assert result.success is True
    assert result.candidate_id == DEFAULT_OPN_CANDIDATE
    assert result.objective_value == "-1.07457773e+00"
    assert result.lp_file is not None
    assert result.output_file is not None
