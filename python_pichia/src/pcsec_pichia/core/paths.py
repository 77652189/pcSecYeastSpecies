from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class ProjectPaths:
    repo_root: Path

    @classmethod
    def discover(cls, start: Path | None = None) -> "ProjectPaths":
        current = (start or Path(__file__)).resolve()
        candidates = [current] if current.is_dir() else [current.parent]
        candidates.extend(candidates[0].parents)
        for candidate in candidates:
            if (candidate / "README.md").exists() and (candidate / "Results").is_dir():
                return cls(candidate)
        raise FileNotFoundError("Could not locate pcSecYeastSpecies repository root.")

    @property
    def results_dir(self) -> Path:
        return self.repo_root / "Results"

    @property
    def models_dir(self) -> Path:
        return self.repo_root / "Model"

    @property
    def pichia_model_mat(self) -> Path:
        return self.models_dir / "pcSecPichia.mat"

    @property
    def pichia_enzymedata_dir(self) -> Path:
        return self.repo_root / "Enzymedata" / "pcSecPichia"

    @property
    def pichia_enzymedata_sec_mat(self) -> Path:
        return self.pichia_enzymedata_dir / "enzymedataSEC_PP.mat"

    @property
    def pichia_enzymedata_machine_mat(self) -> Path:
        return self.pichia_enzymedata_dir / "enzymedataMachine_PP.mat"

    @property
    def pichia_enzymedata_metabolic_mat(self) -> Path:
        return self.pichia_enzymedata_dir / "enzymedat_PP.mat"

    @property
    def pichia_enzymedata_dummy_er_mat(self) -> Path:
        return self.pichia_enzymedata_dir / "enzymedataDummyER_PP.mat"

    @property
    def local_runs_dir(self) -> Path:
        return self.repo_root / "local_runs"

    @property
    def phase1_png(self) -> Path:
        return self.local_runs_dir / "phase1" / "Fig1b_ModelComparisonPpa.png"

    @property
    def smoke_run_dir(self) -> Path:
        return self.local_runs_dir / "SCE_GLC_smoke"

    @property
    def local_preflight_script(self) -> Path:
        return self.repo_root / "local_preflight.ps1"

    @property
    def run_soplex_script(self) -> Path:
        return self.repo_root / "run_soplex_docker.ps1"

    @property
    def matlab_smoke_script(self) -> Path:
        return self.repo_root / "local_smoke_sce_glc.m"

    @property
    def opn_candidate_csv(self) -> Path:
        return self.repo_root / "Data" / "pcSecPichia" / "TargetProtein_OPN_candidates.csv"

    @property
    def opn_candidate_meta_csv(self) -> Path:
        return self.repo_root / "Data" / "pcSecPichia" / "TargetProtein_OPN_candidates_meta.csv"

    @property
    def pichia_aa_id_xlsx(self) -> Path:
        return self.repo_root / "Data" / "pcSecPichia" / "aa_id_PP.xlsx"

    @property
    def opn_run_dir(self) -> Path:
        return self.local_runs_dir / "OPN_PPA_glc_smoke"

    @property
    def opn_design_dir(self) -> Path:
        return self.local_runs_dir / "OPN_design"

    @property
    def opn_matlab_script(self) -> Path:
        return self.repo_root / "local_opn_pichia_glc.m"

    @property
    def pichia_clm_repo(self) -> Path:
        return self.repo_root.parent / "PichiaCLM"
