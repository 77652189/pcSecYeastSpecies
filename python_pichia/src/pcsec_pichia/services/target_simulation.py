from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from pcsec_pichia.engines.base import PichiaSimulationRunResult


@dataclass
class PichiaTargetSimulationService:
    """pcSecPichia simulation service backed by the validated probe prototype.

    This replaces the legacy ``PythonPichiaEngine``-based dispatch that
    generated route-specific LP files. The new service uses the prototype's
    ``CobraModel`` workflow and solves directly with SciPy HiGHS — no external
    solvers or LP files are needed for production simulation runs.
    """

    paths: Any  # pcsec_pichia.core.paths.ProjectPaths or compatible duck type

    def run_glucose_smoke(
        self,
        plan: Any,
        output_dir: Path | None = None,
        mu: float = 0.10,
        production_ratio: float = 1e-8,
        media_type: int = 4,
        timeout_seconds: int = 300,
        compatibility_mode: str = "corrected",
    ) -> PichiaSimulationRunResult:
        """Run a glucose-condition target-protein smoke using the prototype solver.

        The prototype builds the target-protein model (CobraModel), applies
        pcSec constraints, and solves with SciPy HiGHS. No LP file is generated.
        The *production_ratio* parameter is accepted for API compatibility but
        the prototype maximises secretion directly (it does not fix the exchange).
        """
        import numpy as np

        from pcsec_pichia.core.modes import CompatibilityMode
        from pcsec_pichia.core.target_inputs import TargetProteinInput, LeaderCandidateInput
        from pcsec_pichia.targets import target_spec_from_input
        from pcsec_pichia.loading import load_pcsec_pichia_inputs
        from pcsec_pichia.simulation import solve_secretion_capacity

        # Build a TargetSpec from the build plan
        target = TargetProteinInput(
            target_id=plan.target_id,
            protein_name=plan.protein_id,
            abbreviation=plan.protein_id,
            mature_sequence=plan.mature_sequence or plan.full_sequence,
            through_er=int(plan.through_er),
            localization=str(plan.localization),
            disulfide_sites=int(plan.disulfide_sites),
            n_glycosylation_sites=int(plan.n_glycosylation_sites),
            o_glycosylation_sites=int(plan.o_glycosylation_sites),
            transmembrane=int(plan.transmembrane),
            gpi_sites=int(plan.gpi_sites),
            cotranslation=int(plan.cotranslation),
            signal_peptide=int(bool(plan.signal_peptide_sequence)),
            parameter_status="ready_for_model",
        )
        leader = LeaderCandidateInput(
            candidate_id=plan.protein_id,
            leader_sequence=str(plan.leader_sequence or ""),
            signal_peptide_sequence=str(plan.signal_peptide_sequence or ""),
        )
        resolved_target = target_spec_from_input(target, leader)

        # Load model inputs via the prototype-backed loading module
        root = getattr(self.paths, "repo_root", Path.cwd())
        inputs = load_pcsec_pichia_inputs(
            root=root,
            media_type=media_type,
            compatibility_mode=compatibility_mode,  # type: ignore[arg-type]
        )

        # Solve secretion capacity (maximises target exchange at fixed growth)
        result = solve_secretion_capacity(
            model=inputs.prepared_model,
            target=resolved_target,
            amino_acids=inputs.amino_acids,
            metabolic=inputs.metabolic,
            secretory=inputs.secretory,
            combined=inputs.combined,
            growth_rate=mu,
            write_ribosome_translation_constraint=False,
            write_misfolding_constraints=False,
        )

        objective_text = (
            f"{result.objective_value:.9g}"
            if result.objective_value is not None
            else None
        )

        return PichiaSimulationRunResult(
            success=bool(result.success),
            target_id=plan.target_id,
            candidate_id=plan.protein_id,
            mu=mu,
            production_ratio=production_ratio,
            media_type=media_type,
            message=(
                "Python pcSecPichia engine (prototype) solved successfully."
                if result.success
                else f"Python pcSecPichia engine failed: {result.message}"
            ),
            objective_value=str(objective_text) if objective_text is not None else None,
            result_status=str(result.result_status),  # type: ignore[arg-type]
            matlab_alignment_status=str(result.matlab_alignment_status),
            constraint_counts={str(k): int(v) for k, v in result.constraint_counts.items()},
            command_output="",
        )


__all__ = [
    "PichiaTargetSimulationService",
]
