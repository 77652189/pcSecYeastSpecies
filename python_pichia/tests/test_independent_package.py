from __future__ import annotations

from pcsec_pichia.core.paths import ProjectPaths
from pcsec_pichia.core.target_inputs import LeaderCandidateInput, TargetProteinInput
from pcsec_pichia.core.target_protein_plan import build_target_protein_plan
from pcsec_pichia.services.target_simulation import PichiaTargetSimulationService


def test_pcsec_pichia_package_is_importable() -> None:
    """Verify the refactored service module can be imported and instantiated."""
    service = PichiaTargetSimulationService(ProjectPaths.discover())
    assert service is not None


def test_target_simulation_service_runs_opn_smoke(tmp_path) -> None:
    """OPN smoke simulation via the prototype-backed service (quick sanity check)."""
    paths = ProjectPaths.discover()
    target = TargetProteinInput(
        target_id="OPN_ALPHA_FULL_PROJECT",
        protein_name="OPN alpha full project",
        abbreviation="OPN_ALPHA_FULL_PROJECT",
        mature_sequence="M" + "A" * 297,
        through_er=1,
        signal_peptide=1,
        o_glycosylation_sites=20,
        localization="e",
        parameter_status="ready_for_model",
    )
    leader = LeaderCandidateInput(
        candidate_id="OPN_ALPHA_FULL_PROJECT_LEADER",
        leader_sequence="M",
        signal_peptide_sequence="M",
    )
    plan = build_target_protein_plan(target, leader)

    result = PichiaTargetSimulationService(paths).run_glucose_smoke(
        plan,
        output_dir=tmp_path,
        mu=0.10,
        production_ratio=1e-8,
        media_type=4,
    )

    assert result.success is True
    assert result.target_id == "OPN_ALPHA_FULL_PROJECT"
    assert result.candidate_id == "OPN_ALPHA_FULL_PROJECT_LEADER"
    assert result.objective_value is not None
    assert float(result.objective_value) > 0


def test_target_simulation_service_reports_failed_build_for_unsupported_target(tmp_path) -> None:
    """Unsupported target protein type should report failure gracefully."""
    paths = ProjectPaths.discover()
    target = TargetProteinInput(
        target_id="GPI_TEST",
        protein_name="GPI-anchored test",
        abbreviation="GPI_TEST",
        mature_sequence="A" * 50,
        through_er=1,
        signal_peptide=1,
        gpi_sites=2,
        localization="e",
        parameter_status="ready_for_model",
    )
    leader = LeaderCandidateInput(
        candidate_id="GPI_TEST_LEADER",
        leader_sequence="M",
        signal_peptide_sequence="M",
    )
    plan = build_target_protein_plan(target, leader)

    result = PichiaTargetSimulationService(paths).run_glucose_smoke(
        plan,
        output_dir=tmp_path,
        mu=0.10,
        production_ratio=1e-8,
        media_type=4,
    )

    assert result.success is False
