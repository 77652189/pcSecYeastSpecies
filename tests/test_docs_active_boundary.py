from __future__ import annotations

from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]

ACTIVE_DOCS = {
    "data_and_results_policy.md",
    "pichia_current_architecture_and_requirements.md",
    "pichia_next_plan.md",
    "README.md",
}

NON_ACTIVE_REFERENCE_DOCS = {
    "cobrapy_phase0_baseline_assessment_2026-07-06.md",
    "cobrapy_phase3_installed_shadow_validation_2026-07-06.md",
    "opn_pichia_signal_peptide_candidates.md",
    "pichia_cobrapy_import_qa_shadow_plan.md",
    "pichia_homology_crosswalk_architecture.md",
    "pichia_ko_oe_genome_screen_design_2026-07-02.md",
    "pichia_medium_mixed_carbon_objective_plan_2026-06-30.md",
    "pichia_online_external_reference_architecture.md",
    "pichia_python_hlf_design_decisions.md",
    "pichia_python_hlf_project_710_alignment_status_2026-06-26.md",
    "pichia_sce_homology_feasibility_20260708.md",
}

DELETED_OBSOLETE_MIGRATION_DOCS = {
    "pichia_python_architecture.md",
    "pichia_python_next_development_slices_2026-06-26.md",
    "pichia_python_release_validation_2026-06-25.md",
    "pichia_python_migration_strategy.md",
    "pichia_python_refactor_plan.md",
    "migration_progress.md",
}


def test_docs_root_contains_only_reviewed_active_pichia_docs() -> None:
    docs_root = REPO_ROOT / "docs"
    root_markdown_files = {
        path.name for path in docs_root.glob("*.md") if path.is_file()
    }

    assert root_markdown_files == ACTIVE_DOCS


def test_completed_reference_docs_are_not_active_entries() -> None:
    docs_root = REPO_ROOT / "docs"
    root_markdown_files = {path.name for path in docs_root.glob("*.md")}

    assert root_markdown_files.isdisjoint(NON_ACTIVE_REFERENCE_DOCS)


def test_obsolete_migration_plans_are_deleted_not_kept_as_active_debt() -> None:
    docs_root = REPO_ROOT / "docs"
    archive_root = docs_root / "archive"
    python_pichia_docs_root = REPO_ROOT / "python_pichia" / "docs"
    all_doc_names = {
        path.name
        for root in (docs_root, archive_root, python_pichia_docs_root)
        for path in root.glob("*.md")
    }

    assert all_doc_names.isdisjoint(DELETED_OBSOLETE_MIGRATION_DOCS)


def test_next_plan_records_the_current_phase_2_round() -> None:
    text = (REPO_ROOT / "docs" / "pichia_next_plan.md").read_text(encoding="utf-8")
    current_state = text.split("```yaml", 1)[1].split("```", 1)[0]

    assert "current_phase: phase_2_gene_level_oe" in current_state
    assert "current_round: round_6_hlf_opn_acceptance" in current_state
    assert "round_status: ready" in current_state
    assert "从 Phase 1 Round 0 开始" not in text
    assert "从 Phase 2 Round 6 开始" in text
    assert "本阶段已完成并通过端到端验收" in text
