from __future__ import annotations

from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]

ACTIVE_DOCS = {
    "cobrapy_phase0_baseline_assessment_2026-07-06.md",
    "data_and_results_policy.md",
    "pichia_current_architecture_and_requirements.md",
    "pichia_next_plan.md",
    "pichia_ko_oe_genome_screen_design.md",
}

ARCHIVED_REFERENCE_DOCS = {
    "opn_pichia_signal_peptide_candidates.md",
    "pichia_medium_mixed_carbon_objective_plan_2026-06-30.md",
    "pichia_python_hlf_design_decisions.md",
    "pichia_python_hlf_project_710_alignment_status_2026-06-26.md",
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


def test_reference_docs_are_archived_not_active() -> None:
    docs_root = REPO_ROOT / "docs"
    archive_root = docs_root / "archive"
    root_markdown_files = {path.name for path in docs_root.glob("*.md")}
    archived_markdown_files = {path.name for path in archive_root.glob("*.md")}

    assert root_markdown_files.isdisjoint(ARCHIVED_REFERENCE_DOCS)
    assert ARCHIVED_REFERENCE_DOCS.issubset(archived_markdown_files)


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
