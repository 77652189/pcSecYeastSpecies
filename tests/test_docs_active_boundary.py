from __future__ import annotations

from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]

ACTIVE_DOCS = {
    "pichia_python_architecture.md",
    "pichia_python_next_development_slices_2026-06-26.md",
    "pichia_python_release_validation_2026-06-25.md",
}

ARCHIVED_MIGRATION_DOCS = {
    "opn_pichia_signal_peptide_candidates.md",
    "pichia_python_migration_strategy.md",
    "pichia_python_refactor_plan.md",
}


def test_docs_root_contains_only_reviewed_active_pichia_docs() -> None:
    docs_root = REPO_ROOT / "docs"
    root_markdown_files = {
        path.name for path in docs_root.glob("*.md") if path.is_file()
    }

    assert root_markdown_files == ACTIVE_DOCS


def test_historical_migration_docs_are_archived_not_active() -> None:
    docs_root = REPO_ROOT / "docs"
    archive_root = docs_root / "archive"
    root_markdown_files = {path.name for path in docs_root.glob("*.md")}
    archived_markdown_files = {path.name for path in archive_root.glob("*.md")}

    assert root_markdown_files.isdisjoint(ARCHIVED_MIGRATION_DOCS)
    assert ARCHIVED_MIGRATION_DOCS.issubset(archived_markdown_files)
