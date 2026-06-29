from __future__ import annotations

from pathlib import Path
from typing import Any

from app import ensure_python_pichia_on_path

ensure_python_pichia_on_path()

from pcsec_pichia.core.paths import ProjectPaths


def list_pichia_gene_options(
    query: str = "",
    limit: int = 50,
    paths: ProjectPaths | None = None,
) -> list[dict[str, Any]]:
    resolved_paths = paths or ProjectPaths.discover(Path(__file__))
    ensure_python_pichia_on_path()

    from pcsec_pichia.loading import load_pcsec_pichia_inputs
    from pcsec_pichia.screens import default_ko_genes
    from pcsec_pichia.services.gene_catalog import search_full_catalog

    inputs = load_pcsec_pichia_inputs(resolved_paths.repo_root)
    model = inputs.prepared_model
    max_rows = max(1, min(int(limit), 200))
    query_text = str(query or "").strip().lower()
    preferred = set(default_ko_genes(model, 20)) if not query_text else set()

    rows = [
        _gene_option_row(row)
        for row in search_full_catalog(query_text, model=model)
        if query_text or row.get("gene_id") in preferred or row.get("primary_category") == "分泌相关"
    ]
    rows.sort(key=lambda row: _gene_lookup_sort_key(row, query_text, preferred))
    return rows[:max_rows]


def list_curated_pichia_gene_catalog(query: str = "") -> list[Any]:
    """Return curated secretion-pathway KO/OE entries through the app facade."""
    ensure_python_pichia_on_path()

    from pcsec_pichia.services.gene_catalog import search_catalog

    return search_catalog(query)


def list_curated_pichia_gene_catalog_by_category() -> dict[str, list[Any]]:
    """Return curated secretion-pathway KO/OE entries grouped by category."""
    ensure_python_pichia_on_path()

    from pcsec_pichia.services.gene_catalog import get_catalog_by_category

    return get_catalog_by_category()


def load_pichia_full_model_gene_catalog() -> list[dict[str, object]]:
    """Return all model genes with lightweight reaction/process metadata."""
    ensure_python_pichia_on_path()

    from pcsec_pichia.services.gene_catalog import load_full_model_genes

    return load_full_model_genes()


def get_pichia_ko_genes_for_selection(selected_names: list[str]) -> list[str]:
    ensure_python_pichia_on_path()

    from pcsec_pichia.services.gene_catalog import get_ko_genes_for_selection

    return get_ko_genes_for_selection(selected_names)


def get_pichia_ko_reactions_for_selection(selected_names: list[str]) -> list[str]:
    ensure_python_pichia_on_path()

    from pcsec_pichia.services.gene_catalog import get_ko_reactions_for_selection

    return get_ko_reactions_for_selection(selected_names)


def get_pichia_oe_reactions_for_selection(selected_names: list[str]) -> list[str]:
    ensure_python_pichia_on_path()

    from pcsec_pichia.services.gene_catalog import get_oe_reactions_for_selection

    return get_oe_reactions_for_selection(selected_names)


def _gene_option_row(row: dict[str, object]) -> dict[str, Any]:
    sample_reactions = [str(item) for item in row.get("sample_reactions") or []]
    processes = [item.strip() for item in str(row.get("processes") or "").split(",") if item.strip()]
    return {
        "gene_id": str(row.get("gene_id") or ""),
        "reaction_count": int(row.get("n_reactions") or 0),
        "reactions_preview": sample_reactions[:5],
        "secretory_processes": processes,
        "gr_rule_preview": "",
        "suggested_use": "KO / OE gene proxy",
        "primary_category": row.get("primary_category"),
    }


def _gene_lookup_sort_key(row: dict[str, Any], query: str, preferred: set[str]) -> tuple[int, int, str]:
    gene_id = str(row.get("gene_id") or "")
    processes = " ".join(str(item) for item in row.get("secretory_processes") or [])
    exact_rank = 0 if query and gene_id.lower().startswith(query) else 1
    preferred_rank = 0 if gene_id in preferred or _is_secretory_lookup_process(processes) else 1
    return (exact_rank, preferred_rank, gene_id)


def _is_secretory_lookup_process(process: str) -> bool:
    return any(token in process for token in ("翻译", "ER", "DSB", "糖基化", "错误折叠", "Golgi", "分泌"))


__all__ = [
    "get_pichia_ko_genes_for_selection",
    "get_pichia_ko_reactions_for_selection",
    "get_pichia_oe_reactions_for_selection",
    "list_curated_pichia_gene_catalog",
    "list_curated_pichia_gene_catalog_by_category",
    "list_pichia_gene_options",
    "load_pichia_full_model_gene_catalog",
]
