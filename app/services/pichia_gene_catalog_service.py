from __future__ import annotations

import re
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
    from pcsec_pichia.screens import classify_secretory_process, default_ko_genes

    inputs = load_pcsec_pichia_inputs(resolved_paths.repo_root)
    model = inputs.prepared_model
    max_rows = max(1, min(int(limit), 200))
    query_text = str(query or "").strip().lower()
    preferred = set(default_ko_genes(model, 20)) if not query_text else set()
    gene_rows: dict[str, dict[str, Any]] = {}

    for reaction_index, reaction_id in enumerate(model.rxns):
        rule = model.rules[reaction_index] if reaction_index < len(model.rules) else ""
        gr_rule = model.gr_rules[reaction_index] if reaction_index < len(model.gr_rules) else ""
        process = _secretory_process_label(classify_secretory_process(reaction_id))
        for gene_id in _genes_for_reaction(model, str(rule or ""), str(gr_rule or "")):
            if query_text and not _gene_lookup_matches(query_text, gene_id, reaction_id, gr_rule, process):
                continue
            if not query_text and gene_id not in preferred and not _is_secretory_lookup_process(process):
                continue
            row = gene_rows.setdefault(
                gene_id,
                {
                    "gene_id": gene_id,
                    "reaction_count": 0,
                    "reactions_preview": [],
                    "secretory_processes": [],
                    "gr_rule_preview": "",
                    "suggested_use": "KO / OE gene proxy",
                },
            )
            row["reaction_count"] = int(row["reaction_count"]) + 1
            if len(row["reactions_preview"]) < 5:
                row["reactions_preview"].append(reaction_id)
            if process not in row["secretory_processes"]:
                row["secretory_processes"].append(process)
            if not row["gr_rule_preview"] and gr_rule:
                row["gr_rule_preview"] = str(gr_rule)[:160]

    rows = list(gene_rows.values())
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


def _secretory_process_label(process_code: str) -> str:
    labels = {
        "ribosome": "翻译",
        "proteasome_degradation": "蛋白降解",
        "disulfide_folding": "ER 折叠 / DSB",
        "n_glycan_processing": "N-糖基化 NG",
        "o_glycan_processing": "O-糖基化 OG",
        "chaperone_folding": "ER 折叠 / 分子伴侣",
        "erad_misfolding": "错误折叠 / ERAD",
        "er_translocation": "ER 转运",
        "er_to_golgi_transport": "ER 到 Golgi 转运",
        "golgi_surface_transport": "Golgi 到胞外运输",
        "secretory_capacity": "分泌容量",
        "metabolic_or_other": "代谢或其它反应",
        "unknown": "未解析",
    }
    return labels.get(process_code, process_code)


def _genes_for_reaction(model: Any, rule: str, gr_rule: str) -> list[str]:
    genes: list[str] = []
    for match in re.finditer(r"x\((\d+)\)", rule):
        gene_index = int(match.group(1)) - 1
        if 0 <= gene_index < len(model.genes):
            gene_id = str(model.genes[gene_index])
            if gene_id not in genes:
                genes.append(gene_id)
    for gene_id in re.findall(r"PAS_[A-Za-z0-9_\-]+", gr_rule):
        if gene_id in model.gene_index and gene_id not in genes:
            genes.append(gene_id)
    return genes


def _gene_lookup_matches(query: str, gene_id: str, reaction_id: str, gr_rule: str, process: str) -> bool:
    haystack = " ".join([gene_id, reaction_id, gr_rule, process]).lower()
    return all(token in haystack for token in re.split(r"[\s,;，]+", query) if token)


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
