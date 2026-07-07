from __future__ import annotations

import re

from pcsec_pichia.homology.cache_schema import CatalogHomologyQuery
from pcsec_pichia.services.gene_catalog import SECRETION_GENE_CATALOG, SecretionGeneEntry


SCE_SYMBOL_ALIASES: dict[str, tuple[str, ...]] = {
    "DOA10": ("SSM4",),
}


def secretion_catalog_sce_queries() -> tuple[CatalogHomologyQuery, ...]:
    """Build SCE query symbols and aliases from the curated secretion catalog."""

    queries: list[CatalogHomologyQuery] = []
    seen: set[tuple[str, str]] = set()
    for entry in SECRETION_GENE_CATALOG:
        for symbol in expand_catalog_aliases(entry):
            normalized = normalize_sce_symbol(symbol)
            if not normalized:
                continue
            key = (entry.common_name, normalized)
            if key in seen:
                continue
            seen.add(key)
            aliases = SCE_SYMBOL_ALIASES.get(normalized, ())
            queries.append(
                CatalogHomologyQuery(
                    internal_common_name=entry.common_name,
                    query_symbol=normalized,
                    aliases=aliases,
                    internal_gene_id=entry.gene_id,
                )
            )
    return tuple(queries)


def normalize_sce_symbol(symbol: str) -> str:
    """Normalize common-name casing and separators without guessing biology."""

    value = symbol.strip().upper()
    value = re.sub(r"（.*?）|\(.*?\)", "", value)
    value = value.replace("P", "") if value.endswith("P") and len(value) > 3 else value
    return re.sub(r"[^A-Z0-9-]", "", value)


def expand_catalog_aliases(entry: SecretionGeneEntry) -> tuple[str, ...]:
    """Extract plausible SCE symbols from a curated catalog entry."""

    raw = entry.common_name
    raw = re.sub(r"（.*?）|\(.*?\)", "", raw)
    parts = re.split(r"[/,+\s]+", raw)
    aliases = [normalize_sce_symbol(part) for part in parts if part.strip()]
    return tuple(alias for alias in aliases if alias)
