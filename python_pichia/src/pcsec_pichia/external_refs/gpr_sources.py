from __future__ import annotations

from pathlib import Path
from typing import Iterable
from xml.etree import ElementTree

from pcsec_pichia.external_refs.cache_io import DEFAULT_RECORDS_FILENAME, load_external_reference_cache
from pcsec_pichia.external_refs.clients import ExternalFetchConfig
from pcsec_pichia.external_refs.schema import (
    ExternalReactionAssociation,
    ExternalReferenceProvenance,
    sha256_text,
    utc_now_iso,
)


def fetch_external_model_reaction_associations(
    *,
    source_database: str,
    model_id: str,
    gene_or_reaction_query: str,
    config: ExternalFetchConfig,
) -> tuple[ExternalReactionAssociation, ...]:
    """Read external reaction/GPR associations from an explicit offline cache.

    Round 6 intentionally does not perform live API calls. Online refreshers can
    populate the cache in later rounds; this function only filters local records.
    """

    if not config.offline_cache_dir:
        return ()
    records = _load_association_cache_records(Path(config.offline_cache_dir))
    source = source_database.strip().lower()
    model = _normalize_token(model_id)
    query_tokens = _match_tokens(gene_or_reaction_query)
    matched = [
        record
        for record in records
        if record.provenance.source_database.strip().lower() == source
        and _normalize_token(record.external_model_id) == model
        and query_tokens
        and query_tokens & _association_tokens(record)
    ]
    return tuple(sorted(matched, key=lambda record: record.cache_key))


def parse_sbml_gpr_associations(
    sbml_path: Path,
    *,
    source_database: str,
    source_model_id: str,
) -> tuple[ExternalReactionAssociation, ...]:
    """Parse reaction-level GPR associations from a local SBML file."""

    text = sbml_path.read_text(encoding="utf-8")
    root = ElementTree.fromstring(text)
    gene_labels = _gene_product_labels(root)
    retrieved_at = utc_now_iso()
    raw_hash = sha256_text(text)
    associations: list[ExternalReactionAssociation] = []
    for reaction in _iter_elements(root, "reaction"):
        reaction_id = _attr(reaction, "id")
        if not reaction_id:
            continue
        association_node = _first_descendant(reaction, "geneProductAssociation")
        if association_node is None:
            continue
        gene_rule = _rule_from_association(association_node, gene_labels)
        gene_ids = _gene_ids_from_association(association_node, gene_labels)
        if not gene_rule and not gene_ids:
            continue
        provenance = ExternalReferenceProvenance(
            source_database=source_database,
            source_version=source_model_id,
            source_url=str(sbml_path),
            source_query=source_model_id,
            retrieved_at=retrieved_at,
            raw_record_sha256=raw_hash,
        )
        associations.append(
            ExternalReactionAssociation(
                provenance=provenance,
                external_model_id=source_model_id,
                external_reaction_id=reaction_id,
                external_reaction_name=_attr(reaction, "name") or None,
                external_gene_ids=gene_ids,
                gene_rule=gene_rule,
                association_status="external_gpr_candidate",
            )
        )
    return tuple(sorted(associations, key=lambda record: record.cache_key))


def _load_association_cache_records(cache_path: Path) -> tuple[ExternalReactionAssociation, ...]:
    paths: tuple[Path, ...]
    if cache_path.is_file():
        paths = (cache_path,)
    else:
        paths = tuple(
            path
            for path in (
                cache_path / DEFAULT_RECORDS_FILENAME,
                cache_path / "external_reaction_associations.jsonl",
            )
            if path.exists()
        )
    records: list[ExternalReactionAssociation] = []
    for path in paths:
        records.extend(
            record
            for record in load_external_reference_cache(path)
            if isinstance(record, ExternalReactionAssociation)
        )
    return tuple(records)


def _gene_product_labels(root: ElementTree.Element) -> dict[str, str]:
    labels: dict[str, str] = {}
    for element in _iter_elements(root, "geneProduct"):
        gene_id = _attr(element, "id")
        if not gene_id:
            continue
        labels[gene_id] = _attr(element, "label") or _attr(element, "name") or gene_id
    return labels


def _rule_from_association(
    element: ElementTree.Element,
    gene_labels: dict[str, str],
) -> str | None:
    local_name = _local_name(element.tag)
    if local_name == "geneProductRef":
        gene_product = _attr(element, "geneProduct")
        return gene_labels.get(gene_product, gene_product) if gene_product else None
    child_rules = tuple(
        rule
        for child in list(element)
        for rule in (_rule_from_association(child, gene_labels),)
        if rule
    )
    if not child_rules:
        return None
    if local_name in {"and", "or"}:
        joined = f" {local_name} ".join(child_rules)
        return f"({joined})" if len(child_rules) > 1 else joined
    return child_rules[0] if len(child_rules) == 1 else f"({' and '.join(child_rules)})"


def _gene_ids_from_association(
    element: ElementTree.Element,
    gene_labels: dict[str, str],
) -> tuple[str, ...]:
    result: list[str] = []
    for child in _iter_elements(element, "geneProductRef"):
        gene_product = _attr(child, "geneProduct")
        gene_id = gene_labels.get(gene_product, gene_product)
        if gene_id and gene_id not in result:
            result.append(gene_id)
    return tuple(result)


def _association_tokens(record: ExternalReactionAssociation) -> set[str]:
    return _match_tokens(
        record.external_reaction_id,
        record.external_reaction_name,
        record.gene_rule,
        *record.external_gene_ids,
        *record.ec_numbers,
    )


def _match_tokens(*values: object) -> set[str]:
    return {
        token
        for value in values
        for token in (_normalize_token(value), *_split_tokens(value))
        if token
    }


def _split_tokens(value: object) -> tuple[str, ...]:
    text = str(value or "")
    for sep in ("/", ",", ";", "|", "(", ")"):
        text = text.replace(sep, " ")
    return tuple(_normalize_token(part) for part in text.split() if _normalize_token(part))


def _normalize_token(value: object) -> str:
    return str(value or "").strip().upper().replace(" ", "").replace("_", "").replace("-", "")


def _iter_elements(element: ElementTree.Element, local_name: str) -> Iterable[ElementTree.Element]:
    for child in element.iter():
        if _local_name(child.tag) == local_name:
            yield child


def _first_descendant(element: ElementTree.Element, local_name: str) -> ElementTree.Element | None:
    return next(_iter_elements(element, local_name), None)


def _attr(element: ElementTree.Element, local_name: str) -> str:
    for key, value in element.attrib.items():
        if _local_name(key) == local_name:
            return str(value)
    return ""


def _local_name(value: str) -> str:
    return value.rsplit("}", 1)[-1] if "}" in value else value


__all__ = [
    "fetch_external_model_reaction_associations",
    "parse_sbml_gpr_associations",
]
