from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any


GENE_ID_STANDARDIZATION_SCHEMA_VERSION = 1


@dataclass(frozen=True)
class PichiaGeneIdStandardName:
    gene_id: str
    display_name: str
    standard_symbol: str = ""
    protein_name: str = ""
    external_ids: dict[str, str] = field(default_factory=dict)
    annotation_sources: tuple[str, ...] = ()
    annotation_confidence: str = ""
    model_operable: bool = True
    gpr_status: str = ""

    def to_dict(self) -> dict[str, Any]:
        payload = asdict(self)
        payload["annotation_sources"] = list(self.annotation_sources)
        payload["external_ids"] = dict(self.external_ids)
        return payload


def build_pichia_gene_id_standardization_rows(
    full_model_rows: list[dict[str, object]] | tuple[dict[str, object], ...],
) -> tuple[PichiaGeneIdStandardName, ...]:
    rows: list[PichiaGeneIdStandardName] = []
    seen: set[str] = set()
    for source in full_model_rows:
        gene_id = str(source.get("gene_id") or "").strip()
        if not gene_id:
            raise ValueError("Every standard naming row requires a non-empty gene_id.")
        if gene_id in seen:
            raise ValueError(f"Found duplicate gene_id in standard naming rows: {gene_id}")
        seen.add(gene_id)
        external_ids = _external_ids(source.get("external_ids"))
        annotation_sources = _annotation_sources(source, external_ids)
        standard_symbol = str(source.get("standard_gene_symbol") or source.get("standard_symbol") or "").strip()
        protein_name = str(source.get("protein_name") or "").strip()
        display_name = str(source.get("display_name") or standard_symbol or protein_name or gene_id).strip()
        rows.append(
            PichiaGeneIdStandardName(
                gene_id=gene_id,
                display_name=display_name,
                standard_symbol=standard_symbol,
                protein_name=protein_name,
                external_ids=external_ids,
                annotation_sources=annotation_sources,
                annotation_confidence=_annotation_confidence(source, annotation_sources),
                model_operable=True,
                gpr_status=_gpr_status(source),
            )
        )
    return tuple(rows)


def summarize_pichia_gene_id_standardization_rows(
    rows: tuple[PichiaGeneIdStandardName, ...] | list[PichiaGeneIdStandardName],
) -> dict[str, Any]:
    total = len(rows)
    model_only = [row.gene_id for row in rows if _is_model_only(row)]
    annotated = total - len(model_only)
    source_counts: dict[str, int] = {}
    gpr_status_counts: dict[str, int] = {}
    for row in rows:
        gpr_status_counts[row.gpr_status] = gpr_status_counts.get(row.gpr_status, 0) + 1
        for source in row.annotation_sources:
            source_counts[source] = source_counts.get(source, 0) + 1
    return {
        "schema_version": GENE_ID_STANDARDIZATION_SCHEMA_VERSION,
        "total_genes": total,
        "annotated_gene_count": annotated,
        "model_only_count": len(model_only),
        "model_only_gene_ids": model_only,
        "annotation_source_counts": dict(sorted(source_counts.items())),
        "gpr_status_counts": dict(sorted(gpr_status_counts.items())),
    }


def write_pichia_gene_id_standardization_cache(
    rows: tuple[PichiaGeneIdStandardName, ...] | list[PichiaGeneIdStandardName],
    output_path: Path,
) -> None:
    payload = {
        "schema_version": GENE_ID_STANDARDIZATION_SCHEMA_VERSION,
        "summary": summarize_pichia_gene_id_standardization_rows(tuple(rows)),
        "rows": [row.to_dict() for row in rows],
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")


def load_pichia_gene_id_standardization_cache(path: Path) -> tuple[PichiaGeneIdStandardName, ...]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict) or payload.get("schema_version") != GENE_ID_STANDARDIZATION_SCHEMA_VERSION:
        return ()
    rows = payload.get("rows")
    if not isinstance(rows, list):
        return ()
    return tuple(_row_from_payload(row) for row in rows if isinstance(row, dict))


def _row_from_payload(payload: dict[str, object]) -> PichiaGeneIdStandardName:
    return PichiaGeneIdStandardName(
        gene_id=str(payload.get("gene_id") or ""),
        display_name=str(payload.get("display_name") or ""),
        standard_symbol=str(payload.get("standard_symbol") or ""),
        protein_name=str(payload.get("protein_name") or ""),
        external_ids=_external_ids(payload.get("external_ids")),
        annotation_sources=_tuple_strings(payload.get("annotation_sources")),
        annotation_confidence=str(payload.get("annotation_confidence") or ""),
        model_operable=bool(payload.get("model_operable", True)),
        gpr_status=str(payload.get("gpr_status") or ""),
    )


def _annotation_sources(source: dict[str, object], external_ids: dict[str, str]) -> tuple[str, ...]:
    sources = list(_tuple_strings(source.get("evidence_sources")))
    if external_ids.get("uniprot") and "UniProt" not in sources:
        sources.append("UniProt")
    if external_ids.get("ncbi_gene") and "NCBI Gene" not in sources:
        sources.append("NCBI Gene")
    if external_ids.get("kegg") and "KEGG" not in sources:
        sources.append("KEGG")
    if external_ids.get("refseq") and "RefSeq" not in sources:
        sources.append("RefSeq")
    if not sources:
        return ("model_only",)
    return tuple(dict.fromkeys(sources))


def _annotation_confidence(source: dict[str, object], annotation_sources: tuple[str, ...]) -> str:
    value = str(source.get("evidence_confidence") or "").strip()
    if value:
        return value
    return "low_model_only" if annotation_sources == ("model_only",) else "unreviewed"


def _gpr_status(source: dict[str, object]) -> str:
    ko_status = str(source.get("ko_support_status") or "")
    oe_status = str(source.get("oe_support_status") or "")
    has_affected_reactions = bool(source.get("affected_reactions") or source.get("gpr_rules"))
    ko_executable = ko_status == "ko_runnable_gpr_gene_deletion"
    oe_executable = oe_status == "oe_runnable_reaction_proxy"
    if ko_executable and oe_executable:
        return "ko_and_oe_model_executable"
    if ko_executable:
        return "ko_model_executable"
    if oe_executable:
        return "oe_proxy_executable"
    if has_affected_reactions:
        return "gpr_linked_review_required"
    return "model_gene_no_gpr_effect"


def _is_model_only(row: PichiaGeneIdStandardName) -> bool:
    return row.annotation_sources == ("model_only",)


def _external_ids(value: object) -> dict[str, str]:
    if not isinstance(value, dict):
        return {}
    return {str(key): str(item) for key, item in value.items() if str(key).strip() and str(item).strip()}


def _tuple_strings(value: object) -> tuple[str, ...]:
    if value is None:
        return ()
    if isinstance(value, str):
        return (value,) if value.strip() else ()
    if isinstance(value, (list, tuple, set)):
        return tuple(str(item).strip() for item in value if str(item).strip())
    return ()


__all__ = [
    "GENE_ID_STANDARDIZATION_SCHEMA_VERSION",
    "PichiaGeneIdStandardName",
    "build_pichia_gene_id_standardization_rows",
    "load_pichia_gene_id_standardization_cache",
    "summarize_pichia_gene_id_standardization_rows",
    "write_pichia_gene_id_standardization_cache",
]
