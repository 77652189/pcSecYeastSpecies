from __future__ import annotations

import csv
import json
from json import JSONDecodeError
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
from collections.abc import Callable
from urllib.parse import urlencode
from urllib.request import Request, urlopen


DEFAULT_GENE_EVIDENCE_CACHE_DIR = Path("local_runs") / "gene_evidence_cache"
DEFAULT_GENE_EVIDENCE_CACHE = DEFAULT_GENE_EVIDENCE_CACHE_DIR / "gene_evidence.json"
DEFAULT_GENE_EVIDENCE_SUMMARY = DEFAULT_GENE_EVIDENCE_CACHE_DIR / "gene_evidence_summary.json"
DEFAULT_UNIPROT_PROTEOME = "UP000000314"
DEFAULT_KEGG_ORGANISM = "ppa"
WET_LAB_READY = "database_supported_experiment_candidate"
WET_LAB_REVIEW = "manual_review_required"
WET_LAB_MODEL_ONLY = "model_only_not_experiment_ready"

TIER_MODEL_EXECUTABLE = "model_executable"
TIER_EVIDENCE_SUPPORTED = "evidence_supported"
TIER_EXPERIMENT_CALIBRATED = "experiment_calibrated"
TIER_MANUAL_REVIEW = "manual_review_required"
TIER_GROWTH_RISK = "not_recommended_growth_risk"

_EXPERIMENT_CALIBRATED_SOURCES = {"internal_curated", "same_target_same_host_same_intervention"}


@dataclass(frozen=True)
class GeneExternalEvidence:
    gene_id: str
    canonical_gene_id: str
    model_gene_id: str = ""
    standard_gene_symbol: str = ""
    display_name: str = ""
    aliases: tuple[str, ...] = ()
    external_ids: dict[str, str] | None = None
    protein_name: str = ""
    function_annotation: str = ""
    subcellular_location: str = ""
    ec_numbers: tuple[str, ...] = ()
    go_terms: tuple[str, ...] = ()
    ortholog_symbol: str = ""
    wet_lab_readiness: str = WET_LAB_MODEL_ONLY
    evidence_sources: tuple[str, ...] = ()
    evidence_confidence: str = "unreviewed"
    last_refreshed: str = ""

    def to_dict(self) -> dict[str, Any]:
        payload = asdict(self)
        payload["aliases"] = list(self.aliases)
        payload["external_ids"] = dict(self.external_ids or {})
        payload["ec_numbers"] = list(self.ec_numbers)
        payload["go_terms"] = list(self.go_terms)
        payload["evidence_sources"] = list(self.evidence_sources)
        return payload


@dataclass(frozen=True)
class GenePhenotypeEvidence:
    intervention_type: str
    essentiality_status: str = ""
    secretion_screen_effect: str = ""
    evidence_source: str = ""
    evidence_confidence: str = "unreviewed"
    target_protein_context: tuple[str, ...] = ()
    recommended_use: str = ""

    def to_dict(self) -> dict[str, Any]:
        payload = asdict(self)
        payload["target_protein_context"] = list(self.target_protein_context)
        return payload


CURATED_PHENOTYPE_EVIDENCE: dict[str, tuple[GenePhenotypeEvidence, ...]] = {
    "ESSENTIAL_TEST": (
        GenePhenotypeEvidence(
            intervention_type="KO",
            essentiality_status="essential",
            secretion_screen_effect="growth_risk",
            evidence_source="internal_curated",
            evidence_confidence="high",
            target_protein_context=("OPN_ALPHA_FULL_PROJECT",),
            recommended_use="do_not_ko_growth_risk",
        ),
    ),
    "POSITIVE_SCREEN_TEST": (
        GenePhenotypeEvidence(
            intervention_type="OE",
            essentiality_status="nonessential",
            secretion_screen_effect="screen_positive",
            evidence_source="public_screen_positive",
            evidence_confidence="medium",
            target_protein_context=("secreted_reporter",),
            recommended_use="evidence_supported_candidate",
        ),
    ),
    "CALIBRATED_TEST": (
        GenePhenotypeEvidence(
            intervention_type="OE",
            essentiality_status="nonessential",
            secretion_screen_effect="screen_positive",
            evidence_source="same_target_same_host_same_intervention",
            evidence_confidence="high",
            target_protein_context=("OPN_ALPHA_FULL_PROJECT",),
            recommended_use="experiment_calibrated_candidate",
        ),
    ),
    "PDI1": (
        GenePhenotypeEvidence(
            intervention_type="OE",
            essentiality_status="nonessential",
            secretion_screen_effect="screen_positive",
            evidence_source="public_screen_positive",
            evidence_confidence="medium",
            target_protein_context=("dsb_rich_secreted_protein",),
            recommended_use="evidence_supported_oe_proxy",
        ),
    ),
}


def load_curated_phenotype_evidence() -> dict[str, tuple[GenePhenotypeEvidence, ...]]:
    return CURATED_PHENOTYPE_EVIDENCE


def phenotype_evidence_for_candidate(
    gene_id: str,
    intervention_type: str,
    target_protein_context: str | None = None,
    evidence_by_gene: dict[str, tuple[GenePhenotypeEvidence, ...]] | None = None,
    aliases: tuple[str, ...] = (),
) -> GenePhenotypeEvidence | None:
    query_keys = _candidate_gene_keys(gene_id, aliases)
    intervention_key = _phenotype_intervention_key(intervention_type)
    if not intervention_key:
        return None
    records_by_gene = load_curated_phenotype_evidence() if evidence_by_gene is None else evidence_by_gene
    candidates: list[GenePhenotypeEvidence] = []
    for key in query_keys:
        candidates.extend(records_by_gene.get(key, ()))
        candidates.extend(records_by_gene.get(key.upper(), ()))
    matched = [
        record
        for record in candidates
        if _phenotype_intervention_key(record.intervention_type) == intervention_key
    ]
    if not matched:
        return None
    if intervention_key == "KO":
        essential = [record for record in matched if record.essentiality_status == "essential"]
        if essential:
            return _best_phenotype_record(essential)
    context = str(target_protein_context or "").strip()
    context_matches = [record for record in matched if _target_context_matches(record, context)]
    if context_matches:
        return _best_phenotype_record(context_matches)
    return _best_phenotype_record(matched)


def recommendation_tier_for_candidate(
    *,
    gene_id: str = "",
    intervention_type: str = "",
    target_protein_context: str | None = None,
    model_gpr_executable: bool = False,
    oe_reaction_proxy: bool = False,
    resolved: bool = True,
    database_annotation_available: bool = False,
    phenotype_evidence: GenePhenotypeEvidence | None = None,
    aliases: tuple[str, ...] = (),
) -> tuple[str, str, GenePhenotypeEvidence | None]:
    intervention_key = _phenotype_intervention_key(intervention_type)
    if phenotype_evidence is not None:
        evidence = (
            phenotype_evidence
            if _phenotype_intervention_key(phenotype_evidence.intervention_type) == intervention_key
            else None
        )
    else:
        evidence = phenotype_evidence_for_candidate(
            gene_id,
            intervention_type,
            target_protein_context=target_protein_context,
            aliases=aliases,
        )
    if evidence and intervention_key == "KO" and evidence.essentiality_status == "essential":
        return (
            TIER_GROWTH_RISK,
            "KO candidate has curated essentiality evidence and is treated as a growth-risk target.",
            evidence,
        )
    if not resolved:
        return (TIER_MANUAL_REVIEW, "Candidate is unresolved in the current model.", evidence)
    if evidence:
        context_matches = _target_context_matches(evidence, str(target_protein_context or "").strip())
        if (
            evidence.evidence_source in _EXPERIMENT_CALIBRATED_SOURCES
            and _is_high_confidence(evidence.evidence_confidence)
            and context_matches
        ):
            return (
                TIER_EXPERIMENT_CALIBRATED,
                "High-confidence phenotype evidence matches intervention and target context.",
                evidence,
            )
        return (
            TIER_EVIDENCE_SUPPORTED,
            "Phenotype evidence supports this intervention, but source or target context is not fully calibrated.",
            evidence,
        )
    if model_gpr_executable or oe_reaction_proxy:
        return (
            TIER_MODEL_EXECUTABLE,
            "Current model can execute the GPR KO or reaction-level OE proxy, but no phenotype evidence matched.",
            None,
        )
    if database_annotation_available:
        return (
            TIER_MANUAL_REVIEW,
            "Only database annotation is available; no intervention-specific phenotype evidence matched.",
            None,
        )
    return (TIER_MANUAL_REVIEW, "No executable model support or phenotype evidence matched.", None)


def load_gene_evidence_cache(path: Path | str | None = None) -> dict[str, GeneExternalEvidence]:
    cache_path = Path(path) if path is not None else DEFAULT_GENE_EVIDENCE_CACHE
    if not cache_path.exists():
        return {}
    try:
        payload = json.loads(cache_path.read_text(encoding="utf-8"))
    except (OSError, JSONDecodeError, UnicodeDecodeError):
        return {}
    records = payload.get("genes", payload) if isinstance(payload, dict) else payload
    evidence: dict[str, GeneExternalEvidence] = {}
    if not isinstance(records, list):
        return evidence
    for item in records:
        if not isinstance(item, dict):
            continue
        record = gene_external_evidence_from_mapping(item)
        if record.canonical_gene_id:
            evidence[record.canonical_gene_id] = record
    return evidence


def gene_external_evidence_from_mapping(item: dict[str, Any]) -> GeneExternalEvidence:
    gene_id = str(item.get("gene_id") or item.get("canonical_gene_id") or "").strip()
    canonical_gene_id = str(item.get("canonical_gene_id") or gene_id).strip()
    aliases = tuple(str(value).strip() for value in item.get("aliases") or () if str(value).strip())
    sources = tuple(str(value).strip() for value in item.get("evidence_sources") or () if str(value).strip())
    ec_numbers = tuple(str(value).strip() for value in item.get("ec_numbers") or () if str(value).strip())
    go_terms = tuple(str(value).strip() for value in item.get("go_terms") or () if str(value).strip())
    external_ids = item.get("external_ids") if isinstance(item.get("external_ids"), dict) else {}
    return GeneExternalEvidence(
        gene_id=gene_id,
        canonical_gene_id=canonical_gene_id,
        model_gene_id=str(item.get("model_gene_id") or canonical_gene_id or gene_id),
        standard_gene_symbol=str(item.get("standard_gene_symbol") or ""),
        display_name=str(item.get("display_name") or item.get("protein_name") or gene_id),
        aliases=aliases,
        external_ids={str(key): str(value) for key, value in external_ids.items()},
        protein_name=str(item.get("protein_name") or ""),
        function_annotation=str(item.get("function_annotation") or ""),
        subcellular_location=str(item.get("subcellular_location") or ""),
        ec_numbers=ec_numbers,
        go_terms=go_terms,
        ortholog_symbol=str(item.get("ortholog_symbol") or ""),
        wet_lab_readiness=str(item.get("wet_lab_readiness") or WET_LAB_MODEL_ONLY),
        evidence_sources=sources,
        evidence_confidence=str(item.get("evidence_confidence") or "unreviewed"),
        last_refreshed=str(item.get("last_refreshed") or ""),
    )


def build_gene_evidence_cache(
    gene_ids: list[str] | tuple[str, ...],
    output_path: Path | str = DEFAULT_GENE_EVIDENCE_CACHE,
    summary_path: Path | str = DEFAULT_GENE_EVIDENCE_SUMMARY,
    *,
    uniprot_proteome: str = DEFAULT_UNIPROT_PROTEOME,
    kegg_organism: str = DEFAULT_KEGG_ORGANISM,
    timeout_seconds: int = 10,
    progress: Callable[[str], None] | None = None,
) -> dict[str, Any]:
    """Build a wet-lab annotation cache for model genes from authoritative online sources."""
    unique_gene_ids = tuple(dict.fromkeys(str(gene_id).strip() for gene_id in gene_ids if str(gene_id).strip()))
    uniprot_by_gene = _fetch_uniprot_proteome_by_gene(
        uniprot_proteome,
        gene_ids=unique_gene_ids,
        timeout_seconds=timeout_seconds,
        progress=progress,
    )
    if progress:
        progress(f"Fetching KEGG `{kegg_organism}` gene list...")
    kegg_by_gene = _fetch_kegg_gene_descriptions(kegg_organism, timeout_seconds=timeout_seconds)
    refreshed = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    records = [
        _evidence_record_for_gene(
            gene_id,
            uniprot_by_gene.get(gene_id),
            kegg_by_gene.get(gene_id),
            refreshed,
            kegg_organism=kegg_organism,
        )
        for gene_id in unique_gene_ids
    ]
    summary = summarize_gene_evidence_records(records)
    payload = {
        "schema_version": 1,
        "source_policy": "uniprot_proteome_plus_kegg_locus_tag_cache",
        "uniprot_proteome": uniprot_proteome,
        "kegg_organism": kegg_organism,
        "last_refreshed": refreshed,
        "summary": summary,
        "genes": records,
    }
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    summary_output = Path(summary_path)
    summary_output.parent.mkdir(parents=True, exist_ok=True)
    summary_output.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    return summary


def summarize_gene_evidence_records(records: list[dict[str, Any]]) -> dict[str, Any]:
    total = len(records)
    readiness_counts: dict[str, int] = {}
    source_counts: dict[str, int] = {}
    prefix_counts: dict[str, int] = {}
    for record in records:
        readiness = str(record.get("wet_lab_readiness") or WET_LAB_MODEL_ONLY)
        readiness_counts[readiness] = readiness_counts.get(readiness, 0) + 1
        prefix = str(record.get("model_gene_id") or record.get("gene_id") or "").split("_", 1)[0]
        prefix_counts[prefix] = prefix_counts.get(prefix, 0) + 1
        for source in record.get("evidence_sources") or ():
            source_counts[str(source)] = source_counts.get(str(source), 0) + 1
    ready = readiness_counts.get(WET_LAB_READY, 0)
    review = readiness_counts.get(WET_LAB_REVIEW, 0)
    return {
        "total_genes": total,
        "database_supported_count": ready,
        "manual_review_count": review,
        "model_only_count": readiness_counts.get(WET_LAB_MODEL_ONLY, 0),
        "database_supported_fraction": round(ready / total, 4) if total else 0.0,
        "readiness_counts": readiness_counts,
        "source_counts": source_counts,
        "prefix_counts": prefix_counts,
    }


def _fetch_uniprot_proteome_by_gene(
    proteome_id: str,
    *,
    gene_ids: tuple[str, ...] | None = None,
    timeout_seconds: int,
    progress: Callable[[str], None] | None = None,
) -> dict[str, dict[str, str]]:
    if gene_ids:
        return _fetch_uniprot_records_for_genes(
            gene_ids,
            proteome_id=proteome_id,
            timeout_seconds=timeout_seconds,
            progress=progress,
        )
    fields = (
        "accession,gene_names,gene_primary,protein_name,organism_name,"
        "xref_geneid,xref_kegg,xref_refseq,go_p,ec"
    )
    params = urlencode(
        {
            "query": f"proteome:{proteome_id}",
            "format": "tsv",
            "fields": fields,
            "size": "500",
        }
    )
    url = f"https://rest.uniprot.org/uniprotkb/search?{params}"
    rows: list[dict[str, str]] = []
    while url:
        response_text, next_url = _read_text_url_with_next_link(url, timeout_seconds=timeout_seconds)
        rows.extend(_tsv_dict_rows(response_text))
        url = next_url
    by_gene: dict[str, dict[str, str]] = {}
    for row in rows:
        gene_names = _split_external_list(row.get("Gene Names", ""))
        for gene_name in gene_names:
            by_gene.setdefault(gene_name, row)
    return by_gene


def _fetch_uniprot_records_for_genes(
    gene_ids: tuple[str, ...],
    *,
    proteome_id: str,
    timeout_seconds: int,
    chunk_size: int = 100,
    progress: Callable[[str], None] | None = None,
) -> dict[str, dict[str, str]]:
    fields = (
        "accession,gene_names,gene_primary,protein_name,organism_name,"
        "xref_geneid,xref_kegg,xref_refseq,go_p,ec"
    )
    by_gene: dict[str, dict[str, str]] = {}
    chunks = _chunks(gene_ids, chunk_size)
    for index, chunk in enumerate(chunks, start=1):
        gene_query = " OR ".join(f"gene_exact:{gene_id}" for gene_id in chunk if gene_id.startswith("PAS_"))
        if not gene_query:
            continue
        if progress:
            progress(f"Fetching UniProt chunk {index}/{len(chunks)} ({len(chunk)} genes)...")
        params = urlencode(
            {
                "query": f"({gene_query}) AND proteome:{proteome_id}",
                "format": "tsv",
                "fields": fields,
                "size": str(max(50, len(chunk) * 2)),
            }
        )
        url = f"https://rest.uniprot.org/uniprotkb/search?{params}"
        try:
            response_text, _ = _read_text_url_with_next_link(url, timeout_seconds=timeout_seconds)
        except OSError:
            continue
        for row in _tsv_dict_rows(response_text):
            for gene_name in _split_external_list(row.get("Gene Names", "")):
                by_gene.setdefault(gene_name, row)
    return by_gene


def _chunks(values: tuple[str, ...], size: int) -> list[tuple[str, ...]]:
    return [values[index:index + size] for index in range(0, len(values), size)]


def _fetch_kegg_gene_descriptions(
    organism: str,
    *,
    timeout_seconds: int,
) -> dict[str, str]:
    url = f"https://rest.kegg.jp/list/{organism}"
    try:
        text, _ = _read_text_url_with_next_link(url, timeout_seconds=timeout_seconds)
    except OSError:
        return {}
    descriptions: dict[str, str] = {}
    prefix = f"{organism}:"
    for line in text.splitlines():
        if not line.strip() or "\t" not in line:
            continue
        key, description = line.split("\t", 1)
        gene_id = key[len(prefix):] if key.startswith(prefix) else key
        descriptions[gene_id] = description.strip()
    return descriptions


def _read_text_url_with_next_link(url: str, *, timeout_seconds: int) -> tuple[str, str | None]:
    request = Request(url, headers={"User-Agent": "pcSecPichia-gene-evidence-cache/1.0"})
    with urlopen(request, timeout=timeout_seconds) as response:
        text = response.read().decode("utf-8")
        link_header = response.headers.get("Link", "")
    return text, _next_link_from_header(link_header)


def _next_link_from_header(header: str) -> str | None:
    for item in header.split(","):
        item = item.strip()
        if 'rel="next"' not in item:
            continue
        start = item.find("<")
        end = item.find(">")
        if start >= 0 and end > start:
            return item[start + 1:end]
    return None


def _tsv_dict_rows(text: str) -> list[dict[str, str]]:
    if not text.strip():
        return []
    return [dict(row) for row in csv.DictReader(text.splitlines(), delimiter="\t")]


def _evidence_record_for_gene(
    gene_id: str,
    uniprot_row: dict[str, str] | None,
    kegg_description: str | None,
    refreshed: str,
    *,
    kegg_organism: str,
) -> dict[str, Any]:
    if uniprot_row:
        external_ids = _external_ids_from_uniprot(uniprot_row)
        protein_name = _clean_uniprot_protein_name(uniprot_row.get("Protein names", ""))
        function_annotation = protein_name or kegg_description or ""
        standard_gene_symbol = _standard_symbol_from_uniprot(uniprot_row, gene_id)
        display_name = standard_gene_symbol or protein_name or gene_id
        return {
            "gene_id": gene_id,
            "canonical_gene_id": gene_id,
            "model_gene_id": gene_id,
            "standard_gene_symbol": standard_gene_symbol,
            "display_name": display_name,
            "aliases": _aliases_from_uniprot(uniprot_row, gene_id),
            "external_ids": external_ids,
            "protein_name": protein_name,
            "function_annotation": function_annotation,
            "subcellular_location": "",
            "ec_numbers": _split_external_list(uniprot_row.get("EC number", "")),
            "go_terms": _split_semicolon_values(uniprot_row.get("Gene Ontology (biological process)", "")),
            "ortholog_symbol": "",
            "wet_lab_readiness": WET_LAB_READY,
            "evidence_sources": _evidence_sources(external_ids, has_uniprot=True, has_kegg=bool(kegg_description)),
            "evidence_confidence": "high_exact_locus_tag",
            "last_refreshed": refreshed,
        }
    if kegg_description:
        return {
            "gene_id": gene_id,
            "canonical_gene_id": gene_id,
            "model_gene_id": gene_id,
            "standard_gene_symbol": "",
            "display_name": kegg_description,
            "aliases": [],
            "external_ids": {"kegg": f"{kegg_organism}:{gene_id}"},
            "protein_name": kegg_description,
            "function_annotation": kegg_description,
            "subcellular_location": "",
            "ec_numbers": [],
            "go_terms": [],
            "ortholog_symbol": "",
            "wet_lab_readiness": WET_LAB_REVIEW,
            "evidence_sources": ["KEGG"],
            "evidence_confidence": "medium_exact_kegg_locus_tag",
            "last_refreshed": refreshed,
        }
    return {
        "gene_id": gene_id,
        "canonical_gene_id": gene_id,
        "model_gene_id": gene_id,
        "standard_gene_symbol": "",
        "display_name": gene_id,
        "aliases": [],
        "external_ids": {},
        "protein_name": "",
        "function_annotation": "",
        "subcellular_location": "",
        "ec_numbers": [],
        "go_terms": [],
        "ortholog_symbol": "",
        "wet_lab_readiness": WET_LAB_MODEL_ONLY,
        "evidence_sources": [],
        "evidence_confidence": "low_model_only",
        "last_refreshed": refreshed,
    }


def _external_ids_from_uniprot(row: dict[str, str]) -> dict[str, str]:
    external_ids = {"uniprot": row.get("Entry", "").strip()}
    if row.get("GeneID"):
        external_ids["ncbi_gene"] = _strip_terminal_semicolon(row["GeneID"])
    if row.get("KEGG"):
        external_ids["kegg"] = _strip_terminal_semicolon(row["KEGG"])
    if row.get("RefSeq"):
        external_ids["refseq"] = _strip_terminal_semicolon(row["RefSeq"])
    return {key: value for key, value in external_ids.items() if value}


def _standard_symbol_from_uniprot(row: dict[str, str], model_gene_id: str) -> str:
    primary = row.get("Gene Names (primary)", "").strip()
    if primary and primary != model_gene_id:
        return primary
    for gene_name in _split_external_list(row.get("Gene Names", "")):
        if gene_name != model_gene_id and not gene_name.startswith(("PAS_", "AT250_")):
            return gene_name
    return ""


def _aliases_from_uniprot(row: dict[str, str], model_gene_id: str) -> list[str]:
    aliases: list[str] = []
    for gene_name in _split_external_list(row.get("Gene Names", "")):
        if gene_name != model_gene_id and gene_name not in aliases:
            aliases.append(gene_name)
    return aliases


def _evidence_sources(external_ids: dict[str, str], *, has_uniprot: bool, has_kegg: bool) -> list[str]:
    sources: list[str] = []
    if has_uniprot:
        sources.append("UniProt")
    if external_ids.get("ncbi_gene"):
        sources.append("NCBI Gene")
    if external_ids.get("refseq"):
        sources.append("RefSeq")
    if has_kegg or external_ids.get("kegg"):
        sources.append("KEGG")
    return sources


def _clean_uniprot_protein_name(value: str) -> str:
    return str(value or "").split(" (", 1)[0].strip()


def _split_external_list(value: str) -> list[str]:
    return [item.strip() for item in str(value or "").replace(";", " ").split() if item.strip()]


def _split_semicolon_values(value: str) -> list[str]:
    return [item.strip() for item in str(value or "").split(";") if item.strip()]


def _strip_terminal_semicolon(value: str) -> str:
    return str(value or "").strip().rstrip(";")


def _candidate_gene_keys(gene_id: str, aliases: tuple[str, ...] = ()) -> tuple[str, ...]:
    keys: list[str] = []
    for value in (gene_id, *aliases):
        text = str(value or "").strip()
        if text and text not in keys:
            keys.append(text)
    return tuple(keys)


def _phenotype_intervention_key(intervention_type: str) -> str:
    text = str(intervention_type or "").strip().upper()
    if text.startswith("KO"):
        return "KO"
    if text.startswith("OE") or "OVEREXPRESSION" in text:
        return "OE"
    return ""


def _target_context_matches(record: GenePhenotypeEvidence, target_protein_context: str) -> bool:
    contexts = tuple(str(item).strip().lower() for item in record.target_protein_context if str(item).strip())
    if not contexts:
        return True
    current = str(target_protein_context or "").strip().lower()
    if not current:
        return False
    return current in contexts


def _is_high_confidence(value: str) -> bool:
    text = str(value or "").strip().lower()
    return text == "high" or text.startswith("high_") or text.endswith("_high")


def _best_phenotype_record(records: list[GenePhenotypeEvidence]) -> GenePhenotypeEvidence:
    def sort_key(record: GenePhenotypeEvidence) -> tuple[int, int]:
        source_rank = 0 if record.evidence_source in _EXPERIMENT_CALIBRATED_SOURCES else 1
        confidence_rank = 0 if _is_high_confidence(record.evidence_confidence) else 1
        return (source_rank, confidence_rank)

    return sorted(records, key=sort_key)[0]


def evidence_for_gene(
    gene_id: str,
    evidence_by_gene: dict[str, GeneExternalEvidence],
) -> GeneExternalEvidence | None:
    query = str(gene_id or "").strip()
    if not query:
        return None
    if query in evidence_by_gene:
        return evidence_by_gene[query]
    query_lower = query.lower()
    for record in evidence_by_gene.values():
        if query_lower == record.gene_id.lower():
            return record
        if any(query_lower == alias.lower() for alias in record.aliases):
            return record
    return None


def resolve_gene_identifier(
    gene_id: str,
    evidence_by_gene: dict[str, GeneExternalEvidence],
    model: Any | None = None,
) -> str:
    query = str(gene_id or "").strip()
    if not query:
        return ""
    gene_index = getattr(model, "gene_index", {}) if model is not None else {}
    if query in gene_index:
        return query
    record = evidence_for_gene(query, evidence_by_gene)
    if record is None:
        return query
    canonical = record.canonical_gene_id or record.gene_id or query
    if model is not None and canonical not in gene_index:
        return query
    return canonical


__all__ = [
    "DEFAULT_GENE_EVIDENCE_CACHE",
    "DEFAULT_GENE_EVIDENCE_CACHE_DIR",
    "DEFAULT_GENE_EVIDENCE_SUMMARY",
    "WET_LAB_MODEL_ONLY",
    "WET_LAB_READY",
    "WET_LAB_REVIEW",
    "TIER_EVIDENCE_SUPPORTED",
    "TIER_EXPERIMENT_CALIBRATED",
    "TIER_GROWTH_RISK",
    "TIER_MANUAL_REVIEW",
    "TIER_MODEL_EXECUTABLE",
    "CURATED_PHENOTYPE_EVIDENCE",
    "GeneExternalEvidence",
    "GenePhenotypeEvidence",
    "build_gene_evidence_cache",
    "evidence_for_gene",
    "gene_external_evidence_from_mapping",
    "load_gene_evidence_cache",
    "load_curated_phenotype_evidence",
    "phenotype_evidence_for_candidate",
    "recommendation_tier_for_candidate",
    "resolve_gene_identifier",
    "summarize_gene_evidence_records",
]
