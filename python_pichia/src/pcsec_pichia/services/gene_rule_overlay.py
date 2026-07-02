from __future__ import annotations

import copy
import json
from dataclasses import asdict, dataclass, replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
from urllib.parse import urlencode
from urllib.request import Request, urlopen


DEFAULT_GENE_RULE_EVIDENCE_CACHE = (
    Path("local_runs") / "gene_rule_evidence_cache" / "gene_rule_evidence.json"
)
DEFAULT_GENE_RULE_EVIDENCE_SUMMARY = (
    Path("local_runs") / "gene_rule_evidence_cache" / "gene_rule_evidence_summary.json"
)
DEFAULT_GENE_RULE_EVIDENCE_REPORT = (
    Path("local_runs") / "gene_rule_evidence_cache" / "GENE_RULE_EVIDENCE_REPORT.md"
)
DEFAULT_UNIPROT_PROTEOME = "UP000000314"
DEFAULT_KEGG_ORGANISM = "ppa"

HIGH_CONFIDENCE = "high_exact_multi_source"
MEDIUM_CONFIDENCE = "medium_exact_single_source"
LOW_CONFIDENCE = "low_homology_or_name_only"
UNRESOLVED_CONFIDENCE = "unresolved"
PDI_COMPLEX_REACTION = "sec_PDI1_ERV2_Ero1p_complex_formation"
PDI_COMPLEX_NAMES = ("PDI1", "ERO1", "ERV2")
DEFAULT_TARGET_NAMES = ("PDI1", "ERO1", "ERV2", "KAR2", "OCH1", "PEP4", "PRB1")


@dataclass(frozen=True)
class GeneRuleEvidence:
    common_name: str
    candidate_locus_tag: str = ""
    external_ids: dict[str, str] | None = None
    protein_name: str = ""
    evidence_sources: tuple[str, ...] = ()
    confidence: str = UNRESOLVED_CONFIDENCE
    mapped_model_gene_present: bool = False
    target_reaction_ids: tuple[str, ...] = ()
    proposed_gr_rule: str = ""
    proposed_rule: str = ""
    rule_status: str = "not_executable"
    recommended_action: str = "manual_review_required"
    notes: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, object]:
        payload = asdict(self)
        payload["external_ids"] = dict(self.external_ids or {})
        payload["evidence_sources"] = list(self.evidence_sources)
        payload["target_reaction_ids"] = list(self.target_reaction_ids)
        payload["notes"] = list(self.notes)
        return payload


@dataclass(frozen=True)
class GprOverlayEntry:
    reaction_id: str
    gene_locus_tags: tuple[str, ...]
    common_names: tuple[str, ...]
    proposed_gr_rule: str
    proposed_rule: str
    rule_status: str
    recommended_action: str

    def to_dict(self) -> dict[str, object]:
        payload = asdict(self)
        payload["gene_locus_tags"] = list(self.gene_locus_tags)
        payload["common_names"] = list(self.common_names)
        return payload


@dataclass(frozen=True)
class GprOverlayResult:
    entries: tuple[GprOverlayEntry, ...]
    skipped: tuple[dict[str, object], ...]
    warnings: tuple[str, ...]
    result_status: str = "draft_external_evidence_overlay"

    def to_dict(self) -> dict[str, object]:
        return {
            "result_status": self.result_status,
            "entries": [entry.to_dict() for entry in self.entries],
            "skipped": [dict(item) for item in self.skipped],
            "warnings": list(self.warnings),
            "entry_count": len(self.entries),
        }


def build_gene_rule_evidence_cache(
    target_names: tuple[str, ...] = DEFAULT_TARGET_NAMES,
    *,
    output_path: Path | str | None = None,
    summary_path: Path | str | None = None,
    report_path: Path | str | None = None,
    model: Any | None = None,
    enable_online: bool = True,
    timeout_seconds: int = 10,
    progress: Any | None = None,
) -> dict[str, object]:
    """Build a conservative local evidence cache for secretion-engineering GPR overlay.

    The seed evidence is intentionally conservative. It records candidate locus tags
    and source hints, but executable overlay remains limited to high-confidence
    records with target reactions present in the current model.
    """
    requested = _dedupe_names(target_names)
    records = [_seed_evidence_for_name(name) for name in requested]
    source_rows = _fetch_source_rows_for_records(
        records,
        enable_online=enable_online,
        timeout_seconds=timeout_seconds,
        progress=progress,
    )
    records = [_enrich_record_with_source_evidence(record, source_rows) for record in records]
    overlay = build_gpr_overlay(model, {record.common_name: record for record in records}) if model is not None else None
    refreshed = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    payload = {
        "schema_version": 1,
        "species": "Komagataella phaffii GS115 / Pichia pastoris GS115",
        "last_refreshed": refreshed,
        "source_references": {
            "uniprot_proteome": f"https://www.uniprot.org/proteomes/{DEFAULT_UNIPROT_PROTEOME}",
            "uniprot_rest": "https://rest.uniprot.org/uniprotkb/search",
            "kegg_organism": f"https://www.kegg.jp/kegg-bin/show_organism?org={DEFAULT_KEGG_ORGANISM}",
            "kegg_rest_list": f"https://rest.kegg.jp/list/{DEFAULT_KEGG_ORGANISM}",
        },
        "records": [record.to_dict() for record in records],
        "source_policy": (
            "Executable overlays require high_exact_multi_source evidence and an "
            "existing target reaction. Lower-confidence records are display-only."
        ),
        "overlay_probe": overlay.to_dict() if overlay is not None else {},
    }
    output = Path(output_path) if output_path is not None else DEFAULT_GENE_RULE_EVIDENCE_CACHE
    summary = Path(summary_path) if summary_path is not None else DEFAULT_GENE_RULE_EVIDENCE_SUMMARY
    report = Path(report_path) if report_path is not None else DEFAULT_GENE_RULE_EVIDENCE_REPORT
    _write_json(output, payload)
    summary_payload = summarize_gene_rule_evidence_records(records, overlay=overlay)
    summary_payload.update({"cache_path": str(output), "report_path": str(report), "last_refreshed": refreshed})
    _write_json(summary, summary_payload)
    _write_report(report, _build_gene_rule_evidence_report(records, summary_payload, overlay))
    return summary_payload


def summarize_gene_rule_evidence_records(
    records: list[GeneRuleEvidence] | tuple[GeneRuleEvidence, ...],
    *,
    overlay: GprOverlayResult | None = None,
) -> dict[str, object]:
    confidence_counts: dict[str, int] = {}
    rule_status_counts: dict[str, int] = {}
    for record in records:
        confidence_counts[record.confidence] = confidence_counts.get(record.confidence, 0) + 1
        rule_status_counts[record.rule_status] = rule_status_counts.get(record.rule_status, 0) + 1
    total = len(records)
    high = confidence_counts.get(HIGH_CONFIDENCE, 0)
    return {
        "total_records": total,
        "high_confidence_count": high,
        "medium_confidence_count": confidence_counts.get(MEDIUM_CONFIDENCE, 0),
        "low_confidence_count": confidence_counts.get(LOW_CONFIDENCE, 0),
        "unresolved_count": confidence_counts.get(UNRESOLVED_CONFIDENCE, 0),
        "display_only_count": total - (len(overlay.entries) if overlay is not None else 0),
        "executable_overlay_entry_count": len(overlay.entries) if overlay is not None else 0,
        "high_confidence_fraction": round(high / total, 4) if total else 0.0,
        "confidence_counts": confidence_counts,
        "rule_status_counts": rule_status_counts,
    }


def load_gene_rule_evidence_cache(
    cache_path: Path | str | None = None,
) -> dict[str, GeneRuleEvidence]:
    path = Path(cache_path) if cache_path is not None else DEFAULT_GENE_RULE_EVIDENCE_CACHE
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (FileNotFoundError, OSError, json.JSONDecodeError, UnicodeDecodeError):
        return {
            record.common_name: record
            for record in (_seed_evidence_for_name(name) for name in DEFAULT_TARGET_NAMES)
        }
    records = payload.get("records") if isinstance(payload, dict) else None
    if not isinstance(records, list):
        return {}
    result: dict[str, GeneRuleEvidence] = {}
    for item in records:
        if not isinstance(item, dict):
            continue
        record = gene_rule_evidence_from_mapping(item)
        if record.common_name:
            normalised_name = record.common_name.upper()
            result[normalised_name] = replace(record, common_name=normalised_name)
    return result


def gene_rule_evidence_from_mapping(data: dict[str, object]) -> GeneRuleEvidence:
    return GeneRuleEvidence(
        common_name=str(data.get("common_name") or ""),
        candidate_locus_tag=str(data.get("candidate_locus_tag") or ""),
        external_ids={str(k): str(v) for k, v in dict(data.get("external_ids") or {}).items()},
        protein_name=str(data.get("protein_name") or ""),
        evidence_sources=tuple(str(item) for item in data.get("evidence_sources") or ()),
        confidence=str(data.get("confidence") or UNRESOLVED_CONFIDENCE),
        mapped_model_gene_present=bool(data.get("mapped_model_gene_present", False)),
        target_reaction_ids=tuple(str(item) for item in data.get("target_reaction_ids") or ()),
        proposed_gr_rule=str(data.get("proposed_gr_rule") or ""),
        proposed_rule=str(data.get("proposed_rule") or ""),
        rule_status=str(data.get("rule_status") or "not_executable"),
        recommended_action=str(data.get("recommended_action") or "manual_review_required"),
        notes=tuple(str(item) for item in data.get("notes") or ()),
    )


def _fetch_source_rows_for_records(
    records: list[GeneRuleEvidence],
    *,
    enable_online: bool,
    timeout_seconds: int,
    progress: Any | None,
) -> dict[str, dict[str, object]]:
    if not enable_online:
        return {}
    loci = tuple(
        dict.fromkeys(record.candidate_locus_tag for record in records if record.candidate_locus_tag)
    )
    if not loci:
        return {}
    if progress:
        progress(f"Fetching exact UniProt/KEGG evidence for {len(loci)} candidate loci...")
    uniprot_by_locus = _fetch_uniprot_rows_by_locus(loci, timeout_seconds=timeout_seconds)
    kegg_by_locus = _fetch_kegg_descriptions_by_locus(DEFAULT_KEGG_ORGANISM, timeout_seconds=timeout_seconds)
    return {
        locus: {
            "uniprot": uniprot_by_locus.get(locus),
            "kegg_description": kegg_by_locus.get(locus),
        }
        for locus in loci
    }


def _fetch_uniprot_rows_by_locus(
    loci: tuple[str, ...],
    *,
    timeout_seconds: int,
    chunk_size: int = 25,
) -> dict[str, dict[str, str]]:
    fields = "accession,gene_names,gene_primary,protein_name,xref_kegg,xref_geneid,xref_refseq,ec"
    results: dict[str, dict[str, str]] = {}
    for chunk in _chunks(loci, chunk_size):
        query = " OR ".join(f"gene_exact:{locus}" for locus in chunk)
        if not query:
            continue
        params = urlencode(
            {
                "query": f"({query}) AND proteome:{DEFAULT_UNIPROT_PROTEOME}",
                "format": "tsv",
                "fields": fields,
                "size": str(max(10, len(chunk) * 2)),
            }
        )
        try:
            text = _read_text_url(
                f"https://rest.uniprot.org/uniprotkb/search?{params}",
                timeout_seconds=timeout_seconds,
            )
        except OSError:
            continue
        for row in _tsv_dict_rows(text):
            for gene_name in _split_external_list(row.get("Gene Names", "")):
                if gene_name in chunk and gene_name not in results:
                    results[gene_name] = row
    return results


def _fetch_kegg_descriptions_by_locus(
    organism: str,
    *,
    timeout_seconds: int,
) -> dict[str, str]:
    try:
        text = _read_text_url(f"https://rest.kegg.jp/list/{organism}", timeout_seconds=timeout_seconds)
    except OSError:
        return {}
    prefix = f"{organism}:"
    descriptions: dict[str, str] = {}
    for line in text.splitlines():
        if "\t" not in line:
            continue
        key, description = line.split("\t", 1)
        locus = key[len(prefix):] if key.startswith(prefix) else key
        descriptions[locus] = description.strip()
    return descriptions


def _enrich_record_with_source_evidence(
    record: GeneRuleEvidence,
    source_rows: dict[str, dict[str, object]],
) -> GeneRuleEvidence:
    if not record.candidate_locus_tag:
        return record
    source = source_rows.get(record.candidate_locus_tag) or {}
    uniprot_row = source.get("uniprot") if isinstance(source.get("uniprot"), dict) else None
    kegg_description = str(source.get("kegg_description") or "")
    if not uniprot_row and not kegg_description:
        return record

    external_ids = dict(record.external_ids or {})
    evidence_sources = list(record.evidence_sources)
    if uniprot_row:
        external_ids.update(_external_ids_from_uniprot(uniprot_row))
        _append_unique(evidence_sources, "UniProt GS115 proteome exact locus")
        if uniprot_row.get("GeneID"):
            _append_unique(evidence_sources, "NCBI Gene cross-reference")
        if uniprot_row.get("RefSeq"):
            _append_unique(evidence_sources, "RefSeq cross-reference")
    if kegg_description:
        external_ids.setdefault("kegg", f"{DEFAULT_KEGG_ORGANISM}:{record.candidate_locus_tag}")
        _append_unique(evidence_sources, "KEGG ppa exact locus")

    protein_name = _clean_protein_name(str(uniprot_row.get("Protein names", ""))) if uniprot_row else ""
    protein_name = protein_name or record.protein_name or kegg_description
    confidence = record.confidence
    rule_status = record.rule_status
    recommended_action = record.recommended_action
    notes = list(record.notes)
    if (
        record.confidence == MEDIUM_CONFIDENCE
        and uniprot_row
        and kegg_description
        and _record_description_matches_common_name(record.common_name, protein_name, kegg_description)
    ):
        confidence = HIGH_CONFIDENCE
        rule_status = "high_confidence_locus_candidate"
        recommended_action = "eligible_for_overlay_if_all_complex_subunits_are_confirmed"
    elif record.confidence != HIGH_CONFIDENCE:
        rule_status = record.rule_status or "display_only_requires_manual_review"
        recommended_action = record.recommended_action or "manual_review_required"
    if record.common_name == "KAR2" and "SIL1" in protein_name.upper():
        confidence = LOW_CONFIDENCE
        rule_status = "display_only_name_context_not_exact_kar2_locus"
        recommended_action = "manual_locus_review_required"
        _append_unique(notes, "Exact locus evidence points to SIL1/Kar2p nucleotide exchange factor, not KAR2 itself.")
    return replace(
        record,
        external_ids=external_ids,
        protein_name=protein_name,
        evidence_sources=tuple(evidence_sources),
        confidence=confidence,
        rule_status=rule_status,
        recommended_action=recommended_action,
        notes=tuple(notes),
    )


def _record_description_matches_common_name(common_name: str, protein_name: str, kegg_description: str) -> bool:
    text = f"{protein_name} {kegg_description}".upper()
    name = common_name.upper()
    if name == "PDI1":
        return "DISULFIDE" in text or "PDI" in text
    if name == "ERO1":
        return "THIOL OXIDASE" in text or "OXIDATIVE PROTEIN FOLDING" in text
    return name in text


def _external_ids_from_uniprot(row: dict[str, str]) -> dict[str, str]:
    external_ids = {"uniprot": str(row.get("Entry") or "").strip()}
    if row.get("KEGG"):
        external_ids["kegg"] = _strip_terminal_semicolon(row["KEGG"])
    if row.get("GeneID"):
        external_ids["ncbi_gene"] = _strip_terminal_semicolon(row["GeneID"])
    if row.get("RefSeq"):
        external_ids["refseq"] = _strip_terminal_semicolon(row["RefSeq"])
    return {key: value for key, value in external_ids.items() if value}


def _read_text_url(url: str, *, timeout_seconds: int) -> str:
    request = Request(url, headers={"User-Agent": "pcSecPichia-gene-rule-evidence-cache/1.0"})
    with urlopen(request, timeout=timeout_seconds) as response:
        return response.read().decode("utf-8")


def _tsv_dict_rows(text: str) -> list[dict[str, str]]:
    if not text.strip():
        return []
    lines = text.splitlines()
    header = lines[0].split("\t")
    rows: list[dict[str, str]] = []
    for line in lines[1:]:
        values = line.split("\t")
        rows.append({key: values[index] if index < len(values) else "" for index, key in enumerate(header)})
    return rows


def _clean_protein_name(value: str) -> str:
    return str(value or "").split(" (", 1)[0].strip()


def _split_external_list(value: str) -> list[str]:
    return [item.strip() for item in str(value or "").replace(";", " ").split() if item.strip()]


def _strip_terminal_semicolon(value: str) -> str:
    return str(value or "").strip().rstrip(";")


def _append_unique(values: list[str], value: str) -> None:
    if value and value not in values:
        values.append(value)


def build_gpr_overlay(
    model: Any,
    evidence: dict[str, GeneRuleEvidence],
) -> GprOverlayResult:
    evidence = _normalise_evidence_by_name(evidence)
    reaction_ids = {str(item) for item in getattr(model, "rxns", ())}
    reactions_with_gpr = _reactions_with_existing_gpr(model)
    entries: list[GprOverlayEntry] = []
    skipped: list[dict[str, object]] = []
    warnings: list[str] = []

    pdi_entry = _complex_overlay_entry(
        PDI_COMPLEX_REACTION,
        PDI_COMPLEX_NAMES,
        evidence,
        reaction_ids,
        reactions_with_gpr,
    )
    if pdi_entry is not None:
        entries.append(pdi_entry)
    elif any(name in evidence for name in PDI_COMPLEX_NAMES):
        skipped.append(
            {
                "common_name": "/".join(PDI_COMPLEX_NAMES),
                "reaction_id": PDI_COMPLEX_REACTION,
                "reason": _complex_overlay_skip_reason(
                    PDI_COMPLEX_REACTION,
                    PDI_COMPLEX_NAMES,
                    evidence,
                    reaction_ids,
                    reactions_with_gpr,
                ),
            }
        )

    handled_names = set(PDI_COMPLEX_NAMES)
    for record in evidence.values():
        if record.common_name.upper() in handled_names:
            continue
        for reaction_id in record.target_reaction_ids:
            if reaction_id == PDI_COMPLEX_REACTION:
                skipped.append(_skip_row(record, reaction_id, "protected_complex_reaction_requires_complex_rule"))
                continue
            if reaction_id not in reaction_ids:
                skipped.append(_skip_row(record, reaction_id, "target_reaction_missing_in_model"))
                continue
            if reaction_id in reactions_with_gpr:
                skipped.append(_skip_row(record, reaction_id, "target_reaction_already_has_model_gpr_rule"))
                continue
            if not _is_executable_evidence(record):
                skipped.append(_skip_row(record, reaction_id, "confidence_below_executable_threshold"))
                continue
            entries.append(
                GprOverlayEntry(
                    reaction_id=reaction_id,
                    gene_locus_tags=(record.candidate_locus_tag,),
                    common_names=(record.common_name,),
                    proposed_gr_rule=record.candidate_locus_tag,
                    proposed_rule="",
                    rule_status="overlay_executable",
                    recommended_action="enable_only_for_explicit_analysis",
                )
            )

    if entries:
        warnings.append(
            "External evidence GPR overlay is experimental and is not part of the original MATLAB/model GPR."
        )
    return GprOverlayResult(entries=tuple(entries), skipped=tuple(skipped), warnings=tuple(warnings))


def apply_gpr_overlay_for_analysis(model: Any, overlay: GprOverlayResult) -> Any:
    """Return a model copy with overlay GPR rules for explicit analysis only."""
    if not overlay.entries:
        return copy.deepcopy(model)
    genes = [str(item) for item in getattr(model, "genes", ())]
    rxns = [str(item) for item in getattr(model, "rxns", ())]
    rules = list(getattr(model, "rules", ()) or [])
    gr_rules = list(getattr(model, "gr_rules", ()) or [])
    rules = _pad_list(rules, len(rxns), "")
    gr_rules = _pad_list(gr_rules, len(rxns), "")

    for entry in overlay.entries:
        for locus_tag in entry.gene_locus_tags:
            if locus_tag not in genes:
                genes.append(locus_tag)
        indices = [genes.index(locus_tag) + 1 for locus_tag in entry.gene_locus_tags]
        rule = " & ".join(f"x({index})" for index in indices)
        gr_rule = " and ".join(entry.gene_locus_tags)
        try:
            rxn_index = rxns.index(entry.reaction_id)
        except ValueError:
            continue
        rules[rxn_index] = rule
        gr_rules[rxn_index] = gr_rule

    try:
        return replace(model, genes=genes, rules=rules, gr_rules=gr_rules)
    except TypeError:
        copied = copy.deepcopy(model)
        copied.genes = genes
        copied.rules = rules
        copied.gr_rules = gr_rules
        if hasattr(copied, "gene_index") and isinstance(getattr(copied, "gene_index"), dict):
            copied.gene_index = {gene: idx for idx, gene in enumerate(genes)}
        return copied


def overlay_aliases_for_executable_rules(
    evidence: dict[str, GeneRuleEvidence],
    overlay: GprOverlayResult | None = None,
) -> dict[str, str]:
    """Return common-name aliases that are safe to use for explicit overlay analysis."""
    evidence = _normalise_evidence_by_name(evidence)
    executable_loci: set[str] = set()
    if overlay is not None:
        executable_loci = {
            locus_tag
            for entry in overlay.entries
            for locus_tag in entry.gene_locus_tags
            if locus_tag
        }
    if not executable_loci:
        return {}
    aliases: dict[str, str] = {}
    for record in evidence.values():
        if not _is_executable_evidence(record):
            continue
        if record.candidate_locus_tag not in executable_loci:
            continue
        aliases[record.common_name.upper()] = record.candidate_locus_tag
    return aliases


def _complex_overlay_entry(
    reaction_id: str,
    names: tuple[str, ...],
    evidence: dict[str, GeneRuleEvidence],
    reaction_ids: set[str],
    reactions_with_gpr: set[str],
) -> GprOverlayEntry | None:
    if reaction_id not in reaction_ids:
        return None
    if reaction_id in reactions_with_gpr:
        return None
    records = [evidence.get(name) for name in names]
    if any(
        record is None
        or not _is_executable_evidence(record)
        or reaction_id not in record.target_reaction_ids
        for record in records
    ):
        return None
    locus_tags = tuple(record.candidate_locus_tag for record in records if record is not None)
    if len(locus_tags) != len(names):
        return None
    return GprOverlayEntry(
        reaction_id=reaction_id,
        gene_locus_tags=locus_tags,
        common_names=names,
        proposed_gr_rule=" and ".join(locus_tags),
        proposed_rule="",
        rule_status="overlay_executable_complex_rule",
        recommended_action="enable_only_for_explicit_analysis",
    )


def _complex_overlay_skip_reason(
    reaction_id: str,
    names: tuple[str, ...],
    evidence: dict[str, GeneRuleEvidence],
    reaction_ids: set[str],
    reactions_with_gpr: set[str],
) -> str:
    if reaction_id not in reaction_ids:
        return "target_reaction_missing_in_model"
    if reaction_id in reactions_with_gpr:
        return "target_reaction_already_has_model_gpr_rule"
    records = [evidence.get(name) for name in names]
    if any(record is None for record in records):
        return "shared_complex_requires_locus_for_all_subunits"
    if any(record and reaction_id not in record.target_reaction_ids for record in records):
        return "shared_complex_subunit_targets_different_reaction"
    return "shared_complex_requires_high_confidence_locus_for_all_subunits"


def _normalise_evidence_by_name(evidence: dict[str, GeneRuleEvidence]) -> dict[str, GeneRuleEvidence]:
    normalised: dict[str, GeneRuleEvidence] = {}
    for key, record in evidence.items():
        name = str(record.common_name or key or "").strip().upper()
        if not name:
            continue
        normalised[name] = replace(record, common_name=name)
    return normalised


def _is_executable_evidence(record: GeneRuleEvidence | None) -> bool:
    return bool(
        record
        and record.confidence == HIGH_CONFIDENCE
        and record.candidate_locus_tag
        and record.target_reaction_ids
    )


def _skip_row(record: GeneRuleEvidence, reaction_id: str, reason: str) -> dict[str, object]:
    return {
        "common_name": record.common_name,
        "candidate_locus_tag": record.candidate_locus_tag,
        "reaction_id": reaction_id,
        "confidence": record.confidence,
        "reason": reason,
    }


def _reactions_with_existing_gpr(model: Any) -> set[str]:
    rxns = [str(item) for item in getattr(model, "rxns", ())]
    rules = list(getattr(model, "rules", ()) or [])
    gr_rules = list(getattr(model, "gr_rules", ()) or [])
    reactions: set[str] = set()
    for idx, reaction_id in enumerate(rxns):
        rule = str(rules[idx] if idx < len(rules) else "").strip()
        gr_rule = str(gr_rules[idx] if idx < len(gr_rules) else "").strip()
        if rule and rule != "[]" or gr_rule and gr_rule != "[]":
            reactions.add(reaction_id)
    return reactions


def _seed_evidence_for_name(name: str) -> GeneRuleEvidence:
    seeds = _seed_evidence()
    key = str(name or "").strip().upper()
    if key in seeds:
        return seeds[key]
    return GeneRuleEvidence(common_name=key or str(name or ""), confidence=UNRESOLVED_CONFIDENCE)


def _seed_evidence() -> dict[str, GeneRuleEvidence]:
    return {
        "PDI1": GeneRuleEvidence(
            common_name="PDI1",
            candidate_locus_tag="PAS_chr1-1_0160",
            external_ids={"kegg": "ppa:PAS_chr1-1_0160"},
            protein_name="Protein disulfide isomerase family protein",
            evidence_sources=("KEGG ppa function search", "model reaction name evidence"),
            confidence=MEDIUM_CONFIDENCE,
            target_reaction_ids=(PDI_COMPLEX_REACTION,),
            rule_status="display_only_requires_multi_source_confirmation",
            recommended_action="keep_reaction_level_proxy_until_locus_is_confirmed",
            notes=("KEGG also lists other PDI-family candidates; do not treat as final PDI1 locus without review.",),
        ),
        "ERO1": GeneRuleEvidence(
            common_name="ERO1",
            candidate_locus_tag="PAS_chr1-1_0011",
            external_ids={"kegg": "ppa:PAS_chr1-1_0011"},
            protein_name="Thiol oxidase required for oxidative protein folding in the ER",
            evidence_sources=("KEGG ppa find ERO1", "model reaction name evidence"),
            confidence=MEDIUM_CONFIDENCE,
            target_reaction_ids=(PDI_COMPLEX_REACTION,),
            rule_status="display_only_requires_multi_source_confirmation",
            recommended_action="keep_reaction_level_proxy_until_locus_is_confirmed",
        ),
        "ERV2": GeneRuleEvidence(
            common_name="ERV2",
            candidate_locus_tag="",
            evidence_sources=("KEGG ppa name search",),
            confidence=UNRESOLVED_CONFIDENCE,
            target_reaction_ids=(PDI_COMPLEX_REACTION,),
            rule_status="not_executable",
            recommended_action="manual_locus_review_required",
            notes=("KEGG ERV2 name search returned ERV25-like p24 family evidence, not a reliable ERV2 locus.",),
        ),
        "KAR2": GeneRuleEvidence(
            common_name="KAR2",
            candidate_locus_tag="PAS_chr1-1_0237",
            external_ids={"kegg": "ppa:PAS_chr1-1_0237"},
            protein_name="Nucleotide exchange factor for ER lumenal Hsp70 chaperone Kar2p",
            evidence_sources=("KEGG ppa find KAR2",),
            confidence=LOW_CONFIDENCE,
            target_reaction_ids=("sec_Kar2p_complex_formation",),
            rule_status="display_only_name_context_not_exact_kar2_locus",
            recommended_action="manual_locus_review_required",
        ),
        "OCH1": GeneRuleEvidence(
            common_name="OCH1",
            confidence=UNRESOLVED_CONFIDENCE,
            target_reaction_ids=("sec_Och1p_complex_formation",),
            rule_status="not_executable",
            recommended_action="manual_locus_review_required",
        ),
        "PEP4": GeneRuleEvidence(
            common_name="PEP4",
            candidate_locus_tag="PAS_chr3_1087",
            external_ids={"kegg": "ppa:PAS_chr3_1087"},
            protein_name="Vacuolar aspartyl protease / proteinase A candidate",
            evidence_sources=("KEGG ppa function search",),
            confidence=LOW_CONFIDENCE,
            rule_status="display_only_conflicts_with_existing_model_gene_annotation",
            recommended_action="do_not_replace_existing_model_gene_without_review",
            notes=("Existing curated model KO uses PAS_chr2-2_0107; external function search suggests a different locus.",),
        ),
        "PRB1": GeneRuleEvidence(
            common_name="PRB1",
            candidate_locus_tag="",
            evidence_sources=("KEGG ppa find PRB1", "KEGG ppa vacuolar protease search"),
            confidence=LOW_CONFIDENCE,
            rule_status="display_only_multiple_or_indirect_candidates",
            recommended_action="manual_locus_review_required",
            notes=("Search results include indirect 'processed by Prb1p' or multiple proteinase B candidates.",),
        ),
    }


def _dedupe_names(names: tuple[str, ...]) -> tuple[str, ...]:
    result: list[str] = []
    for name in names:
        text = str(name or "").strip().upper()
        if text and text not in result:
            result.append(text)
    return tuple(result)


def _chunks(values: tuple[str, ...], size: int) -> list[tuple[str, ...]]:
    return [values[index:index + size] for index in range(0, len(values), size)]


def _pad_list(values: list[Any], length: int, fill: Any) -> list[Any]:
    if len(values) >= length:
        return values
    return [*values, *([fill] * (length - len(values)))]


def _write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")


def _write_report(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _build_gene_rule_evidence_report(
    records: list[GeneRuleEvidence],
    summary: dict[str, object],
    overlay: GprOverlayResult | None,
) -> str:
    lines = [
        "# 毕赤酵母分泌工程 GPR 证据缓存报告",
        "",
        "## 结论",
        "",
        "- 这是 Python 侧 evidence overlay，不是原始 MATLAB/model GPR。",
        "- 默认 pipeline 不启用 overlay；只有用户显式勾选后才用于 KO/OE 预检和小批量分析。",
        "- 只有高置信度、目标反应存在、且不覆盖原始模型 GPR 的条目才可能变成可执行 overlay。",
        "- PDI1/ERO1/ERV2 共享复合体反应，需要三个亚基 locus 都确认后才允许生成复合体 GPR。",
        "",
        "## 覆盖率摘要",
        "",
        f"- total records: `{summary.get('total_records')}`",
        f"- high confidence: `{summary.get('high_confidence_count')}`",
        f"- medium confidence: `{summary.get('medium_confidence_count')}`",
        f"- low confidence: `{summary.get('low_confidence_count')}`",
        f"- unresolved: `{summary.get('unresolved_count')}`",
        f"- executable overlay entries: `{summary.get('executable_overlay_entry_count')}`",
        "",
        "## 记录明细",
        "",
        "| common name | candidate locus | confidence | rule status | reactions | sources | recommended action |",
        "| --- | --- | --- | --- | --- | --- | --- |",
    ]
    for record in records:
        lines.append(
            "| "
            + " | ".join(
                [
                    record.common_name,
                    record.candidate_locus_tag,
                    record.confidence,
                    record.rule_status,
                    ", ".join(record.target_reaction_ids),
                    ", ".join(record.evidence_sources),
                    record.recommended_action,
                ]
            )
            + " |"
        )
    if overlay is not None:
        lines.extend(
            [
                "",
                "## Overlay Probe",
                "",
                f"- executable entries: `{len(overlay.entries)}`",
                f"- skipped rows: `{len(overlay.skipped)}`",
            ]
        )
        if overlay.skipped:
            lines.extend(["", "| common name | reaction | reason |", "| --- | --- | --- |"])
            for row in overlay.skipped:
                lines.append(
                    f"| {row.get('common_name', '')} | {row.get('reaction_id', '')} | {row.get('reason', '')} |"
                )
    lines.extend(
        [
            "",
            "## 数据源",
            "",
            f"- UniProt GS115 proteome: `{DEFAULT_UNIPROT_PROTEOME}`",
            f"- KEGG organism: `{DEFAULT_KEGG_ORGANISM}`",
            "- 常用名搜索结果只作为人工复核线索；自动升级只使用 locus tag 精确匹配证据。",
        ]
    )
    return "\n".join(lines) + "\n"


__all__ = [
    "DEFAULT_GENE_RULE_EVIDENCE_CACHE",
    "DEFAULT_GENE_RULE_EVIDENCE_REPORT",
    "DEFAULT_GENE_RULE_EVIDENCE_SUMMARY",
    "GeneRuleEvidence",
    "GprOverlayEntry",
    "GprOverlayResult",
    "apply_gpr_overlay_for_analysis",
    "build_gene_rule_evidence_cache",
    "build_gpr_overlay",
    "gene_rule_evidence_from_mapping",
    "load_gene_rule_evidence_cache",
    "overlay_aliases_for_executable_rules",
    "summarize_gene_rule_evidence_records",
]
