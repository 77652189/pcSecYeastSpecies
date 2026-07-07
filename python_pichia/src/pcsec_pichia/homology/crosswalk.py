from __future__ import annotations

import csv
import json
from dataclasses import asdict
from pathlib import Path
from typing import Any

from pcsec_pichia.homology.cache_schema import (
    BlastConfig,
    CatalogHomologyQuery,
    CacheWriteResult,
    HomologyAuditSummary,
    HomologyCrosswalkRow,
    NameAuditRow,
    ProteinRecord,
    ReciprocalBestHit,
    RuleTransferAuditRow,
)
from pcsec_pichia.homology.review_rules import (
    classify_homology_review_status,
    classify_name_consistency,
    classify_rule_transfer_status,
)


def build_homology_crosswalk(
    sce_queries: tuple[CatalogHomologyQuery, ...],
    sce_records: tuple[ProteinRecord, ...],
    pichia_records: tuple[ProteinRecord, ...],
    model_gene_index: set[str],
    rbh_calls: tuple[ReciprocalBestHit, ...],
    blast_config: BlastConfig,
) -> tuple[HomologyCrosswalkRow, ...]:
    """Build the offline crosswalk from local sequence assets and RBH calls."""

    sce_by_symbol = _sce_symbol_lookup(sce_records)
    rbh_by_query = {call.query_id: call for call in rbh_calls}
    pichia_by_gene = {record.gene_id: record for record in pichia_records}
    rows: list[HomologyCrosswalkRow] = []
    for query in sce_queries:
        sce_record = _resolve_sce_record(query, sce_by_symbol)
        if sce_record is None:
            status, warnings = classify_homology_review_status(
                is_rbh=False,
                identity_pct=None,
                query_coverage=None,
                subject_coverage=None,
                in_model_gene_index=False,
                min_identity=blast_config.min_identity,
                min_coverage=blast_config.min_coverage,
                unresolved_query=True,
            )
            rows.append(
                HomologyCrosswalkRow(
                    internal_common_name=query.internal_common_name,
                    query_symbol=query.query_symbol,
                    sce_orf=None,
                    pichia_gene_id=None,
                    pichia_model_gene_id=None,
                    is_rbh=False,
                    identity_pct=None,
                    evalue=None,
                    query_coverage=None,
                    subject_coverage=None,
                    in_model_gene_index=False,
                    review_status=status,
                    warnings=warnings,
                )
            )
            continue
        call = rbh_by_query.get(sce_record.gene_id)
        hit = call.forward_hit if call else None
        pichia_gene_id = call.subject_id if call else None
        in_model = bool(pichia_gene_id and pichia_gene_id in model_gene_index)
        status, warnings = classify_homology_review_status(
            is_rbh=bool(call and call.is_rbh),
            identity_pct=hit.identity_pct if hit else None,
            query_coverage=hit.query_coverage if hit else None,
            subject_coverage=hit.subject_coverage if hit else None,
            in_model_gene_index=in_model,
            min_identity=blast_config.min_identity,
            min_coverage=blast_config.min_coverage,
        )
        if call and call.failure_reason:
            warnings = (*warnings, call.failure_reason)
        pichia_record = pichia_by_gene.get(pichia_gene_id or "")
        rows.append(
            HomologyCrosswalkRow(
                internal_common_name=query.internal_common_name,
                query_symbol=query.query_symbol,
                sce_orf=sce_record.gene_id,
                pichia_gene_id=pichia_gene_id,
                pichia_model_gene_id=pichia_gene_id if in_model else None,
                is_rbh=bool(call and call.is_rbh),
                identity_pct=hit.identity_pct if hit else None,
                evalue=hit.evalue if hit else None,
                query_coverage=hit.query_coverage if hit else None,
                subject_coverage=hit.subject_coverage if hit else None,
                in_model_gene_index=in_model,
                review_status=status,
                warnings=warnings,
                external_accession=pichia_record.accession if pichia_record and pichia_record.accession else "",
                external_gene_name=pichia_record.symbol if pichia_record and pichia_record.symbol else "",
                external_locus_tag=pichia_record.gene_id if pichia_record else "",
                external_aliases=pichia_record.aliases if pichia_record else (),
            )
        )
    return tuple(rows)


def build_name_audit_rows(crosswalk: tuple[HomologyCrosswalkRow, ...]) -> tuple[NameAuditRow, ...]:
    """Build name-audit rows while keeping model operability separate."""

    rows: list[NameAuditRow] = []
    for row in crosswalk:
        name_status = classify_name_consistency(
            internal_common_name=row.query_symbol,
            external_gene_name=row.external_gene_name,
            external_aliases=row.external_aliases,
            is_rbh=row.is_rbh,
        )
        rows.append(
            NameAuditRow(
                internal_gene_id=row.pichia_model_gene_id or "",
                internal_common_name=row.internal_common_name,
                internal_sequence_id=row.sce_orf or "",
                external_accession=row.external_accession,
                external_gene_name=row.external_gene_name,
                external_locus_tag=row.external_locus_tag,
                external_aliases=row.external_aliases,
                identity_pct=row.identity_pct,
                query_coverage=row.query_coverage,
                subject_coverage=row.subject_coverage,
                evalue=row.evalue,
                is_rbh=row.is_rbh,
                in_model_gene_index=row.in_model_gene_index,
                name_consistency_status=name_status,
                review_status=row.review_status,
                warnings=row.warnings,
            )
        )
    return tuple(rows)


def build_rule_transfer_audit_rows(crosswalk: tuple[HomologyCrosswalkRow, ...]) -> tuple[RuleTransferAuditRow, ...]:
    """Build rule-transfer audit rows without changing phenotype evidence."""

    rows: list[RuleTransferAuditRow] = []
    for row in crosswalk:
        status, warnings = classify_rule_transfer_status(
            homology_review_status=row.review_status,
            is_rbh=row.is_rbh,
            in_model_gene_index=row.in_model_gene_index,
        )
        rows.append(
            RuleTransferAuditRow(
                internal_common_name=row.internal_common_name,
                query_symbol=row.query_symbol,
                sce_orf=row.sce_orf or "",
                pichia_gene_id=row.pichia_gene_id or "",
                pichia_model_gene_id=row.pichia_model_gene_id or "",
                is_rbh=row.is_rbh,
                in_model_gene_index=row.in_model_gene_index,
                identity_pct=row.identity_pct,
                query_coverage=row.query_coverage,
                subject_coverage=row.subject_coverage,
                evalue=row.evalue,
                homology_review_status=row.review_status,
                rule_transfer_status=status,
                warnings=(*row.warnings, *warnings),
            )
        )
    return tuple(rows)


def summarize_homology_audits(
    *,
    blast_status: str,
    homology_rows: tuple[HomologyCrosswalkRow, ...],
    name_audit_rows: tuple[NameAuditRow, ...],
    rule_transfer_rows: tuple[RuleTransferAuditRow, ...],
) -> HomologyAuditSummary:
    return HomologyAuditSummary(
        blast_status=blast_status,
        homology_row_count=len(homology_rows),
        name_audit_row_count=len(name_audit_rows),
        rule_transfer_row_count=len(rule_transfer_rows),
        homology_review_status_counts=_count_by(homology_rows, "review_status"),
        name_consistency_status_counts=_count_by(name_audit_rows, "name_consistency_status"),
        rule_transfer_status_counts=_count_by(rule_transfer_rows, "rule_transfer_status"),
    )


def write_homology_cache(
    crosswalk: tuple[HomologyCrosswalkRow, ...],
    jsonl_path: Path,
    tsv_path: Path,
) -> CacheWriteResult:
    """Write deterministic JSONL and TSV cache outputs."""

    jsonl_path.parent.mkdir(parents=True, exist_ok=True)
    tsv_path.parent.mkdir(parents=True, exist_ok=True)
    ordered = sorted(crosswalk, key=lambda row: (row.internal_common_name, row.query_symbol, row.sce_orf or ""))
    with jsonl_path.open("w", encoding="utf-8", newline="\n") as handle:
        for row in ordered:
            handle.write(json.dumps(_json_row(row), ensure_ascii=False, sort_keys=True) + "\n")
    fieldnames = list(_tsv_row(ordered[0]).keys()) if ordered else list(_tsv_row(None).keys())
    with tsv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in ordered:
            writer.writerow(_tsv_row(row))
    return CacheWriteResult(jsonl_path=jsonl_path, tsv_path=tsv_path, row_count=len(ordered))


def write_name_audit_cache(
    rows: tuple[NameAuditRow, ...],
    jsonl_path: Path,
    tsv_path: Path,
) -> CacheWriteResult:
    ordered = tuple(sorted(rows, key=lambda row: (row.internal_common_name, row.internal_sequence_id)))
    return _write_dataclass_cache(ordered, jsonl_path, tsv_path)


def write_rule_transfer_audit_cache(
    rows: tuple[RuleTransferAuditRow, ...],
    jsonl_path: Path,
    tsv_path: Path,
) -> CacheWriteResult:
    ordered = tuple(sorted(rows, key=lambda row: (row.internal_common_name, row.query_symbol, row.sce_orf)))
    return _write_dataclass_cache(ordered, jsonl_path, tsv_path)


def load_homology_cache(path: Path) -> tuple[HomologyCrosswalkRow, ...]:
    """Load a previously generated JSONL cache without rerunning BLAST."""

    rows: list[HomologyCrosswalkRow] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            payload = json.loads(line)
            rows.append(
                HomologyCrosswalkRow(
                    **{
                        **payload,
                        "warnings": tuple(payload.get("warnings", ())),
                        "external_aliases": tuple(payload.get("external_aliases", ())),
                    }
                )
            )
    return tuple(rows)


def load_name_audit_cache(path: Path) -> tuple[NameAuditRow, ...]:
    rows: list[NameAuditRow] = []
    for payload in _read_jsonl(path):
        rows.append(
            NameAuditRow(
                **{
                    **payload,
                    "external_aliases": tuple(payload.get("external_aliases", ())),
                    "warnings": tuple(payload.get("warnings", ())),
                }
            )
        )
    return tuple(rows)


def load_rule_transfer_audit_cache(path: Path) -> tuple[RuleTransferAuditRow, ...]:
    rows: list[RuleTransferAuditRow] = []
    for payload in _read_jsonl(path):
        rows.append(
            RuleTransferAuditRow(
                **{
                    **payload,
                    "warnings": tuple(payload.get("warnings", ())),
                }
            )
        )
    return tuple(rows)


def _sce_symbol_lookup(records: tuple[ProteinRecord, ...]) -> dict[str, ProteinRecord]:
    lookup: dict[str, ProteinRecord] = {}
    for record in records:
        keys = [record.gene_id]
        if record.symbol:
            keys.append(record.symbol)
        keys.extend(record.aliases)
        for key in keys:
            lookup.setdefault(key.upper(), record)
    return lookup


def _resolve_sce_record(
    query: CatalogHomologyQuery,
    lookup: dict[str, ProteinRecord],
) -> ProteinRecord | None:
    for key in (query.query_symbol, *query.aliases):
        record = lookup.get(key.upper())
        if record:
            return record
    return None


def _json_row(row: HomologyCrosswalkRow) -> dict[str, object]:
    payload = asdict(row)
    payload["warnings"] = list(row.warnings)
    payload["external_aliases"] = list(row.external_aliases)
    return payload


def _tsv_row(row: HomologyCrosswalkRow | None) -> dict[str, object]:
    if row is None:
        return {
            "internal_common_name": "",
            "query_symbol": "",
            "sce_orf": "",
            "pichia_gene_id": "",
            "pichia_model_gene_id": "",
            "is_rbh": False,
            "identity_pct": "",
            "evalue": "",
            "query_coverage": "",
            "subject_coverage": "",
            "in_model_gene_index": False,
            "review_status": "",
            "warnings": "",
            "external_accession": "",
            "external_gene_name": "",
            "external_locus_tag": "",
            "external_aliases": "",
        }
    payload = asdict(row)
    payload["warnings"] = ";".join(row.warnings)
    payload["external_aliases"] = ";".join(row.external_aliases)
    return payload


def _write_dataclass_cache(rows: tuple[Any, ...], jsonl_path: Path, tsv_path: Path) -> CacheWriteResult:
    jsonl_path.parent.mkdir(parents=True, exist_ok=True)
    tsv_path.parent.mkdir(parents=True, exist_ok=True)
    with jsonl_path.open("w", encoding="utf-8", newline="\n") as handle:
        for row in rows:
            handle.write(json.dumps(_generic_json_row(row), ensure_ascii=False, sort_keys=True) + "\n")
    fieldnames = list(_generic_tsv_row(rows[0]).keys()) if rows else []
    with tsv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(_generic_tsv_row(row))
    return CacheWriteResult(jsonl_path=jsonl_path, tsv_path=tsv_path, row_count=len(rows))


def _generic_json_row(row: Any) -> dict[str, object]:
    payload = asdict(row)
    for key, value in list(payload.items()):
        if isinstance(value, tuple):
            payload[key] = list(value)
    return payload


def _generic_tsv_row(row: Any) -> dict[str, object]:
    payload = asdict(row)
    for key, value in list(payload.items()):
        if isinstance(value, tuple):
            payload[key] = ";".join(str(item) for item in value)
    return payload


def _read_jsonl(path: Path) -> list[dict[str, object]]:
    payloads: list[dict[str, object]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                payloads.append(json.loads(line))
    return payloads


def _count_by(rows: tuple[Any, ...], attribute: str) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in rows:
        value = str(getattr(row, attribute))
        counts[value] = counts.get(value, 0) + 1
    return dict(sorted(counts.items()))
