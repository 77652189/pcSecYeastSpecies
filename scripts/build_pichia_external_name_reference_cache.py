from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict, dataclass, replace
from datetime import datetime
from pathlib import Path
from typing import Callable

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "python_pichia" / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "python_pichia" / "src"))

from pcsec_pichia.homology.cache_schema import ExternalNameReference, NameAuditRow
from pcsec_pichia.homology.crosswalk import load_external_name_reference_cache, load_name_audit_cache
from pcsec_pichia.homology.external_fetch import (
    ExternalFetchConfig,
    ExternalFetchResult,
    fetch_external_name_references,
)


NAME_AUDIT_JSONL = "sce_to_pichia_name_audit.jsonl"
EXTERNAL_REFERENCE_JSONL = "external_name_references.jsonl"
EXTERNAL_REFERENCE_SUMMARY_MD = "external_name_reference_summary.md"


@dataclass(frozen=True)
class ExternalReferenceQuery:
    query: str
    match_key: str
    internal_common_name: str
    internal_sequence_id: str


FetchMany = Callable[[str, ExternalFetchConfig], tuple[ExternalFetchResult, ...]]


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    root = Path(args.root).resolve()
    name_audit_path = _resolve_name_audit_path(args, root)
    if not name_audit_path.exists():
        print(f"Name audit cache not found: {name_audit_path}", file=sys.stderr)
        return 2

    output_path = _resolve_output_path(args.output_path, root)
    summary_path = output_path.with_name(EXTERNAL_REFERENCE_SUMMARY_MD)
    if output_path.exists() and not args.overwrite and not args.dry_run:
        print(f"Output already exists; pass --overwrite to replace: {output_path}", file=sys.stderr)
        return 2

    name_rows = load_name_audit_cache(name_audit_path)
    sources = _parse_sources(args.sources)
    config = ExternalFetchConfig.from_env(
        timeout_seconds=args.timeout,
        retry_count=args.retry_count,
        delay_seconds=args.delay,
        enabled_sources=sources,
    )
    queries = collect_external_reference_queries(name_rows, limit=args.limit)
    if args.dry_run:
        print(
            json.dumps(
                {
                    "name_audit_path": str(name_audit_path),
                    "output_path": str(output_path),
                    "query_count": len(queries),
                    "queries": [asdict(query) for query in queries],
                },
                ensure_ascii=False,
                indent=2,
            )
        )
        return 0

    fixture_references = _load_offline_fixture(args.offline_fixture, root)
    references, summary = build_external_reference_cache(
        name_rows,
        config=config,
        sources=sources,
        limit=args.limit,
        fixture_references=fixture_references,
    )
    write_external_reference_cache(references, output_path)
    write_external_reference_summary(summary_path, summary, output_path=output_path)
    print(
        json.dumps(
            {
                "output_path": str(output_path),
                "summary_path": str(summary_path),
                "reference_count": len(references),
                "query_count": summary["query_count"],
                "sources": sorted(summary["source_counts"]),
                "warning_count": len(summary["warnings"]),
            },
            ensure_ascii=False,
            indent=2,
        )
    )
    return 0


def collect_external_reference_queries(
    name_rows: tuple[NameAuditRow, ...],
    *,
    limit: int = 0,
) -> tuple[ExternalReferenceQuery, ...]:
    queries: list[ExternalReferenceQuery] = []
    seen: set[str] = set()
    for row in name_rows:
        for match_key, value in (
            ("external_accession", row.external_accession),
            ("external_locus_tag", row.external_locus_tag),
            ("internal_gene_id", row.internal_gene_id),
        ):
            query = str(value or "").strip()
            if not query or query.lower() in seen:
                continue
            seen.add(query.lower())
            queries.append(
                ExternalReferenceQuery(
                    query=query,
                    match_key=match_key,
                    internal_common_name=row.internal_common_name,
                    internal_sequence_id=row.internal_sequence_id,
                )
            )
            if limit and len(queries) >= limit:
                return tuple(queries)
    return tuple(queries)


def build_external_reference_cache(
    name_rows: tuple[NameAuditRow, ...],
    *,
    config: ExternalFetchConfig,
    sources: tuple[str, ...],
    limit: int = 0,
    fetch_many: FetchMany = fetch_external_name_references,
    fixture_references: tuple[ExternalNameReference, ...] = (),
) -> tuple[tuple[ExternalNameReference, ...], dict[str, object]]:
    queries = collect_external_reference_queries(name_rows, limit=limit)
    fetch_config = replace(config, enabled_sources=sources)
    references: list[ExternalNameReference] = []
    warnings: list[str] = []
    source_counts: dict[str, int] = {}
    success_count = 0
    failure_count = 0
    seen_refs: set[tuple[str, str, str, str]] = set()

    for query in queries:
        results = (
            _fixture_results_for_query(query, sources, fixture_references)
            if fixture_references
            else fetch_many(query.query, fetch_config)
        )
        query_success = False
        for result in results:
            if result.warnings:
                warnings.extend(f"{result.source_database}:{query.query}: {warning}" for warning in result.warnings)
            if not result.success:
                continue
            query_success = True
            for reference in result.references:
                enriched = _with_query_warning(reference, query)
                ref_key = _reference_key(enriched)
                if ref_key in seen_refs:
                    continue
                seen_refs.add(ref_key)
                references.append(enriched)
                source_counts[enriched.source_database] = source_counts.get(enriched.source_database, 0) + 1
        if query_success:
            success_count += 1
        else:
            failure_count += 1
            warnings.append(f"no reference built for {query.match_key}={query.query}")

    summary: dict[str, object] = {
        "query_count": len(queries),
        "reference_count": len(references),
        "success_count": success_count,
        "failure_count": failure_count,
        "source_counts": dict(sorted(source_counts.items())),
        "warnings": warnings,
    }
    return tuple(sorted(references, key=_reference_key)), summary


def write_external_reference_cache(rows: tuple[ExternalNameReference, ...], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="\n") as handle:
        for row in rows:
            payload = asdict(row)
            payload["aliases"] = list(row.aliases)
            payload["warnings"] = list(row.warnings)
            handle.write(json.dumps(payload, ensure_ascii=False, sort_keys=True) + "\n")


def write_external_reference_summary(summary_path: Path, summary: dict[str, object], *, output_path: Path) -> None:
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# External Name Reference Cache Summary",
        "",
        f"- output_path: `{output_path}`",
        f"- query_count: {summary['query_count']}",
        f"- reference_count: {summary['reference_count']}",
        f"- success_count: {summary['success_count']}",
        f"- failure_count: {summary['failure_count']}",
        "",
        "## Source Counts",
        "",
    ]
    source_counts = summary.get("source_counts")
    if isinstance(source_counts, dict) and source_counts:
        for source, count in source_counts.items():
            lines.append(f"- {source}: {count}")
    else:
        lines.append("- none: 0")
    warnings = summary.get("warnings")
    if isinstance(warnings, list) and warnings:
        lines.extend(["", "## Warnings", ""])
        for warning in warnings:
            lines.append(f"- {warning}")
    summary_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _parse_args(argv: list[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build offline external name reference cache for Pichia homology audits.")
    parser.add_argument("--root", default=str(REPO_ROOT))
    parser.add_argument("--homology-cache-dir", default="")
    parser.add_argument("--name-audit-jsonl", default="")
    parser.add_argument("--output-path", default="")
    parser.add_argument("--sources", default="uniprot,ncbi,sgd")
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--timeout", type=float, default=10.0)
    parser.add_argument("--delay", type=float, default=0.34)
    parser.add_argument("--retry-count", type=int, default=1)
    parser.add_argument("--offline-fixture", default="")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args(argv)


def _resolve_name_audit_path(args: argparse.Namespace, root: Path) -> Path:
    if args.name_audit_jsonl:
        path = Path(args.name_audit_jsonl)
        return path if path.is_absolute() else root / path
    if args.homology_cache_dir:
        path = Path(args.homology_cache_dir)
        run_dir = path if path.is_absolute() else root / path
        return run_dir / NAME_AUDIT_JSONL
    return root / "local_runs" / "pichia_homology_cache" / "smoke" / NAME_AUDIT_JSONL


def _resolve_output_path(value: str, root: Path) -> Path:
    if value:
        path = Path(value)
        resolved = path if path.is_absolute() else root / path
        return resolved if resolved.suffix else resolved / EXTERNAL_REFERENCE_JSONL
    return root / "local_runs" / "pichia_homology_cache" / _run_name() / EXTERNAL_REFERENCE_JSONL


def _parse_sources(value: str) -> tuple[str, ...]:
    sources = tuple(source.strip().lower() for source in value.split(",") if source.strip())
    return sources or ("uniprot", "ncbi", "sgd")


def _load_offline_fixture(value: str, root: Path) -> tuple[ExternalNameReference, ...]:
    if not value:
        return ()
    path = Path(value)
    resolved = path if path.is_absolute() else root / path
    return load_external_name_reference_cache(resolved)


def _fixture_results_for_query(
    query: ExternalReferenceQuery,
    sources: tuple[str, ...],
    fixture_references: tuple[ExternalNameReference, ...],
) -> tuple[ExternalFetchResult, ...]:
    results: list[ExternalFetchResult] = []
    for source in sources:
        refs = tuple(
            reference
            for reference in fixture_references
            if reference.source_database.lower() == source.lower() and _reference_matches_query(reference, query.query)
        )
        results.append(
            ExternalFetchResult(
                source_database=source,
                query=query.query,
                success=bool(refs),
                references=refs,
                warnings=() if refs else (f"offline fixture had no {source} reference for `{query.query}`",),
                error_summary="" if refs else "no_fixture_record",
            )
        )
    return tuple(results)


def _reference_matches_query(reference: ExternalNameReference, query: str) -> bool:
    normalized = _normalize(query)
    values = (
        reference.accession,
        reference.gene_name,
        reference.locus_tag,
        *reference.aliases,
    )
    return normalized in {_normalize(value) for value in values if value}


def _with_query_warning(reference: ExternalNameReference, query: ExternalReferenceQuery) -> ExternalNameReference:
    note = f"query={query.query}; match_key={query.match_key}"
    return replace(reference, warnings=tuple(dict.fromkeys((*reference.warnings, note))))


def _reference_key(reference: ExternalNameReference) -> tuple[str, str, str, str]:
    return (
        reference.source_database,
        reference.accession,
        reference.gene_name,
        reference.locus_tag,
    )


def _normalize(value: object) -> str:
    return str(value or "").strip().lower()


def _run_name() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")


if __name__ == "__main__":
    raise SystemExit(main())
