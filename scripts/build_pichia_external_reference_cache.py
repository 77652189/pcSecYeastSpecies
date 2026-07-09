from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
PYTHON_SRC = REPO_ROOT / "python_pichia" / "src"
if str(PYTHON_SRC) not in sys.path:
    sys.path.insert(0, str(PYTHON_SRC))

from pcsec_pichia.external_refs.clients import ExternalFetchConfig
from pcsec_pichia.external_refs.queries import (
    ExternalReferenceQuery,
    build_external_reference_queries,
)
from pcsec_pichia.external_refs.refresh import (
    FAILED_QUERIES_FILENAME,
    SUMMARY_FILENAME,
    build_external_reference_cache,
)
from pcsec_pichia.services.gene_catalog import SECRETION_GENE_CATALOG


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    output_dir = _resolve_output_dir(args.output_dir)
    queries = _collect_queries(args)
    if args.limit and args.limit > 0:
        queries = queries[: args.limit]
    sources = _parse_sources(args.sources)
    config = ExternalFetchConfig.from_env(
        sources=sources,
        timeout_seconds=args.timeout,
        retry_attempts=args.retry_attempts,
        min_interval_seconds=args.min_interval,
    )
    if args.dry_run:
        print(
            json.dumps(
                {
                    "output_dir": str(output_dir),
                    "sources": sources,
                    "query_count": len(queries),
                    "queries": [query.to_dict() for query in queries],
                },
                ensure_ascii=False,
                indent=2,
            )
        )
        return 0
    manifest = build_external_reference_cache(queries, output_dir, config=config)
    print(
        json.dumps(
            {
                "output_dir": str(output_dir),
                "records_path": str(output_dir / "external_reference_records.jsonl"),
                "manifest_path": str(output_dir / "external_reference_manifest.json"),
                "failed_queries_path": str(output_dir / FAILED_QUERIES_FILENAME),
                "summary_path": str(output_dir / SUMMARY_FILENAME),
                "query_count": manifest.query_count,
                "record_count": manifest.record_count,
                "failed_query_count": manifest.failed_query_count,
                "source_counts": dict(manifest.source_counts),
            },
            ensure_ascii=False,
            indent=2,
        )
    )
    return 0


def _parse_args(argv: list[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build controlled external reference cache for pcSecPichia.")
    parser.add_argument("--queries-jsonl", default="", help="Optional JSONL of ExternalReferenceQuery.to_dict() rows.")
    parser.add_argument("--homology-cache", default="", help="Optional homology cache JSONL path.")
    parser.add_argument("--name-audit-jsonl", default="", help="Optional name audit JSONL path.")
    parser.add_argument("--sources", default="uniprot,sgd,ncbi")
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--output-dir", default="local_runs/pichia_external_reference_cache/manual")
    parser.add_argument("--timeout", type=float, default=20.0)
    parser.add_argument("--retry-attempts", type=int, default=3)
    parser.add_argument("--min-interval", type=float, default=0.2)
    parser.add_argument("--no-catalog-default", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args(argv)


def _collect_queries(args: argparse.Namespace) -> tuple[ExternalReferenceQuery, ...]:
    queries: list[ExternalReferenceQuery] = []
    if args.queries_jsonl:
        queries.extend(_load_query_jsonl(_resolve_path(args.queries_jsonl)))
    homology_cache = _resolve_path(args.homology_cache) if args.homology_cache else None
    name_audit = _resolve_path(args.name_audit_jsonl) if args.name_audit_jsonl else None
    if homology_cache is not None or name_audit is not None:
        queries.extend(
            build_external_reference_queries(
                homology_cache=homology_cache,
                name_audit=name_audit,
            )
        )
    if not queries and not args.no_catalog_default:
        queries.extend(build_external_reference_queries(gene_catalog_rows=SECRETION_GENE_CATALOG))
    return tuple(queries)


def _load_query_jsonl(path: Path) -> tuple[ExternalReferenceQuery, ...]:
    rows: list[ExternalReferenceQuery] = []
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            payload = json.loads(line)
            if not isinstance(payload, dict):
                raise ValueError(f"Query JSONL row must be an object at {path}:{line_number}.")
            rows.append(
                ExternalReferenceQuery(
                    query_type=payload["query_type"],
                    query_value=payload["query_value"],
                    source_context=payload.get("source_context", "queries_jsonl"),
                    source_id=payload.get("source_id", "queries_jsonl"),
                    source_row_id=payload.get("source_row_id", ""),
                    preferred_sources=tuple(payload.get("preferred_sources") or ()),
                    warnings=tuple(payload.get("warnings") or ()),
                    metadata=dict(payload.get("metadata") or {}),
                )
            )
    return tuple(rows)


def _parse_sources(value: str) -> tuple[str, ...]:
    sources = tuple(source.strip().lower() for source in value.split(",") if source.strip())
    return sources or ("uniprot", "sgd", "ncbi")


def _resolve_output_dir(value: str) -> Path:
    return _resolve_path(value)


def _resolve_path(value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else REPO_ROOT / path


if __name__ == "__main__":
    raise SystemExit(main())
