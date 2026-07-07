from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict
from datetime import datetime
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "python_pichia" / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "python_pichia" / "src"))

from pcsec_pichia.homology.blast_runner import find_blastp_executable, make_blast_db, parse_blast_tsv, run_blastp
from pcsec_pichia.homology.cache_schema import BlastConfig, CatalogHomologyQuery, ProteinRecord
from pcsec_pichia.homology.catalog_inputs import secretion_catalog_sce_queries
from pcsec_pichia.homology.crosswalk import (
    build_homology_crosswalk,
    build_name_audit_rows,
    build_rule_transfer_audit_rows,
    load_external_name_reference_cache,
    summarize_homology_audits,
    write_homology_cache,
    write_name_audit_cache,
    write_rule_transfer_audit_cache,
)
from pcsec_pichia.homology.rbh import compute_reciprocal_best_hits
from pcsec_pichia.homology.sequence_sources import (
    load_pichia_model_gene_index,
    load_protein_sequences_from_mat,
    write_fasta,
)


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    root = Path(args.root).resolve()
    output_dir = Path(args.output_dir)
    if not output_dir.is_absolute():
        output_dir = root / output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    config = BlastConfig(
        max_evalue=args.max_evalue,
        min_identity=args.min_identity,
        min_coverage=args.min_coverage,
        max_target_seqs=args.max_target_seqs,
        threads=args.threads,
        blast_bin=Path(args.blast_bin) if args.blast_bin else None,
    )

    sce_records = _with_sce_annotations(
        load_protein_sequences_from_mat(root / "Data" / "pcSecYeast" / "Protein_Sequence.mat", "sce"),
        root,
    )
    pichia_records = _with_pichia_annotations(
        load_protein_sequences_from_mat(
            root / "Data" / "pcSecPichia" / "Protein_Sequence.mat",
            "pichia",
        ),
        root,
    )
    model_gene_index = load_pichia_model_gene_index(root)
    queries = secretion_catalog_sce_queries() if args.catalog_only else _all_sce_queries(sce_records)

    fasta_dir = output_dir / "fasta"
    sce_query_fasta = write_fasta(_records_for_queries(sce_records, queries), fasta_dir / "sce_catalog_queries.fasta")
    sce_proteome_fasta = write_fasta(sce_records, fasta_dir / "sce_proteome.fasta")
    pichia_fasta = write_fasta(pichia_records, fasta_dir / "pichia_proteome.fasta")
    forward_tsv = output_dir / "blast_forward.tsv"
    reverse_tsv = output_dir / "blast_reverse.tsv"

    blastp = find_blastp_executable(config)
    blast_status = "unavailable"
    if args.parse_existing:
        blast_status = "parsed_existing"
    elif args.skip_blast or blastp is None:
        blast_status = "skipped_unavailable" if blastp is None else "skipped_by_request"
        forward_tsv.write_text("", encoding="ascii")
        reverse_tsv.write_text("", encoding="ascii")
    else:
        sce_db = output_dir / "blastdb" / "sce"
        pichia_db = output_dir / "blastdb" / "pichia"
        db_results = [
            make_blast_db(sce_proteome_fasta, sce_db, blastp),
            make_blast_db(pichia_fasta, pichia_db, blastp),
        ]
        if not all(result.available for result in db_results):
            blast_status = "unavailable"
            forward_tsv.write_text("", encoding="ascii")
            reverse_tsv.write_text("", encoding="ascii")
        else:
            forward = run_blastp(sce_query_fasta, pichia_db, forward_tsv, config)
            reverse = run_blastp(pichia_fasta, sce_db, reverse_tsv, config)
            if not (forward.available and reverse.available):
                blast_status = "failed"
                return 2
            blast_status = "completed"

    forward_hits = parse_blast_tsv(forward_tsv) if forward_tsv.exists() else ()
    reverse_hits = parse_blast_tsv(reverse_tsv) if reverse_tsv.exists() else ()
    rbh_calls = compute_reciprocal_best_hits(forward_hits, reverse_hits)
    crosswalk = build_homology_crosswalk(queries, sce_records, pichia_records, model_gene_index, rbh_calls, config)
    external_references = _load_external_references(args.external_name_reference_cache, root)
    name_audit = build_name_audit_rows(crosswalk, external_references=external_references)
    rule_transfer_audit = build_rule_transfer_audit_rows(crosswalk)
    jsonl_path = output_dir / "sce_to_pichia_homology_cache.jsonl"
    tsv_path = output_dir / "sce_to_pichia_homology_cache.tsv"
    write_result = write_homology_cache(crosswalk, jsonl_path, tsv_path)
    name_jsonl_path = output_dir / "sce_to_pichia_name_audit.jsonl"
    name_tsv_path = output_dir / "sce_to_pichia_name_audit.tsv"
    name_write_result = write_name_audit_cache(name_audit, name_jsonl_path, name_tsv_path)
    rule_jsonl_path = output_dir / "sce_to_pichia_rule_transfer_audit.jsonl"
    rule_tsv_path = output_dir / "sce_to_pichia_rule_transfer_audit.tsv"
    rule_write_result = write_rule_transfer_audit_cache(rule_transfer_audit, rule_jsonl_path, rule_tsv_path)
    summary = summarize_homology_audits(
        blast_status=blast_status,
        homology_rows=crosswalk,
        name_audit_rows=name_audit,
        rule_transfer_rows=rule_transfer_audit,
    )
    summary_json_path = output_dir / "homology_audit_summary.json"
    summary_md_path = output_dir / "homology_audit_summary.md"
    _write_summary(summary_json_path, summary_md_path, summary)
    print(
        json.dumps(
            {
                "blast_status": blast_status,
                "row_count": write_result.row_count,
                "jsonl_path": str(jsonl_path),
                "tsv_path": str(tsv_path),
                "name_audit_row_count": name_write_result.row_count,
                "name_audit_jsonl_path": str(name_jsonl_path),
                "name_audit_tsv_path": str(name_tsv_path),
                "external_name_reference_count": len(external_references),
                "rule_transfer_row_count": rule_write_result.row_count,
                "rule_transfer_jsonl_path": str(rule_jsonl_path),
                "rule_transfer_tsv_path": str(rule_tsv_path),
                "summary_json_path": str(summary_json_path),
                "summary_md_path": str(summary_md_path),
            },
            ensure_ascii=False,
            indent=2,
        )
    )
    return 0


def _parse_args(argv: list[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build offline SCE-to-Pichia homology/name-audit cache.")
    parser.add_argument("--root", default=str(REPO_ROOT))
    parser.add_argument("--output-dir", default=str(Path("local_runs") / "pichia_homology_cache" / _run_name()))
    parser.add_argument("--blast-bin", default="")
    parser.add_argument("--min-identity", type=float, default=30.0)
    parser.add_argument("--min-coverage", type=float, default=50.0)
    parser.add_argument("--max-evalue", type=float, default=1e-10)
    parser.add_argument("--max-target-seqs", type=int, default=5)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--catalog-only", action="store_true")
    parser.add_argument("--skip-blast", action="store_true")
    parser.add_argument("--parse-existing", action="store_true")
    parser.add_argument(
        "--external-name-reference-cache",
        default="",
        help="Optional offline JSONL/TSV cache for future UniProt/NCBI/SGD/KEGG name crosschecks.",
    )
    return parser.parse_args(argv)


def _run_name() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def _with_sce_annotations(records: tuple[ProteinRecord, ...], root: Path) -> tuple[ProteinRecord, ...]:
    try:
        import pandas as pd
    except Exception:
        return records
    path = root / "Data" / "pcSecYeast" / "Protein_annotation_uniprot.xlsx"
    if not path.exists():
        return records
    annotation = pd.read_excel(path)
    by_orf: dict[str, tuple[str, tuple[str, ...]]] = {}
    for _, row in annotation.iterrows():
        orf = str(row.get("Gene names  (ordered locus )", "")).strip()
        if not orf:
            continue
        symbol = str(row.get("Entry name", "")).strip().upper()
        protein_names = str(row.get("Protein names", "")).upper()
        aliases = tuple(alias for alias in (symbol, "DOA10" if "DOA10" in protein_names else "") if alias)
        by_orf[orf] = (symbol, aliases)
    return tuple(
        ProteinRecord(
            organism=record.organism,
            gene_id=record.gene_id,
            sequence=record.sequence,
            symbol=by_orf.get(record.gene_id, ("", ()))[0] or record.symbol,
            aliases=by_orf.get(record.gene_id, ("", ()))[1] or record.aliases,
            accession=record.accession,
            source=record.source,
        )
        for record in records
    )


def _with_pichia_annotations(records: tuple[ProteinRecord, ...], root: Path) -> tuple[ProteinRecord, ...]:
    try:
        import pandas as pd
    except Exception:
        return records
    path = root / "Data" / "pcSecPichia" / "Protein_annotation_uniprot.xlsx"
    if not path.exists():
        return records
    annotation = pd.read_excel(path)
    by_locus: dict[str, tuple[str, tuple[str, ...], str]] = {}
    for _, row in annotation.iterrows():
        locus = str(row.get("Gene Names (ordered locus)", "")).strip()
        if not locus:
            continue
        entry_name = str(row.get("Entry Name", "")).strip()
        accession = str(row.get("Entry", "")).strip()
        aliases = tuple(alias for alias in (entry_name, accession) if alias and alias.lower() != "nan")
        by_locus[locus] = (entry_name if entry_name.lower() != "nan" else "", aliases, accession)
    return tuple(
        ProteinRecord(
            organism=record.organism,
            gene_id=record.gene_id,
            sequence=record.sequence,
            symbol=by_locus.get(record.gene_id, ("", (), ""))[0] or record.symbol,
            aliases=by_locus.get(record.gene_id, ("", (), ""))[1] or record.aliases,
            accession=by_locus.get(record.gene_id, ("", (), ""))[2] or record.accession,
            source=record.source,
        )
        for record in records
    )


def _records_for_queries(
    records: tuple[ProteinRecord, ...],
    queries: tuple[CatalogHomologyQuery, ...],
) -> tuple[ProteinRecord, ...]:
    lookup: dict[str, ProteinRecord] = {}
    for record in records:
        keys = [record.gene_id]
        if record.symbol:
            keys.append(record.symbol)
        keys.extend(record.aliases)
        for key in keys:
            lookup.setdefault(key.upper(), record)
    selected: dict[str, ProteinRecord] = {}
    for query in queries:
        for key in (query.query_symbol, *query.aliases):
            record = lookup.get(key.upper())
            if record:
                selected[record.gene_id] = record
                break
    return tuple(selected.values())


def _all_sce_queries(records: tuple[ProteinRecord, ...]) -> tuple[CatalogHomologyQuery, ...]:
    return tuple(
        CatalogHomologyQuery(
            internal_common_name=record.symbol or record.gene_id,
            query_symbol=record.symbol or record.gene_id,
            aliases=(record.gene_id, *record.aliases),
            source="sce_proteome",
        )
        for record in records
    )


def _load_external_references(cache_path: str, root: Path):
    if not cache_path:
        return ()
    path = Path(cache_path)
    if not path.is_absolute():
        path = root / path
    if not path.exists():
        return ()
    return load_external_name_reference_cache(path)


def _write_summary(json_path: Path, md_path: Path, summary: object) -> None:
    payload = asdict(summary)
    json_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    lines = [
        "# Pichia Homology Audit Summary",
        "",
        f"- blast_status: {payload['blast_status']}",
        f"- homology_row_count: {payload['homology_row_count']}",
        f"- name_audit_row_count: {payload['name_audit_row_count']}",
        f"- rule_transfer_row_count: {payload['rule_transfer_row_count']}",
        "",
        "## Homology Review Status Counts",
        "",
    ]
    for status, count in payload["homology_review_status_counts"].items():
        lines.append(f"- {status}: {count}")
    lines.extend(["", "## Name Consistency Status Counts", ""])
    for status, count in payload["name_consistency_status_counts"].items():
        lines.append(f"- {status}: {count}")
    lines.extend(["", "## Rule Transfer Status Counts", ""])
    for status, count in payload["rule_transfer_status_counts"].items():
        lines.append(f"- {status}: {count}")
    lines.extend(
        [
            "",
            "## Boundary",
            "",
            "- BLAST/RBH is sequence-level homology evidence only.",
            "- RBH does not automatically update SECRETION_GENE_CATALOG.",
            "- RBH does not automatically make a Pichia candidate KO/OE-operable.",
            "- Model operability is reported separately via in_model_gene_index.",
        ]
    )
    md_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
