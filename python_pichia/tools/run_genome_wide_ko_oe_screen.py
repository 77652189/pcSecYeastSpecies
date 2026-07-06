"""Standalone entry point for the genome-wide KO/OE growth-secretion tradeoff screen.

Usage (from the python_pichia/ directory, with src/ on PYTHONPATH):

    python tools/run_genome_wide_ko_oe_screen.py --mode fast --limit 50 --run-name pilot_50
    python tools/run_genome_wide_ko_oe_screen.py --mode fast --run-name overnight_full

Design rationale: ../docs/pichia_ko_oe_genome_screen_design.md (in the repo root docs/ folder).
Writes CSV + a Markdown summary to local_runs/<run-name>/.
"""

from __future__ import annotations

import argparse
import csv
import sys
import time
from pathlib import Path

_SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if str(_SRC_DIR) not in sys.path:
    sys.path.insert(0, str(_SRC_DIR))

from pcsec_pichia.loading import load_pcsec_pichia_inputs, repo_root  # noqa: E402
from pcsec_pichia.screens.genome_wide_tradeoff import (  # noqa: E402
    DEFAULT_REFERENCE_GROWTH_RATE,
    run_genome_wide_tradeoff_screen,
)
from pcsec_pichia.targets import load_builtin_targets  # noqa: E402

CSV_FIELDS = [
    "target_id",
    "gene_id",
    "intervention_type",
    "support_status",
    "secretory_process",
    "gpr_role",
    "mapping_confidence",
    "max_feasible_mu",
    "secretion_at_max_feasible_mu",
    "wildtype_max_feasible_mu",
    "wildtype_secretion_at_max_feasible_mu",
    "growth_retention_ratio",
    "secretion_ratio_vs_wildtype",
    "affected_reactions",
    "skipped_reason",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--targets", default="hLF,OPN_ALPHA_FULL_PROJECT", help="Comma-separated builtin target_ids.")
    parser.add_argument("--mode", choices=["fast", "precise"], default="fast")
    parser.add_argument("--reference-growth-rate", type=float, default=DEFAULT_REFERENCE_GROWTH_RATE)
    parser.add_argument("--limit", type=int, default=None, help="Only screen the first N genes (for pilot runs).")
    parser.add_argument("--genes", default=None, help="Comma-separated explicit gene_id list; overrides --limit.")
    parser.add_argument("--run-name", required=True, help="Output subfolder under local_runs/.")
    args = parser.parse_args()
    return args


def main() -> None:
    args = parse_args()
    root = repo_root()
    out_dir = root / "local_runs" / args.run_name
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[{time.strftime('%H:%M:%S')}] loading model and inputs...")
    inputs = load_pcsec_pichia_inputs(root)
    model = inputs.prepared_model
    all_gene_ids = [str(gene_id) for gene_id in model.genes]
    print(f"[{time.strftime('%H:%M:%S')}] model loaded: {len(all_gene_ids)} genes")

    if args.genes:
        gene_ids = [gene.strip() for gene in args.genes.split(",") if gene.strip()]
    elif args.limit is not None:
        gene_ids = all_gene_ids[: args.limit]
    else:
        gene_ids = all_gene_ids
    print(f"[{time.strftime('%H:%M:%S')}] screening {len(gene_ids)} genes, mode={args.mode}")

    targets_by_id = {target.target_id: target for target in load_builtin_targets(root)}
    target_ids = [target_id.strip() for target_id in args.targets.split(",") if target_id.strip()]

    csv_path = out_dir / "gene_tradeoff_rows.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=CSV_FIELDS)
        writer.writeheader()

        summary_lines: list[str] = [f"# Genome-wide KO/OE tradeoff screen: {args.run_name}", ""]
        summary_lines.append(f"- mode: {args.mode}")
        summary_lines.append(f"- reference growth rate: {args.reference_growth_rate}")
        summary_lines.append(f"- genes screened: {len(gene_ids)}")
        summary_lines.append(f"- targets: {', '.join(target_ids)}")
        summary_lines.append("")

        for target_id in target_ids:
            target = targets_by_id.get(target_id)
            if target is None:
                print(f"[WARN] unknown target_id {target_id!r}, skipping. Known: {sorted(targets_by_id)}")
                continue

            t0 = time.time()
            last_report = t0

            def progress_callback(done: int, total: int, gene_id: str, _t0: float = t0) -> None:
                nonlocal last_report
                now = time.time()
                if now - last_report >= 30 or done == total:
                    elapsed = now - _t0
                    rate = done / elapsed if elapsed > 0 else 0.0
                    eta_seconds = (total - done) / rate if rate > 0 else float("inf")
                    print(
                        f"[{time.strftime('%H:%M:%S')}] {target_id}: {done}/{total} genes "
                        f"({elapsed/60:.1f} min elapsed, ~{eta_seconds/60:.1f} min remaining), last={gene_id}"
                    )
                    last_report = now

            result = run_genome_wide_tradeoff_screen(
                model,
                target,
                inputs.amino_acids,
                inputs.metabolic,
                inputs.secretory,
                inputs.combined,
                gene_ids,
                mode=args.mode,
                reference_growth_rate=args.reference_growth_rate,
                progress_callback=progress_callback,
            )
            elapsed_total = time.time() - t0
            print(
                f"[{time.strftime('%H:%M:%S')}] {target_id} done in {elapsed_total/60:.1f} min "
                f"({len(result['rows'])} rows)"
            )

            for row in result["rows"]:
                writer.writerow(
                    {
                        "target_id": row["target_id"],
                        "gene_id": row["gene_id"],
                        "intervention_type": row["intervention_type"],
                        "support_status": row["support_status"],
                        "secretory_process": row["secretory_process"],
                        "gpr_role": row["gpr_role"],
                        "mapping_confidence": row["mapping_confidence"],
                        "max_feasible_mu": row["max_feasible_mu"],
                        "secretion_at_max_feasible_mu": row["secretion_at_max_feasible_mu"],
                        "wildtype_max_feasible_mu": row["wildtype_max_feasible_mu"],
                        "wildtype_secretion_at_max_feasible_mu": row["wildtype_secretion_at_max_feasible_mu"],
                        "growth_retention_ratio": row["growth_retention_ratio"],
                        "secretion_ratio_vs_wildtype": row["secretion_ratio_vs_wildtype"],
                        "affected_reactions": ";".join(row["affected_reactions"]),
                        "skipped_reason": row["skipped_reason"],
                    }
                )
            csv_file.flush()

            skipped = sum(1 for row in result["rows"] if row["skipped_reason"])
            summary_lines.append(f"## {target_id}")
            summary_lines.append(f"- wall time: {elapsed_total/60:.1f} minutes")
            summary_lines.append(f"- mu points: {result['mu_points']}")
            summary_lines.append(f"- wildtype max feasible mu: {result['wildtype_max_feasible_mu']}")
            summary_lines.append(f"- rows: {len(result['rows'])} (skipped/no structural effect: {skipped})")
            summary_lines.append("")

    report_path = out_dir / "SUMMARY.md"
    report_path.write_text("\n".join(summary_lines), encoding="utf-8")
    print(f"[{time.strftime('%H:%M:%S')}] wrote {csv_path}")
    print(f"[{time.strftime('%H:%M:%S')}] wrote {report_path}")


if __name__ == "__main__":
    main()
