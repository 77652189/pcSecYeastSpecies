"""Parallel (multiprocessing) entry point for the genome-wide KO/OE tradeoff screen.

Each gene's KO/OE solve is independent, so this splits (target, gene) pairs
across worker processes. Each worker loads its own copy of the model once and
caches the per-target preparation (target-specific enzyme data, baseline
wildtype secretion by mu) so it is only computed once per worker per target,
not once per gene.

Usage (from python_pichia/):

    python tools/run_genome_wide_ko_oe_screen_parallel.py --mode fast --workers 6 --run-name overnight_full

Design rationale: ../docs/pichia_ko_oe_genome_screen_design.md
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path
from typing import Any

_SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if str(_SRC_DIR) not in sys.path:
    sys.path.insert(0, str(_SRC_DIR))

from pcsec_pichia.loading import load_pcsec_pichia_inputs, repo_root  # noqa: E402
from pcsec_pichia.screens.genome_wide_tradeoff import (  # noqa: E402
    COMPLEX_OE_HYPOTHESIS_ASSUMPTION,
    DEFAULT_REFERENCE_GROWTH_RATE,
    catalog_reaction_candidates,
    gene_ko_tradeoff,
    gene_oe_tradeoff,
    mu_points_for_mode,
    reaction_ko_tradeoff,
    reaction_oe_tradeoff,
    resolve_complex_subunit_oe_hypothesis_candidates,
    wildtype_secretion_by_mu,
)
from pcsec_pichia.screens.gene_interventions import plan_gene_knockout, plan_gene_overexpression  # noqa: E402
from pcsec_pichia.screens._prototype_adapter import build_supported_target_model, build_target_enzymedata  # noqa: E402
from pcsec_pichia.targets import load_builtin_targets  # noqa: E402

CSV_FIELDS = [
    "target_id",
    "gene_id",
    "common_name",
    "candidate_kind",
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
    "solve_outcome_counts",
    "has_timeout",
    "timeout_mu_points",
    "proven_infeasible_mu_points",
    "other_solver_failure_mu_points",
    "feasibility_interpretation",
    "affected_reactions",
    "skipped_reason",
    "hypothesis_note",
]

_WORKER: dict[str, Any] = {}


def _worker_init(repo_root_str: str) -> None:
    """Runs once per worker process: load the model a single time."""
    root = Path(repo_root_str)
    inputs = load_pcsec_pichia_inputs(root)
    _WORKER["root"] = root
    _WORKER["inputs"] = inputs
    _WORKER["targets_by_id"] = {target.target_id: target for target in load_builtin_targets(root)}
    _WORKER["target_cache"] = {}


def _prepare_target(target_id: str, mode: str, reference_growth_rate: float) -> dict[str, Any]:
    """Cached per worker per (target, mode, reference_growth_rate); avoids rebuilding per gene."""
    key = (target_id, mode, reference_growth_rate)
    cache = _WORKER["target_cache"]
    if key in cache:
        return cache[key]

    inputs = _WORKER["inputs"]
    target = _WORKER["targets_by_id"][target_id]
    build = build_supported_target_model(inputs.prepared_model, target, inputs.amino_acids)
    if not build.supported or build.model is None or build.exchange_reaction_id is None:
        raise ValueError(f"Target {target_id!r} unsupported (build_status={build.build_status!r})")

    target_enzymedata = build_target_enzymedata(target, build.model, inputs.secretory)
    target_secretory = inputs.secretory.with_reaction_coefficients(target_enzymedata.reaction_coefficients)
    target_combined = inputs.combined.with_target(target_enzymedata)
    complex_subunits = getattr(inputs.secretory, "complex_subunits", None)
    mu_points = mu_points_for_mode(reference_growth_rate, mode)
    baseline_by_mu = wildtype_secretion_by_mu(
        build.model, build.exchange_reaction_id, inputs.metabolic, target_secretory, target_combined, mu_points
    )
    wt_feasible = [
        {"mu": mu, "success": entry["success"], "secretion_flux": entry["objective_value"]}
        for mu, entry in baseline_by_mu.items()
        if entry["success"]
    ]
    wildtype_best = max(wt_feasible, key=lambda point: point["mu"]) if wt_feasible else None

    prepared = {
        "model": build.model,
        "exchange_reaction_id": build.exchange_reaction_id,
        "metabolic": inputs.metabolic,
        "target_secretory": target_secretory,
        "target_combined": target_combined,
        "complex_subunits": complex_subunits,
        "mu_points": mu_points,
        "baseline_by_mu": baseline_by_mu,
        "wildtype_best": wildtype_best,
    }
    cache[key] = prepared
    return prepared


def _attach_wildtype(row: dict[str, Any], wildtype_best: dict[str, Any] | None) -> dict[str, Any]:
    wt_mu = wildtype_best["mu"] if wildtype_best else None
    wt_secretion = wildtype_best["secretion_flux"] if wildtype_best else None
    row["wildtype_max_feasible_mu"] = wt_mu
    row["wildtype_secretion_at_max_feasible_mu"] = wt_secretion
    row["growth_retention_ratio"] = row["max_feasible_mu"] / wt_mu if row["max_feasible_mu"] is not None and wt_mu else None
    row["secretion_ratio_vs_wildtype"] = (
        row["secretion_at_max_feasible_mu"] / wt_secretion
        if row["secretion_at_max_feasible_mu"] is not None and wt_secretion
        else None
    )
    return row


def _csv_value(value: object) -> object:
    if isinstance(value, (dict, list, tuple)):
        return json.dumps(value, ensure_ascii=False)
    return value


def _row_to_csv_record(row: dict[str, Any]) -> dict[str, object]:
    return {
        "target_id": row["target_id"],
        "gene_id": row["gene_id"],
        "common_name": row.get("common_name", ""),
        "candidate_kind": row.get("candidate_kind", "gene"),
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
        "solve_outcome_counts": _csv_value(row.get("solve_outcome_counts", {})),
        "has_timeout": row.get("has_timeout", False),
        "timeout_mu_points": _csv_value(row.get("timeout_mu_points", ())),
        "proven_infeasible_mu_points": _csv_value(row.get("proven_infeasible_mu_points", ())),
        "other_solver_failure_mu_points": _csv_value(row.get("other_solver_failure_mu_points", ())),
        "feasibility_interpretation": row.get("feasibility_interpretation", ""),
        "affected_reactions": ";".join(row.get("affected_reactions", ())),
        "skipped_reason": row["skipped_reason"],
        "hypothesis_note": row.get("hypothesis_note", ""),
    }


def _run_one_gene(target_id: str, gene_id: str, mode: str, reference_growth_rate: float) -> list[dict[str, Any]]:
    """Task function executed in a worker process for one (target, gene) pair."""
    prepared = _prepare_target(target_id, mode, reference_growth_rate)
    model = prepared["model"]
    ko_plan = plan_gene_knockout(model, gene_id)
    oe_plan = plan_gene_overexpression(model, gene_id, complex_subunits=prepared["complex_subunits"])

    ko_row = gene_ko_tradeoff(
        model, gene_id, ko_plan, prepared["exchange_reaction_id"], prepared["metabolic"],
        prepared["target_secretory"], prepared["target_combined"], prepared["mu_points"], prepared["complex_subunits"],
    )
    oe_row = gene_oe_tradeoff(
        model, gene_id, oe_plan, prepared["exchange_reaction_id"], prepared["metabolic"],
        prepared["target_secretory"], prepared["target_combined"], prepared["mu_points"], prepared["baseline_by_mu"],
        prepared["complex_subunits"],
    )
    rows = []
    for row in (ko_row, oe_row):
        row["target_id"] = target_id
        row.setdefault("common_name", "")
        row.setdefault("candidate_kind", "gene")
        _attach_wildtype(row, prepared["wildtype_best"])
        rows.append(row)
    return rows


def _run_one_catalog_reaction(target_id: str, candidate: dict[str, Any], mode: str, reference_growth_rate: float) -> list[dict[str, Any]]:
    """Task function executed in a worker process for one (target, curated catalog reaction) pair.

    Unlike _run_one_gene, a catalog candidate is already intervention-type-specific
    (catalog_reaction_candidates() dedupes by (intervention_type, reaction_id)), so this
    produces one row per task instead of a KO+OE pair.
    """
    prepared = _prepare_target(target_id, mode, reference_growth_rate)
    model = prepared["model"]
    common_args = (
        candidate["reaction_id"], candidate["common_name"], candidate["category"],
        prepared["exchange_reaction_id"], prepared["metabolic"], prepared["target_secretory"],
        prepared["target_combined"], prepared["mu_points"],
    )
    if candidate["intervention_type"] == "KO":
        row = reaction_ko_tradeoff(model, *common_args, prepared["complex_subunits"])
    else:
        row = reaction_oe_tradeoff(model, *common_args, prepared["baseline_by_mu"], prepared["complex_subunits"])
    row["target_id"] = target_id
    _attach_wildtype(row, prepared["wildtype_best"])
    return [row]


def _run_one_complex_hypothesis_target(target_id: str, gene_ids: list[str], mode: str, reference_growth_rate: float) -> list[dict[str, Any]]:
    """Task function executed in a worker process for one target's complex-OE-hypothesis gene list.

    Unlike _run_one_gene/_run_one_catalog_reaction (one task per gene/reaction), this takes the
    whole per-target gene list at once: resolving genes to reactions needs a built model
    (resolve_complex_subunit_oe_hypothesis_candidates), unlike catalog_reaction_candidates()
    which is static and can be split into per-reaction tasks before any model exists. The
    hypothesis scope is minute-scale (a few dozen genes collapsing to a handful of reactions
    per target after dedup), so running one target's reactions sequentially within a task
    costs little and avoids a second model build just to pre-split the work.
    """
    prepared = _prepare_target(target_id, mode, reference_growth_rate)
    model = prepared["model"]
    candidates = resolve_complex_subunit_oe_hypothesis_candidates(model, gene_ids, prepared["complex_subunits"])
    rows = []
    for candidate in candidates:
        row = reaction_oe_tradeoff(
            model, candidate["reaction_id"], candidate["common_name"], candidate["category"],
            prepared["exchange_reaction_id"], prepared["metabolic"], prepared["target_secretory"],
            prepared["target_combined"], prepared["mu_points"], prepared["baseline_by_mu"], prepared["complex_subunits"],
            candidate_kind="complex_oe_hypothesis", hypothesis_note=COMPLEX_OE_HYPOTHESIS_ASSUMPTION,
        )
        row["target_id"] = target_id
        _attach_wildtype(row, prepared["wildtype_best"])
        rows.append(row)
    return rows


def _task_item_label(item: Any) -> str:
    if isinstance(item, str):
        return item
    if isinstance(item, dict):
        return item["reaction_id"]
    if isinstance(item, list):
        return f"{len(item)}_hypothesis_genes"
    return str(item)


def _write_status(status_path: Path, **fields: Any) -> None:
    payload = {"updated_at": datetime.now().isoformat(), **fields}
    tmp_path = status_path.with_suffix(".json.tmp")
    tmp_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    tmp_path.replace(status_path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--targets", default="hLF,OPN_ALPHA_FULL_PROJECT")
    parser.add_argument("--mode", choices=["fast", "precise"], default="fast")
    parser.add_argument("--reference-growth-rate", type=float, default=DEFAULT_REFERENCE_GROWTH_RATE)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--genes", default=None)
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--run-name", required=True)
    parser.add_argument(
        "--scope",
        choices=["gene", "catalog", "complex_hypothesis"],
        default="gene",
        help="gene=all ~1025 model genes (default, hour-scale); catalog=~30 unique reactions from "
        "the curated SECRETION_GENE_CATALOG literature shortlist (minute-scale); "
        "complex_hypothesis=OE-only hypothetical-whole-complex-overexpression test for "
        "complex_subunit genes whose KO already showed a real secretion decrease in a prior "
        "gene-scope run (minute-scale; requires --source-run). See COMPLEX_OE_HYPOTHESIS_ASSUMPTION "
        "in genome_wide_tradeoff.py for what this scope does and does not prove.",
    )
    parser.add_argument(
        "--source-run",
        default=None,
        help="Run name under local_runs/ whose gene_tradeoff_rows.csv supplies the qualifying gene "
        "list for --scope complex_hypothesis, via genome_wide_screen_analysis."
        "complex_subunit_oe_hypothesis_candidates() (per target). Required for that scope.",
    )
    return parser.parse_args()


def main(args: argparse.Namespace | None = None) -> None:
    args = args or parse_args()
    root = repo_root()
    out_dir = root / "local_runs" / args.run_name
    out_dir.mkdir(parents=True, exist_ok=True)
    status_path = out_dir / "status.json"
    errors: list[str] = []

    _write_status(status_path, status="starting", pid=os.getpid(), done=0, total=0, message="loading model...")
    target_ids = [target_id.strip() for target_id in args.targets.split(",") if target_id.strip()]

    if args.scope == "catalog":
        print(f"[{time.strftime('%H:%M:%S')}] loading model and curated catalog candidates...")
        load_pcsec_pichia_inputs(root)  # warms up any shared caches before forking workers
        catalog_candidates = catalog_reaction_candidates()
        tasks = [(target_id, candidate) for target_id in target_ids for candidate in catalog_candidates]
        total = len(tasks)
        print(
            f"[{time.strftime('%H:%M:%S')}] {len(catalog_candidates)} curated reactions x {len(target_ids)} targets "
            f"= {total} tasks, mode={args.mode}, workers={args.workers}"
        )
    elif args.scope == "complex_hypothesis":
        if not args.source_run:
            raise SystemExit("--scope complex_hypothesis requires --source-run <run_name> (a completed gene-scope run)")
        app_root = str(root)
        if app_root not in sys.path:
            sys.path.insert(0, app_root)
        from app.services.genome_wide_screen_analysis import (
            complex_subunit_oe_hypothesis_candidates,
            load_gene_tradeoff_csv,
        )

        source_csv = root / "local_runs" / args.source_run / "gene_tradeoff_rows.csv"
        print(f"[{time.strftime('%H:%M:%S')}] deriving hypothesis candidate genes from {source_csv}...")
        source_frame = load_gene_tradeoff_csv(str(source_csv))
        load_pcsec_pichia_inputs(root)  # warms up any shared caches before forking workers

        gene_ids_by_target = {
            target_id: complex_subunit_oe_hypothesis_candidates(source_frame, target_id) for target_id in target_ids
        }
        for target_id, gene_ids in gene_ids_by_target.items():
            print(f"[{time.strftime('%H:%M:%S')}]   {target_id}: {len(gene_ids)} hypothesis candidate genes")
        tasks = [(target_id, gene_ids) for target_id, gene_ids in gene_ids_by_target.items() if gene_ids]
        total = len(tasks)
        print(
            f"[{time.strftime('%H:%M:%S')}] {total} target(s) with hypothesis candidates, "
            f"mode={args.mode}, workers={args.workers}"
        )
    else:
        print(f"[{time.strftime('%H:%M:%S')}] loading model to discover gene list...")
        inputs = load_pcsec_pichia_inputs(root)
        all_gene_ids = [str(gene_id) for gene_id in inputs.prepared_model.genes]

        if args.genes:
            gene_ids = [gene.strip() for gene in args.genes.split(",") if gene.strip()]
        elif args.limit is not None:
            gene_ids = all_gene_ids[: args.limit]
        else:
            gene_ids = all_gene_ids

        tasks = [(target_id, gene_id) for target_id in target_ids for gene_id in gene_ids]
        total = len(tasks)
        print(
            f"[{time.strftime('%H:%M:%S')}] {len(gene_ids)} genes x {len(target_ids)} targets = {total} tasks, "
            f"mode={args.mode}, workers={args.workers}"
        )

    csv_filenames = {
        "catalog": "catalog_reaction_tradeoff_rows.csv",
        "complex_hypothesis": "complex_hypothesis_tradeoff_rows.csv",
    }
    csv_path = out_dir / csv_filenames.get(args.scope, "gene_tradeoff_rows.csv")
    t0 = time.time()
    last_report = t0
    done = 0
    _write_status(
        status_path, status="running", pid=os.getpid(), done=0, total=total,
        started_at=datetime.now().isoformat(), targets=target_ids, mode=args.mode, scope=args.scope,
    )
    task_fns = {"catalog": _run_one_catalog_reaction, "complex_hypothesis": _run_one_complex_hypothesis_target}
    task_fn = task_fns.get(args.scope, _run_one_gene)
    with csv_path.open("w", newline="", encoding="utf-8") as csv_file, ProcessPoolExecutor(
        max_workers=args.workers, initializer=_worker_init, initargs=(str(root),)
    ) as pool:
        writer = csv.DictWriter(csv_file, fieldnames=CSV_FIELDS)
        writer.writeheader()

        futures = {
            pool.submit(task_fn, target_id, item, args.mode, args.reference_growth_rate): (
                target_id,
                _task_item_label(item),
            )
            for target_id, item in tasks
        }
        for future in as_completed(futures):
            target_id, item_label = futures[future]
            try:
                rows = future.result()
            except Exception as exc:  # noqa: BLE001 - record and continue; one bad candidate should not kill the run
                print(f"[ERROR] {target_id}/{item_label}: {exc!r}")
                errors.append(f"{target_id}/{item_label}: {exc!r}")
                done += 1
                continue
            for row in rows:
                writer.writerow(_row_to_csv_record(row))
            done += 1
            now = time.time()
            if now - last_report >= 60 or done == total:
                elapsed = now - t0
                rate = done / elapsed if elapsed > 0 else 0.0
                eta_seconds = (total - done) / rate if rate > 0 else float("inf")
                print(
                    f"[{time.strftime('%H:%M:%S')}] {done}/{total} tasks "
                    f"({elapsed/60:.1f} min elapsed, ~{eta_seconds/60:.1f} min remaining), last={target_id}/{item_label}"
                )
                _write_status(
                    status_path, status="running", pid=os.getpid(), done=done, total=total,
                    elapsed_minutes=round(elapsed / 60, 1),
                    eta_minutes=(round(eta_seconds / 60, 1) if eta_seconds != float("inf") else None),
                    last_gene=f"{target_id}/{item_label}", error_count=len(errors),
                )
                last_report = now
                csv_file.flush()

    total_elapsed = time.time() - t0
    print(f"[{time.strftime('%H:%M:%S')}] all done in {total_elapsed/60:.1f} minutes. wrote {csv_path}")
    _write_status(
        status_path, status="done", pid=os.getpid(), done=total, total=total,
        elapsed_minutes=round(total_elapsed / 60, 1), csv_path=str(csv_path), error_count=len(errors),
        errors=errors[:50], targets=target_ids, mode=args.mode, scope=args.scope,
    )


def _main_with_error_status() -> None:
    args = parse_args()
    status_path = repo_root() / "local_runs" / args.run_name / "status.json"
    try:
        main(args)
    except Exception as exc:  # noqa: BLE001 - surface fatal errors to the polling UI, then re-raise for the log
        status_path.parent.mkdir(parents=True, exist_ok=True)
        _write_status(status_path, status="error", pid=os.getpid(), message=f"{type(exc).__name__}: {exc}")
        raise


if __name__ == "__main__":
    _main_with_error_status()
