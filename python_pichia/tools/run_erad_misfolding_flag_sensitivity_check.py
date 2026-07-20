"""Direction-3-as-acceptance-dimension check (not a new screening product).

Question: does enabling the currently-disabled ERAD/misfolding constraint
(PichiaSimulationRequest.enable_misfolding_constraint /
genome_wide_tradeoff.write_misfolding_constraints, default False everywhere)
change candidate rankings for the ~10 genes/reactions the project's own
curated catalog already tags as ERAD/proteasome-related (CAT_ERAD,
CAT_PROTEASOME in services/gene_catalog.py)?

Deliberately restricted to this small, already-known candidate set rather
than a genome-wide re-screen: this answers "is it worth investing further
in enabling this constraint" cheaply, without building a new product.
Reuses run_genome_wide_tradeoff_screen (for the 2 candidates with a real
gene_id: PEP4, PRB1) and the same reaction_ko_tradeoff/reaction_oe_tradeoff
primitives run_catalog_reaction_tradeoff_screen itself uses (restricted to
the ERAD/proteasome subset of catalog_reaction_candidates(), instead of all
~30 catalog reactions).

Usage (from python_pichia/, with src/ on PYTHONPATH):
    python tools/run_erad_misfolding_flag_sensitivity_check.py --target-id hLF --output-root ../local_runs/erad_misfolding_sensitivity
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any

_SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if str(_SRC_DIR) not in sys.path:
    sys.path.insert(0, str(_SRC_DIR))

from pcsec_pichia.loading import load_pcsec_pichia_inputs, repo_root  # noqa: E402
from pcsec_pichia.screens import build_supported_target_model, build_target_enzymedata  # noqa: E402
from pcsec_pichia.screens.genome_wide_tradeoff import (  # noqa: E402
    DEFAULT_OE_FACTOR,
    DEFAULT_REFERENCE_GROWTH_RATE,
    catalog_reaction_candidates,
    mu_points_for_mode,
    reaction_ko_tradeoff,
    reaction_oe_tradeoff,
    run_genome_wide_tradeoff_screen,
    wildtype_secretion_by_mu,
)
from pcsec_pichia.services.gene_catalog import CAT_ERAD, CAT_PROTEASOME  # noqa: E402
from pcsec_pichia.targets import load_builtin_targets  # noqa: E402

# PEP4, PRB1: the only two CAT_PROTEASOME/CAT_ERAD entries with a real
# gene_id (vacuolar proteases, routinely knocked out in Pichia expression
# strains). Every other ERAD/proteasome catalog entry is a multi-subunit
# complex the model only supports at reaction level (see
# catalog_reaction_candidates()), handled separately below.
GENE_LEVEL_CANDIDATES: tuple[str, ...] = ("PAS_chr2-2_0107", "PAS_chr2-1_0785")


def _erad_proteasome_reaction_candidates() -> tuple[dict[str, Any], ...]:
    return tuple(
        candidate
        for candidate in catalog_reaction_candidates()
        if candidate["category"] in (CAT_ERAD, CAT_PROTEASOME)
    )


def _attach_wildtype_comparison(row: dict[str, Any], wildtype_best: dict[str, Any] | None) -> dict[str, Any]:
    wt_mu = wildtype_best["mu"] if wildtype_best else None
    wt_secretion = wildtype_best["secretion_flux"] if wildtype_best else None
    row["wildtype_max_feasible_mu"] = wt_mu
    row["wildtype_secretion_at_max_feasible_mu"] = wt_secretion
    row["secretion_ratio_vs_wildtype"] = (
        row["secretion_at_max_feasible_mu"] / wt_secretion
        if row["secretion_at_max_feasible_mu"] is not None and wt_secretion
        else None
    )
    return row


def _max_feasible_point(points: list[dict[str, Any]]) -> dict[str, Any] | None:
    feasible = [point for point in points if point["success"]]
    if not feasible:
        return None
    return max(feasible, key=lambda point: point["mu"])


def _run_reaction_level_rows(
    inputs: Any,
    target: Any,
    *,
    mode: str,
    reference_growth_rate: float,
    factor: float,
    write_misfolding_constraints: bool,
) -> list[dict[str, Any]]:
    build = build_supported_target_model(inputs.prepared_model, target, inputs.amino_acids)
    target_enzymedata = build_target_enzymedata(target, build.model, inputs.secretory)
    target_secretory = inputs.secretory.with_reaction_coefficients(target_enzymedata.reaction_coefficients)
    target_combined = inputs.combined.with_target(target_enzymedata)
    complex_subunits = getattr(inputs.secretory, "complex_subunits", None)
    mu_points = mu_points_for_mode(reference_growth_rate, mode)
    baseline_by_mu = wildtype_secretion_by_mu(
        build.model,
        build.exchange_reaction_id,
        inputs.metabolic,
        target_secretory,
        target_combined,
        mu_points,
        False,
        write_misfolding_constraints,
    )
    wildtype_best = _max_feasible_point(
        [{"mu": mu, "success": entry["success"], "secretion_flux": entry["objective_value"]} for mu, entry in baseline_by_mu.items()]
    )

    rows: list[dict[str, Any]] = []
    for candidate in _erad_proteasome_reaction_candidates():
        if candidate["intervention_type"] == "KO":
            row = reaction_ko_tradeoff(
                build.model,
                candidate["reaction_id"],
                candidate["common_name"],
                candidate["category"],
                build.exchange_reaction_id,
                inputs.metabolic,
                target_secretory,
                target_combined,
                mu_points,
                complex_subunits,
                False,
                write_misfolding_constraints,
            )
        else:
            row = reaction_oe_tradeoff(
                build.model,
                candidate["reaction_id"],
                candidate["common_name"],
                candidate["category"],
                build.exchange_reaction_id,
                inputs.metabolic,
                target_secretory,
                target_combined,
                mu_points,
                baseline_by_mu,
                complex_subunits,
                factor,
                False,
                write_misfolding_constraints,
            )
        row["target_id"] = target.target_id
        _attach_wildtype_comparison(row, wildtype_best)
        rows.append(row)
    return rows


def _run_for_flag_state(inputs: Any, target: Any, write_misfolding_constraints: bool, mode: str) -> list[dict[str, Any]]:
    gene_result = run_genome_wide_tradeoff_screen(
        inputs.prepared_model,
        target,
        inputs.amino_acids,
        inputs.metabolic,
        inputs.secretory,
        inputs.combined,
        list(GENE_LEVEL_CANDIDATES),
        mode=mode,
        reference_growth_rate=DEFAULT_REFERENCE_GROWTH_RATE,
        write_misfolding_constraints=write_misfolding_constraints,
    )
    reaction_rows = _run_reaction_level_rows(
        inputs,
        target,
        mode=mode,
        reference_growth_rate=DEFAULT_REFERENCE_GROWTH_RATE,
        factor=DEFAULT_OE_FACTOR,
        write_misfolding_constraints=write_misfolding_constraints,
    )
    return [*gene_result["rows"], *reaction_rows]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-id", default="hLF")
    parser.add_argument("--mode", choices=["fast", "precise"], default="fast")
    parser.add_argument("--output-root", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = repo_root()
    inputs = load_pcsec_pichia_inputs(root)
    targets_by_id = {target.target_id: target for target in load_builtin_targets(root)}
    target = targets_by_id[args.target_id]

    results_by_flag: dict[bool, list[dict[str, Any]]] = {}
    for flag_value in (False, True):
        t0 = time.time()
        results_by_flag[flag_value] = _run_for_flag_state(inputs, target, flag_value, args.mode)
        print(f"[{time.strftime('%H:%M:%S')}] flag={flag_value}: {len(results_by_flag[flag_value])} rows in {time.time() - t0:.1f}s")

    off_by_key = {(str(row["gene_id"]), str(row["intervention_type"])): row for row in results_by_flag[False]}
    on_by_key = {(str(row["gene_id"]), str(row["intervention_type"])): row for row in results_by_flag[True]}
    comparison = []
    for key in sorted(set(off_by_key) | set(on_by_key)):
        off_row = off_by_key.get(key)
        on_row = on_by_key.get(key)
        comparison.append(
            {
                "gene_or_reaction_id": key[0],
                "common_name": (off_row or on_row or {}).get("common_name"),
                "intervention_type": key[1],
                "secretion_ratio_vs_wildtype_flag_off": off_row.get("secretion_ratio_vs_wildtype") if off_row else None,
                "secretion_ratio_vs_wildtype_flag_on": on_row.get("secretion_ratio_vs_wildtype") if on_row else None,
                "max_feasible_mu_flag_off": off_row.get("max_feasible_mu") if off_row else None,
                "max_feasible_mu_flag_on": on_row.get("max_feasible_mu") if on_row else None,
                "skipped_reason_flag_off": off_row.get("skipped_reason") if off_row else "missing",
                "skipped_reason_flag_on": on_row.get("skipped_reason") if on_row else "missing",
            }
        )

    report = {"target_id": args.target_id, "mode": args.mode, "comparison": comparison}
    text = json.dumps(report, ensure_ascii=False, indent=2, default=str)
    print(text)
    out_path = Path(args.output_root) / f"erad_misfolding_flag_sensitivity_{args.target_id}.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text, encoding="utf-8")
    print(f"[{time.strftime('%H:%M:%S')}] wrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
