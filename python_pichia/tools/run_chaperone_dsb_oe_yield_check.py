"""Direction-3-as-acceptance-dimension check (not a new screening product).

Question: does overexpressing the ER folding/chaperone (KAR2/BiP, SSA1,
YDJ1, BIP/NEFS) and disulfide-bond (PDI1, ERO1, ERV2 and their complexes)
catalog reactions already show a predicted secretion benefit for hLF/OPN
with the *current, live* default settings (misfolding constraint off) -
i.e. is this classic Pichia secretion-engineering strategy (co-expressing
folding helpers) already something the existing OE screen would surface,
without needing any new development?

Deliberately restricted to the CAT_ER_FOLDING/CAT_DSB subset of the
project's own curated catalog (services/gene_catalog.py) rather than a
genome-wide re-screen. Mirrors
run_erad_misfolding_flag_sensitivity_check.py's reaction-level tradeoff
helpers; also reports the misfolding-flag-on state for reference, since
folding/DSB capacity is conceptually adjacent to that constraint even
though it is not gated by it.

Usage (from python_pichia/, with src/ on PYTHONPATH):
    python tools/run_chaperone_dsb_oe_yield_check.py --target-id hLF --output-root ../local_runs/chaperone_dsb_oe_check
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
    wildtype_secretion_by_mu,
)
from pcsec_pichia.services.gene_catalog import CAT_DSB, CAT_ER_FOLDING  # noqa: E402
from pcsec_pichia.targets import load_builtin_targets  # noqa: E402


def _folding_dsb_reaction_candidates() -> tuple[dict[str, Any], ...]:
    return tuple(
        candidate
        for candidate in catalog_reaction_candidates()
        if candidate["category"] in (CAT_ER_FOLDING, CAT_DSB)
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
    for candidate in _folding_dsb_reaction_candidates():
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

    rows_by_flag: dict[bool, list[dict[str, Any]]] = {}
    for flag_value in (False, True):
        t0 = time.time()
        rows_by_flag[flag_value] = _run_reaction_level_rows(
            inputs,
            target,
            mode=args.mode,
            reference_growth_rate=DEFAULT_REFERENCE_GROWTH_RATE,
            factor=DEFAULT_OE_FACTOR,
            write_misfolding_constraints=flag_value,
        )
        print(f"[{time.strftime('%H:%M:%S')}] flag={flag_value}: {len(rows_by_flag[flag_value])} rows in {time.time() - t0:.1f}s")

    off_by_key = {str(row["gene_id"]): row for row in rows_by_flag[False]}
    on_by_key = {str(row["gene_id"]): row for row in rows_by_flag[True]}
    comparison = []
    for key in sorted(set(off_by_key) | set(on_by_key)):
        off_row = off_by_key.get(key)
        on_row = on_by_key.get(key)
        comparison.append(
            {
                "reaction_id": key,
                "common_name": (off_row or on_row or {}).get("common_name"),
                "category": (off_row or on_row or {}).get("category"),
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
    out_path = Path(args.output_root) / f"chaperone_dsb_oe_yield_check_{args.target_id}.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text, encoding="utf-8")
    print(f"[{time.strftime('%H:%M:%S')}] wrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
