"""App service (C2 编排): run the modified-strain → 下一步 OE 候选 two-pass and assemble the readout.

Thin facade over the engine orchestration `pcsec_pichia.next_oe_candidates.analyze_next_oe_candidates`
(solve modified strain → bottleneck attribution → bounded OE dose-response) followed by the pure
app-layer ranking `per_strain_oe_candidates.build_next_oe_candidates_readout` (C1). Keeps the engine
call + reload guard in the service layer (mirrors `pichia_secretion_runner`) and never raises into the
UI: engine failures come back as a stable readout dict with `error` set and no candidates.
"""

from __future__ import annotations

from typing import Any, Sequence

from app.services.per_strain_oe_candidates import DEFAULT_TOP_N, build_next_oe_candidates_readout
from app.services.pichia_secretion_runner import _ensure_pcsec_pichia_analysis_api


def run_next_oe_candidate_analysis(
    *,
    target_id: str,
    ko_reaction_ids: Sequence[str] = (),
    oe_reaction_ids: Sequence[str] = (),
    oe_factor: float = 2.0,
    mu: float = 0.10,
    media_type: int = 4,
    carbon_source_id: str = "glucose",
    enable_ribosome_translation_constraint: bool = False,
    enable_misfolding_constraint: bool = False,
    target_input: Any | None = None,
    leader_candidate: Any | None = None,
    top_n: int = DEFAULT_TOP_N,
) -> dict[str, Any]:
    """Solve the modified strain and return the ranked next-OE-candidate readout (opt-in, extra solves).

    Returns C1's readout dict (`candidates` / `dose_response_available` / `caveats` / ...) enriched with
    the modified-solve status, applied modifications, and honest warnings. On any engine error, returns
    an empty-candidate readout with `error` set instead of raising.
    """

    try:
        _ensure_pcsec_pichia_analysis_api()
        from pcsec_pichia.next_oe_candidates import analyze_next_oe_candidates

        raw = analyze_next_oe_candidates(
            target_id=target_id,
            ko_reaction_ids=tuple(ko_reaction_ids),
            oe_reaction_ids=tuple(oe_reaction_ids),
            oe_factor=float(oe_factor),
            mu=float(mu),
            media_type=int(media_type),
            carbon_source_id=carbon_source_id,
            enable_ribosome_translation_constraint=bool(enable_ribosome_translation_constraint),
            enable_misfolding_constraint=bool(enable_misfolding_constraint),
            target_input=target_input,
            leader_candidate=leader_candidate,
            top_n=int(top_n),
        )
    except Exception as exc:  # keep the UI facade stable; engine detail lands in `error`.
        readout = build_next_oe_candidates_readout([], None, top_n=int(top_n))
        readout.update(
            {
                "target_id": target_id,
                "modified_solve_success": False,
                "error": f"{type(exc).__name__}: {exc}",
            }
        )
        return readout

    readout = build_next_oe_candidates_readout(
        raw.get("oe_actionable_bottlenecks") or [],
        raw.get("dose_response"),
        top_n=int(top_n),
    )
    readout.update(
        {
            "target_id": raw.get("target_id"),
            "modified_solve_success": bool(raw.get("modified_solve_success")),
            "modified_objective_value": raw.get("modified_objective_value"),
            "lp_attribution_status": raw.get("lp_attribution_status"),
            "applied_modifications": raw.get("applied_modifications"),
            "modification_warnings": list(raw.get("modification_warnings") or []),
            "floor_constraints_not_oe_addressable": list(raw.get("floor_constraints_not_oe_addressable") or []),
            "carbon_source_id": raw.get("carbon_source_id"),
            "medium_condition_id": raw.get("medium_condition_id"),
            "mu": raw.get("mu"),
            "error": None,
        }
    )
    return readout


__all__ = ["run_next_oe_candidate_analysis"]
