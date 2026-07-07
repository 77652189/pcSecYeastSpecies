"""Genome-wide KO/OE growth-vs-secretion tradeoff screen.

Historical design rationale and rejected alternatives:
docs/archive/pichia_ko_oe_genome_screen_design_2026-07-02.md

For each gene's KO (or OE) plan, sweeps a set of fixed growth rates (mu) and
reports the max feasible growth rate plus the secretion achieved there. This
reuses the existing run_pcsec_growth_tradeoff / run_pcsec_oe_screen
perturbation machinery (fixed-mu solves) instead of solving for the growth
rate directly, which avoids the mu / protein-budget circular dependency
entirely (see design doc, decision 4).
"""

from __future__ import annotations

from types import SimpleNamespace
from typing import Any, Callable

from pcsec_pichia.screens._prototype_adapter import (
    AminoAcidStoichiometry,
    CobraModel,
    CombinedEnzymeData,
    DEFAULT_SOLVER_TIME_LIMIT_SECONDS,
    MetabolicEnzymeData,
    SecretoryEnzymeData,
    TargetSpec,
    build_supported_target_model,
    build_target_enzymedata,
    run_pcsec_oe_screen,
    solve_pcsec_maximize,
)
from pcsec_pichia.screens.gene_interventions import (
    GeneInterventionPlan,
    plan_gene_knockout,
    plan_gene_overexpression,
)
from pcsec_pichia.screens.gene_perturbation_map import build_reaction_perturbation_mapping

# Coarse sweep for a first-pass genome-wide screen; fine sweep for a
# follow-up shortlist. Both are fractions of the reference growth rate.
FAST_MU_FRACTIONS: tuple[float, ...] = (1.0, 0.5, 0.1)
PRECISE_MU_FRACTIONS: tuple[float, ...] = (1.0, 0.9, 0.75, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.02)
DEFAULT_REFERENCE_GROWTH_RATE = 0.10
DEFAULT_OE_FACTOR = 2.0
DEFAULT_TIMEOUT_RETRY_TIME_LIMIT_SECONDS = DEFAULT_SOLVER_TIME_LIMIT_SECONDS * 3.0
SOLVE_OUTCOME_SUCCESS = "success"
SOLVE_OUTCOME_TIMEOUT = "time_limit_reached"
SOLVE_OUTCOME_PROVEN_INFEASIBLE = "proven_infeasible"
SOLVE_OUTCOME_OTHER_FAILURE = "other_solver_failure"
SOLVE_OUTCOMES: tuple[str, ...] = (
    SOLVE_OUTCOME_SUCCESS,
    SOLVE_OUTCOME_TIMEOUT,
    SOLVE_OUTCOME_PROVEN_INFEASIBLE,
    SOLVE_OUTCOME_OTHER_FAILURE,
)


def mu_points_for_mode(reference_growth_rate: float, mode: str) -> list[float]:
    """Fast mode = 3 coarse points; precise mode = 11 finer points."""
    fractions = FAST_MU_FRACTIONS if mode == "fast" else PRECISE_MU_FRACTIONS
    points = {round(reference_growth_rate * fraction, 6) for fraction in fractions}
    return sorted(value for value in points if value > 0)


def wildtype_secretion_by_mu(
    model: CobraModel,
    exchange_reaction_id: str,
    metabolic: MetabolicEnzymeData,
    secretory: SecretoryEnzymeData,
    combined: CombinedEnzymeData,
    mu_points: list[float],
    write_ribosome_translation_constraint: bool = False,
    write_misfolding_constraints: bool = False,
    time_limit_seconds: float = DEFAULT_SOLVER_TIME_LIMIT_SECONDS,
    timeout_retry_time_limit_seconds: float | None = DEFAULT_TIMEOUT_RETRY_TIME_LIMIT_SECONDS,
) -> dict[float, dict[str, Any]]:
    """Solve the unperturbed model once per mu point; shared baseline for every gene at that mu."""
    baseline_by_mu: dict[float, dict[str, Any]] = {}
    for mu in mu_points:
        fixed_model = model.with_bounds({"BIOMASS": (mu, mu)})
        solved, _counts, retry_metadata = _solve_pcsec_maximize_with_timeout_retry(
            fixed_model,
            exchange_reaction_id,
            metabolic=metabolic,
            secretory=secretory,
            combined=combined,
            mu=mu,
            key_reactions=("BIOMASS", exchange_reaction_id),
            write_ribosome_translation_constraint=write_ribosome_translation_constraint,
            write_misfolding_constraints=write_misfolding_constraints,
            time_limit_seconds=time_limit_seconds,
            timeout_retry_time_limit_seconds=timeout_retry_time_limit_seconds,
        )
        baseline_by_mu[mu] = {
            "success": solved.success,
            "objective_value": solved.objective_value,
            "solver_retry_count": retry_metadata.get("solver_retry_count", 0),
        }
    return baseline_by_mu


def _max_feasible_point(points: list[dict[str, Any]]) -> dict[str, Any] | None:
    feasible = [point for point in points if point.get("success")]
    if not feasible:
        return None
    return max(feasible, key=lambda point: point["mu"])


def _classify_solve_outcome(success: bool, status: object = None, message: object = None) -> str:
    if success:
        return SOLVE_OUTCOME_SUCCESS
    status_text = "" if status is None else str(status).strip().lower()
    message_text = "" if message is None else str(message).strip().lower()
    combined = f"{status_text} {message_text}"
    if status_text == "1" or "time limit" in combined or "time_limit" in combined or "highs status 13" in combined:
        return SOLVE_OUTCOME_TIMEOUT
    if status_text == "2" or "infeasible" in combined:
        return SOLVE_OUTCOME_PROVEN_INFEASIBLE
    return SOLVE_OUTCOME_OTHER_FAILURE


def _retry_enabled(time_limit_seconds: float, timeout_retry_time_limit_seconds: float | None) -> bool:
    return timeout_retry_time_limit_seconds is not None and timeout_retry_time_limit_seconds > time_limit_seconds


def _solve_pcsec_maximize_with_timeout_retry(
    model: CobraModel,
    objective_reaction: str,
    metabolic: MetabolicEnzymeData,
    secretory: SecretoryEnzymeData,
    combined: CombinedEnzymeData,
    mu: float,
    key_reactions: tuple[str, ...] = (),
    write_ribosome_translation_constraint: bool = False,
    write_misfolding_constraints: bool = False,
    time_limit_seconds: float = DEFAULT_SOLVER_TIME_LIMIT_SECONDS,
    timeout_retry_time_limit_seconds: float | None = DEFAULT_TIMEOUT_RETRY_TIME_LIMIT_SECONDS,
) -> tuple[Any, dict[str, int], dict[str, Any]]:
    solved, counts = solve_pcsec_maximize(
        model,
        objective_reaction,
        metabolic=metabolic,
        secretory=secretory,
        combined=combined,
        mu=mu,
        key_reactions=key_reactions,
        write_ribosome_translation_constraint=write_ribosome_translation_constraint,
        write_misfolding_constraints=write_misfolding_constraints,
        time_limit_seconds=time_limit_seconds,
    )
    retry_metadata: dict[str, Any] = {"solver_retry_count": 0}
    initial_outcome = _classify_solve_outcome(solved.success, solved.status, solved.message)
    if initial_outcome != SOLVE_OUTCOME_TIMEOUT or not _retry_enabled(time_limit_seconds, timeout_retry_time_limit_seconds):
        return solved, counts, retry_metadata

    retry_solved, retry_counts = solve_pcsec_maximize(
        model,
        objective_reaction,
        metabolic=metabolic,
        secretory=secretory,
        combined=combined,
        mu=mu,
        key_reactions=key_reactions,
        write_ribosome_translation_constraint=write_ribosome_translation_constraint,
        write_misfolding_constraints=write_misfolding_constraints,
        time_limit_seconds=timeout_retry_time_limit_seconds,
    )
    retry_metadata.update(
        {
            "solver_retry_count": 1,
            "initial_solve_outcome": initial_outcome,
            "initial_status": solved.status,
            "initial_message": solved.message,
            "initial_time_limit_seconds": time_limit_seconds,
            "retry_time_limit_seconds": timeout_retry_time_limit_seconds,
        }
    )
    return retry_solved, retry_counts, retry_metadata


def _oe_rows_have_timeout(rows: list[dict[str, Any]]) -> bool:
    return any(
        _classify_solve_outcome(bool(row.get("success")), row.get("status"), row.get("message")) == SOLVE_OUTCOME_TIMEOUT
        for row in rows
    )


def _run_pcsec_oe_screen_with_timeout_retry(
    model: CobraModel,
    baseline: Any,
    reactions: list[str],
    objective: str,
    metabolic: MetabolicEnzymeData,
    secretory: SecretoryEnzymeData,
    combined: CombinedEnzymeData,
    mu: float,
    factor: float = DEFAULT_OE_FACTOR,
    write_ribosome_translation_constraint: bool = False,
    write_misfolding_constraints: bool = False,
    time_limit_seconds: float = DEFAULT_SOLVER_TIME_LIMIT_SECONDS,
    timeout_retry_time_limit_seconds: float | None = DEFAULT_TIMEOUT_RETRY_TIME_LIMIT_SECONDS,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    rows = run_pcsec_oe_screen(
        model,
        baseline,
        reactions,
        objective,
        metabolic=metabolic,
        secretory=secretory,
        combined=combined,
        mu=mu,
        factor=factor,
        write_ribosome_translation_constraint=write_ribosome_translation_constraint,
        write_misfolding_constraints=write_misfolding_constraints,
        time_limit_seconds=time_limit_seconds,
    )
    retry_metadata: dict[str, Any] = {"solver_retry_count": 0}
    if not _oe_rows_have_timeout(rows) or not _retry_enabled(time_limit_seconds, timeout_retry_time_limit_seconds):
        return rows, retry_metadata

    failure_row = _representative_failure_row(rows) or {}
    retry_rows = run_pcsec_oe_screen(
        model,
        baseline,
        reactions,
        objective,
        metabolic=metabolic,
        secretory=secretory,
        combined=combined,
        mu=mu,
        factor=factor,
        write_ribosome_translation_constraint=write_ribosome_translation_constraint,
        write_misfolding_constraints=write_misfolding_constraints,
        time_limit_seconds=timeout_retry_time_limit_seconds,
    )
    retry_metadata.update(
        {
            "solver_retry_count": 1,
            "initial_solve_outcome": SOLVE_OUTCOME_TIMEOUT,
            "initial_status": failure_row.get("status"),
            "initial_message": failure_row.get("message"),
            "initial_time_limit_seconds": time_limit_seconds,
            "retry_time_limit_seconds": timeout_retry_time_limit_seconds,
        }
    )
    return retry_rows, retry_metadata


def _tradeoff_point(
    mu: float,
    success: bool,
    status: object,
    secretion_flux: float | None,
    message: object = None,
    retry_metadata: dict[str, Any] | None = None,
) -> dict[str, Any]:
    point = {
        "mu": mu,
        "success": success,
        "status": status,
        "secretion_flux": secretion_flux if success else None,
        "solve_outcome": _classify_solve_outcome(success, status, message),
    }
    if message is not None:
        point["message"] = str(message)
    if retry_metadata:
        point.update(retry_metadata)
    return point


def _representative_failure_row(rows: list[dict[str, Any]]) -> dict[str, Any] | None:
    priority = {
        SOLVE_OUTCOME_TIMEOUT: 0,
        SOLVE_OUTCOME_OTHER_FAILURE: 1,
        SOLVE_OUTCOME_PROVEN_INFEASIBLE: 2,
    }
    failed_rows = [row for row in rows if not row.get("success")]
    if not failed_rows:
        return None
    return min(
        failed_rows,
        key=lambda row: priority.get(
            _classify_solve_outcome(False, row.get("status"), row.get("message")),
            priority[SOLVE_OUTCOME_OTHER_FAILURE],
        ),
    )


def _solve_outcome_summary(points: list[dict[str, Any]], best: dict[str, Any] | None) -> dict[str, Any]:
    counts = {outcome: 0 for outcome in SOLVE_OUTCOMES}
    classified_points: list[tuple[dict[str, Any], str]] = []
    for point in points:
        outcome = str(point.get("solve_outcome") or _classify_solve_outcome(point.get("success"), point.get("status"), point.get("message")))
        counts[outcome if outcome in counts else SOLVE_OUTCOME_OTHER_FAILURE] += 1
        classified_points.append((point, outcome if outcome in counts else SOLVE_OUTCOME_OTHER_FAILURE))

    timeout_mu_points = tuple(point["mu"] for point, outcome in classified_points if outcome == SOLVE_OUTCOME_TIMEOUT)
    proven_infeasible_mu_points = tuple(
        point["mu"] for point, outcome in classified_points if outcome == SOLVE_OUTCOME_PROVEN_INFEASIBLE
    )
    other_failure_mu_points = tuple(
        point["mu"] for point, outcome in classified_points if outcome == SOLVE_OUTCOME_OTHER_FAILURE
    )
    timeout_retry_mu_points = tuple(
        point["mu"] for point in points if int(point.get("solver_retry_count") or 0) > 0
    )
    best_mu = best["mu"] if best else None
    upper_bound_timeout = best_mu is None or any(mu > best_mu for mu in timeout_mu_points)
    upper_bound_other_failure = best_mu is None or any(mu > best_mu for mu in other_failure_mu_points)
    if not points:
        interpretation = "not_evaluated"
    elif timeout_mu_points and upper_bound_timeout:
        interpretation = "inconclusive_due_to_timeout"
    elif other_failure_mu_points and upper_bound_other_failure:
        interpretation = "inconclusive_due_to_solver_failure"
    else:
        interpretation = "definitive"

    return {
        "solve_outcome_counts": counts,
        "has_timeout": bool(timeout_mu_points),
        "timeout_mu_points": timeout_mu_points,
        "proven_infeasible_mu_points": proven_infeasible_mu_points,
        "other_solver_failure_mu_points": other_failure_mu_points,
        "solver_retry_count": sum(int(point.get("solver_retry_count") or 0) for point in points),
        "timeout_retry_mu_points": timeout_retry_mu_points,
        "feasibility_interpretation": interpretation,
    }


def _skipped_row(gene_id: str, intervention_type: str, plan: GeneInterventionPlan) -> dict[str, Any]:
    support_status = plan.ko_support_status if intervention_type == "KO" else plan.oe_support_status
    # Even though no LP was solved, the plan already knows which reaction(s) this gene
    # would have touched (explain_only_reactions for OE-blocked complex subunits) -
    # keep that so a "skipped" row still explains *why* and *what* was skipped, instead
    # of leaving affected_reactions/secretory_process blank.
    reactions = list(plan.affected_reactions or plan.explain_only_reactions)
    mapping = build_reaction_perturbation_mapping(reactions[0] if reactions else None, None)
    return {
        "gene_id": gene_id,
        "intervention_type": intervention_type,
        "affected_reactions": reactions,
        "secretory_process": mapping.secretory_process,
        "gpr_role": plan.gpr_role,
        "mapping_confidence": plan.mapping_confidence,
        "support_status": support_status,
        "max_feasible_mu": None,
        "secretion_at_max_feasible_mu": None,
        "tradeoff_points": (),
        **_solve_outcome_summary([], None),
        "skipped_reason": "no_structural_effect",
    }


def _summarize_row(
    gene_id: str,
    intervention_type: str,
    plan: GeneInterventionPlan,
    points: list[dict[str, Any]],
    complex_subunits: dict[str, list[dict[str, object]]] | None,
) -> dict[str, Any]:
    reactions = list(plan.inactive_reactions if intervention_type == "KO" else plan.executable_reactions)
    mapping = build_reaction_perturbation_mapping(reactions[0] if reactions else None, complex_subunits)
    best = _max_feasible_point(points)
    support_status = plan.ko_support_status if intervention_type == "KO" else plan.oe_support_status
    outcome_summary = _solve_outcome_summary(points, best)
    return {
        "gene_id": gene_id,
        "intervention_type": intervention_type,
        "affected_reactions": reactions,
        "secretory_process": mapping.secretory_process,
        "gpr_role": plan.gpr_role,
        "mapping_confidence": plan.mapping_confidence,
        "support_status": support_status,
        "max_feasible_mu": best["mu"] if best else None,
        "secretion_at_max_feasible_mu": best["secretion_flux"] if best else None,
        "tradeoff_points": tuple(points),
        **outcome_summary,
        "skipped_reason": None,
    }


def gene_ko_tradeoff(
    model: CobraModel,
    gene_id: str,
    plan: GeneInterventionPlan,
    exchange_reaction_id: str,
    metabolic: MetabolicEnzymeData,
    secretory: SecretoryEnzymeData,
    combined: CombinedEnzymeData,
    mu_points: list[float],
    complex_subunits: dict[str, list[dict[str, object]]] | None = None,
    write_ribosome_translation_constraint: bool = False,
    write_misfolding_constraints: bool = False,
    time_limit_seconds: float = DEFAULT_SOLVER_TIME_LIMIT_SECONDS,
    timeout_retry_time_limit_seconds: float | None = DEFAULT_TIMEOUT_RETRY_TIME_LIMIT_SECONDS,
) -> dict[str, Any]:
    """Growth-vs-secretion tradeoff curve for knocking out one gene."""
    if not plan.inactive_reactions:
        return _skipped_row(gene_id, "KO", plan)

    ko_model = model.with_bounds({reaction_id: (0.0, 0.0) for reaction_id in plan.inactive_reactions})
    points: list[dict[str, Any]] = []
    for mu in mu_points:
        fixed_model = ko_model.with_bounds({"BIOMASS": (mu, mu)})
        solved, _counts, retry_metadata = _solve_pcsec_maximize_with_timeout_retry(
            fixed_model,
            exchange_reaction_id,
            metabolic=metabolic,
            secretory=secretory,
            combined=combined,
            mu=mu,
            key_reactions=("BIOMASS", exchange_reaction_id),
            write_ribosome_translation_constraint=write_ribosome_translation_constraint,
            write_misfolding_constraints=write_misfolding_constraints,
            time_limit_seconds=time_limit_seconds,
            timeout_retry_time_limit_seconds=timeout_retry_time_limit_seconds,
        )
        points.append(
            _tradeoff_point(mu, solved.success, solved.status, solved.objective_value, solved.message, retry_metadata)
        )
    return _summarize_row(gene_id, "KO", plan, points, complex_subunits)


def gene_oe_tradeoff(
    model: CobraModel,
    gene_id: str,
    plan: GeneInterventionPlan,
    exchange_reaction_id: str,
    metabolic: MetabolicEnzymeData,
    secretory: SecretoryEnzymeData,
    combined: CombinedEnzymeData,
    mu_points: list[float],
    baseline_by_mu: dict[float, dict[str, Any]],
    complex_subunits: dict[str, list[dict[str, object]]] | None = None,
    factor: float = DEFAULT_OE_FACTOR,
    write_ribosome_translation_constraint: bool = False,
    write_misfolding_constraints: bool = False,
    time_limit_seconds: float = DEFAULT_SOLVER_TIME_LIMIT_SECONDS,
    timeout_retry_time_limit_seconds: float | None = DEFAULT_TIMEOUT_RETRY_TIME_LIMIT_SECONDS,
) -> dict[str, Any]:
    """Growth-vs-secretion tradeoff curve for overexpressing one gene's reaction(s)."""
    if not plan.executable_reactions:
        return _skipped_row(gene_id, "OE", plan)

    points: list[dict[str, Any]] = []
    for mu in mu_points:
        fixed_model = model.with_bounds({"BIOMASS": (mu, mu)})
        baseline_entry = baseline_by_mu.get(mu, {"objective_value": None})
        baseline_ns = SimpleNamespace(objective_value=baseline_entry.get("objective_value"))
        oe_rows, retry_metadata = _run_pcsec_oe_screen_with_timeout_retry(
            fixed_model,
            baseline_ns,
            list(plan.executable_reactions),
            exchange_reaction_id,
            metabolic=metabolic,
            secretory=secretory,
            combined=combined,
            mu=mu,
            factor=factor,
            write_ribosome_translation_constraint=write_ribosome_translation_constraint,
            write_misfolding_constraints=write_misfolding_constraints,
            time_limit_seconds=time_limit_seconds,
            timeout_retry_time_limit_seconds=timeout_retry_time_limit_seconds,
        )
        best_reaction_row = max(
            (row for row in oe_rows if row.get("success")),
            key=lambda row: row["objective_value"],
            default=None,
        )
        failure_row = _representative_failure_row(oe_rows)
        points.append(
            _tradeoff_point(
                mu,
                best_reaction_row is not None,
                best_reaction_row["status"]
                if best_reaction_row
                else (failure_row["status"] if failure_row else "no_reactions"),
                best_reaction_row["objective_value"] if best_reaction_row else None,
                best_reaction_row.get("message") if best_reaction_row else (failure_row.get("message") if failure_row else None),
                retry_metadata,
            )
        )
    return _summarize_row(gene_id, "OE", plan, points, complex_subunits)


def catalog_reaction_candidates() -> tuple[dict[str, Any], ...]:
    """Unique reaction candidates from the curated secretion gene catalog, tested for BOTH
    KO and OE regardless of which single field (ko_reaction_id/oe_reaction_id) a curator
    happened to fill in.

    Most catalog entries were curated from literature that only reported one direction
    (usually OE) for a given complex-formation reaction - e.g. PDI1/ERO1/ERV2 only ever
    had oe_reaction_id set. But a reaction can always be validly bounded to zero (KO) or
    scaled up (OE) regardless of which direction
    the original paper happened to test; only screening the one direction the literature
    reported undersells what the model can actually tell you. Many entries are subunits of
    the same curated complex and share one reaction id (e.g. PDI1/ERO1/ERV2 all point to
    sec_PDI1_ERV2_Ero1p_complex_formation) - this collapses them to one row per unique
    reaction id per intervention type, keeping every contributing common_name for display.
    """
    from pcsec_pichia.services.gene_catalog import SECRETION_GENE_CATALOG

    by_reaction: dict[str, dict[str, Any]] = {}
    for entry in SECRETION_GENE_CATALOG:
        reaction_id = entry.oe_reaction_id or entry.ko_reaction_id
        if not reaction_id:
            continue
        existing = by_reaction.get(reaction_id)
        if existing:
            if entry.common_name not in existing["common_names"]:
                existing["common_names"].append(entry.common_name)
        else:
            by_reaction[reaction_id] = {
                "reaction_id": reaction_id,
                "category": entry.category,
                "common_names": [entry.common_name],
            }

    candidates: list[dict[str, Any]] = []
    for item in by_reaction.values():
        common_name = "/".join(item.pop("common_names"))
        for intervention_type in ("KO", "OE"):
            candidates.append({**item, "intervention_type": intervention_type, "common_name": common_name})
    return tuple(candidates)


def _summarize_catalog_row(
    reaction_id: str,
    common_name: str,
    category: str,
    intervention_type: str,
    points: list[dict[str, Any]],
    complex_subunits: dict[str, list[dict[str, object]]] | None,
    candidate_kind: str = "catalog_reaction",
    hypothesis_note: str = "",
) -> dict[str, Any]:
    mapping = build_reaction_perturbation_mapping(reaction_id, complex_subunits)
    best = _max_feasible_point(points)
    outcome_summary = _solve_outcome_summary(points, best)
    return {
        "gene_id": reaction_id,
        "common_name": common_name,
        "candidate_kind": candidate_kind,
        "category": category,
        "intervention_type": intervention_type,
        "affected_reactions": (reaction_id,),
        "secretory_process": mapping.secretory_process,
        "gpr_role": "catalog_curated",
        "mapping_confidence": "curated",
        "support_status": f"catalog_{intervention_type.lower()}_direct_reaction",
        "max_feasible_mu": best["mu"] if best else None,
        "secretion_at_max_feasible_mu": best["secretion_flux"] if best else None,
        "tradeoff_points": tuple(points),
        **outcome_summary,
        "skipped_reason": None,
        "hypothesis_note": hypothesis_note,
    }


def reaction_ko_tradeoff(
    model: CobraModel,
    reaction_id: str,
    common_name: str,
    category: str,
    exchange_reaction_id: str,
    metabolic: MetabolicEnzymeData,
    secretory: SecretoryEnzymeData,
    combined: CombinedEnzymeData,
    mu_points: list[float],
    complex_subunits: dict[str, list[dict[str, object]]] | None = None,
    write_ribosome_translation_constraint: bool = False,
    write_misfolding_constraints: bool = False,
    candidate_kind: str = "catalog_reaction",
    hypothesis_note: str = "",
    time_limit_seconds: float = DEFAULT_SOLVER_TIME_LIMIT_SECONDS,
    timeout_retry_time_limit_seconds: float | None = DEFAULT_TIMEOUT_RETRY_TIME_LIMIT_SECONDS,
) -> dict[str, Any]:
    """Growth-vs-secretion tradeoff for directly knocking out one curated reaction.

    Unlike gene_ko_tradeoff, this does not resolve a gene through GPR first - most
    curated catalog entries (e.g. PDI1/ERO1/ERV2) are complex-level MATLAB pseudo-reactions
    with no single resolvable gene_id, so the reaction bound is set directly.
    """
    ko_model = model.with_bounds({reaction_id: (0.0, 0.0)})
    points: list[dict[str, Any]] = []
    for mu in mu_points:
        fixed_model = ko_model.with_bounds({"BIOMASS": (mu, mu)})
        solved, _counts, retry_metadata = _solve_pcsec_maximize_with_timeout_retry(
            fixed_model,
            exchange_reaction_id,
            metabolic=metabolic,
            secretory=secretory,
            combined=combined,
            mu=mu,
            key_reactions=("BIOMASS", exchange_reaction_id),
            write_ribosome_translation_constraint=write_ribosome_translation_constraint,
            write_misfolding_constraints=write_misfolding_constraints,
            time_limit_seconds=time_limit_seconds,
            timeout_retry_time_limit_seconds=timeout_retry_time_limit_seconds,
        )
        points.append(
            _tradeoff_point(mu, solved.success, solved.status, solved.objective_value, solved.message, retry_metadata)
        )
    return _summarize_catalog_row(
        reaction_id, common_name, category, "KO", points, complex_subunits, candidate_kind, hypothesis_note
    )


def reaction_oe_tradeoff(
    model: CobraModel,
    reaction_id: str,
    common_name: str,
    category: str,
    exchange_reaction_id: str,
    metabolic: MetabolicEnzymeData,
    secretory: SecretoryEnzymeData,
    combined: CombinedEnzymeData,
    mu_points: list[float],
    baseline_by_mu: dict[float, dict[str, Any]],
    complex_subunits: dict[str, list[dict[str, object]]] | None = None,
    factor: float = DEFAULT_OE_FACTOR,
    write_ribosome_translation_constraint: bool = False,
    write_misfolding_constraints: bool = False,
    candidate_kind: str = "catalog_reaction",
    hypothesis_note: str = "",
    time_limit_seconds: float = DEFAULT_SOLVER_TIME_LIMIT_SECONDS,
    timeout_retry_time_limit_seconds: float | None = DEFAULT_TIMEOUT_RETRY_TIME_LIMIT_SECONDS,
) -> dict[str, Any]:
    """Growth-vs-secretion tradeoff for directly overexpressing one curated reaction.

    Same rationale as reaction_ko_tradeoff: skips gene/GPR resolution and calls
    run_pcsec_oe_screen directly on the curated reaction_id.
    """
    points: list[dict[str, Any]] = []
    for mu in mu_points:
        fixed_model = model.with_bounds({"BIOMASS": (mu, mu)})
        baseline_entry = baseline_by_mu.get(mu, {"objective_value": None})
        baseline_ns = SimpleNamespace(objective_value=baseline_entry.get("objective_value"))
        oe_rows, retry_metadata = _run_pcsec_oe_screen_with_timeout_retry(
            fixed_model,
            baseline_ns,
            [reaction_id],
            exchange_reaction_id,
            metabolic=metabolic,
            secretory=secretory,
            combined=combined,
            mu=mu,
            factor=factor,
            write_ribosome_translation_constraint=write_ribosome_translation_constraint,
            write_misfolding_constraints=write_misfolding_constraints,
            time_limit_seconds=time_limit_seconds,
            timeout_retry_time_limit_seconds=timeout_retry_time_limit_seconds,
        )
        best_row = oe_rows[0] if oe_rows and oe_rows[0].get("success") else None
        failure_row = _representative_failure_row(oe_rows)
        points.append(
            _tradeoff_point(
                mu,
                best_row is not None,
                best_row["status"] if best_row else (failure_row["status"] if failure_row else "no_reaction"),
                best_row["objective_value"] if best_row else None,
                best_row.get("message") if best_row else (failure_row.get("message") if failure_row else None),
                retry_metadata,
            )
        )
    return _summarize_catalog_row(
        reaction_id, common_name, category, "OE", points, complex_subunits, candidate_kind, hypothesis_note
    )


def run_catalog_reaction_tradeoff_screen(
    model: CobraModel,
    target: TargetSpec,
    amino_acids: AminoAcidStoichiometry,
    metabolic: MetabolicEnzymeData,
    secretory: SecretoryEnzymeData,
    combined: CombinedEnzymeData,
    mode: str = "fast",
    reference_growth_rate: float = DEFAULT_REFERENCE_GROWTH_RATE,
    factor: float = DEFAULT_OE_FACTOR,
    write_ribosome_translation_constraint: bool = False,
    write_misfolding_constraints: bool = False,
    progress_callback: Callable[[int, int, str], None] | None = None,
    time_limit_seconds: float = DEFAULT_SOLVER_TIME_LIMIT_SECONDS,
    timeout_retry_time_limit_seconds: float | None = DEFAULT_TIMEOUT_RETRY_TIME_LIMIT_SECONDS,
) -> dict[str, Any]:
    """Run KO+OE tradeoff screening over the curated catalog's unique reactions for one target.

    Second, much smaller screen entry point alongside run_genome_wide_tradeoff_screen:
    instead of iterating all ~1025 model genes, iterates the ~30 unique reactions named
    in SECRETION_GENE_CATALOG (literature-curated secretion engineering targets). Same
    row shape/columns as the gene-level screen (candidate_kind distinguishes them) so the
    existing analysis/UI code applies unchanged.
    """
    build = build_supported_target_model(model, target, amino_acids)
    if not build.supported or build.model is None or build.exchange_reaction_id is None:
        raise ValueError(f"Target {target.target_id!r} is not supported by this model (build_status={build.build_status!r}).")

    target_enzymedata = build_target_enzymedata(target, build.model, secretory)
    target_secretory = secretory.with_reaction_coefficients(target_enzymedata.reaction_coefficients)
    target_combined = combined.with_target(target_enzymedata)
    complex_subunits = getattr(secretory, "complex_subunits", None)

    mu_points = mu_points_for_mode(reference_growth_rate, mode)
    baseline_by_mu = wildtype_secretion_by_mu(
        build.model,
        build.exchange_reaction_id,
        metabolic,
        target_secretory,
        target_combined,
        mu_points,
        write_ribosome_translation_constraint,
        write_misfolding_constraints,
        time_limit_seconds,
        timeout_retry_time_limit_seconds,
    )
    wildtype_best = _max_feasible_point(
        [{"mu": mu, "success": entry["success"], "secretion_flux": entry["objective_value"]} for mu, entry in baseline_by_mu.items()]
    )

    candidates = catalog_reaction_candidates()
    rows: list[dict[str, Any]] = []
    total = len(candidates)
    for index, candidate in enumerate(candidates):
        if candidate["intervention_type"] == "KO":
            row = reaction_ko_tradeoff(
                build.model,
                candidate["reaction_id"],
                candidate["common_name"],
                candidate["category"],
                build.exchange_reaction_id,
                metabolic,
                target_secretory,
                target_combined,
                mu_points,
                complex_subunits,
                write_ribosome_translation_constraint,
                write_misfolding_constraints,
                time_limit_seconds=time_limit_seconds,
                timeout_retry_time_limit_seconds=timeout_retry_time_limit_seconds,
            )
        else:
            row = reaction_oe_tradeoff(
                build.model,
                candidate["reaction_id"],
                candidate["common_name"],
                candidate["category"],
                build.exchange_reaction_id,
                metabolic,
                target_secretory,
                target_combined,
                mu_points,
                baseline_by_mu,
                complex_subunits,
                factor,
                write_ribosome_translation_constraint,
                write_misfolding_constraints,
                time_limit_seconds=time_limit_seconds,
                timeout_retry_time_limit_seconds=timeout_retry_time_limit_seconds,
            )
        row["target_id"] = target.target_id
        _attach_wildtype_comparison(row, wildtype_best)
        rows.append(row)
        if progress_callback:
            progress_callback(index + 1, total, candidate["reaction_id"])

    return {
        "target_id": target.target_id,
        "mode": mode,
        "reference_growth_rate": reference_growth_rate,
        "mu_points": mu_points,
        "wildtype_max_feasible_mu": wildtype_best["mu"] if wildtype_best else None,
        "wildtype_secretion_at_max_feasible_mu": wildtype_best["secretion_flux"] if wildtype_best else None,
        "gene_count": total,
        "rows": rows,
    }


COMPLEX_OE_HYPOTHESIS_ASSUMPTION = (
    "假设性结果：数值假设该反应涉及的复合体所有亚基能按比例协同过表达，用同一个kcat/反应上限"
    f"乘数（factor={DEFAULT_OE_FACTOR}x）代表整体产能提升；不代表对任何单个基因做过表达就能达到"
    "此效果——单基因过表达默认不increment复合体capacity，正是这些基因原本被跳过OE测试的原因。"
    "没有验证亚基化学计量或哪个亚基是真正限速步骤，这两项信息目前不在模型里，也不在本代码库中，"
    "需要外部文献/蛋白质组学数据（比如亚基丰度数据库）才能确认。"
)


def resolve_complex_subunit_oe_hypothesis_candidates(
    model: CobraModel,
    gene_ids: list[str],
    complex_subunits: dict[str, list[dict[str, object]]] | None = None,
) -> tuple[dict[str, Any], ...]:
    """Re-resolve reaction id(s) for genes whose OE was skipped as complex_subunit.

    Takes gene ids (e.g. from genome_wide_screen_analysis's ko_yield_down rows that also
    have gpr_role == "complex_subunit"), re-derives their GPR-associated reaction(s) via
    plan_gene_overexpression fresh from the live model rather than trusting affected_reactions
    in an existing CSV (older screen runs predate that field being populated for skipped
    rows). Dedupes by reaction id the same way catalog_reaction_candidates() does, since
    multiple subunit genes commonly share one reaction.
    """
    by_reaction: dict[str, dict[str, Any]] = {}
    for gene_id in gene_ids:
        plan = plan_gene_overexpression(model, gene_id, complex_subunits=complex_subunits)
        for reaction_id in plan.explain_only_reactions:
            existing = by_reaction.get(reaction_id)
            if existing:
                if gene_id not in existing["gene_ids"]:
                    existing["gene_ids"].append(gene_id)
            else:
                by_reaction[reaction_id] = {"reaction_id": reaction_id, "gene_ids": [gene_id]}

    candidates = []
    for item in by_reaction.values():
        gene_ids_for_reaction = item.pop("gene_ids")
        candidates.append({**item, "common_name": "/".join(gene_ids_for_reaction), "category": "complex_oe_hypothesis"})
    return tuple(candidates)


def run_complex_subunit_oe_hypothesis_screen(
    model: CobraModel,
    target: TargetSpec,
    amino_acids: AminoAcidStoichiometry,
    metabolic: MetabolicEnzymeData,
    secretory: SecretoryEnzymeData,
    combined: CombinedEnzymeData,
    gene_ids: list[str],
    mode: str = "fast",
    reference_growth_rate: float = DEFAULT_REFERENCE_GROWTH_RATE,
    factor: float = DEFAULT_OE_FACTOR,
    write_ribosome_translation_constraint: bool = False,
    write_misfolding_constraints: bool = False,
    progress_callback: Callable[[int, int, str], None] | None = None,
    time_limit_seconds: float = DEFAULT_SOLVER_TIME_LIMIT_SECONDS,
    timeout_retry_time_limit_seconds: float | None = DEFAULT_TIMEOUT_RETRY_TIME_LIMIT_SECONDS,
) -> dict[str, Any]:
    """OE-only hypothesis screen for complex-subunit genes whose KO hurts secretion.

    These genes already have real KO data (that's how they were selected as candidates -
    see genome_wide_screen_analysis.complex_subunit_oe_hypothesis_candidates); what's
    missing is OE, which plan_gene_overexpression() correctly refuses to run at gene level
    for complex subunits. This tests the reaction directly instead - see
    COMPLEX_OE_HYPOTHESIS_ASSUMPTION for exactly what that does and does not prove.
    """
    build = build_supported_target_model(model, target, amino_acids)
    if not build.supported or build.model is None or build.exchange_reaction_id is None:
        raise ValueError(f"Target {target.target_id!r} is not supported by this model (build_status={build.build_status!r}).")

    target_enzymedata = build_target_enzymedata(target, build.model, secretory)
    target_secretory = secretory.with_reaction_coefficients(target_enzymedata.reaction_coefficients)
    target_combined = combined.with_target(target_enzymedata)
    complex_subunits = getattr(secretory, "complex_subunits", None)

    mu_points = mu_points_for_mode(reference_growth_rate, mode)
    baseline_by_mu = wildtype_secretion_by_mu(
        build.model,
        build.exchange_reaction_id,
        metabolic,
        target_secretory,
        target_combined,
        mu_points,
        write_ribosome_translation_constraint,
        write_misfolding_constraints,
        time_limit_seconds,
        timeout_retry_time_limit_seconds,
    )
    wildtype_best = _max_feasible_point(
        [{"mu": mu, "success": entry["success"], "secretion_flux": entry["objective_value"]} for mu, entry in baseline_by_mu.items()]
    )

    candidates = resolve_complex_subunit_oe_hypothesis_candidates(build.model, gene_ids, complex_subunits)
    rows: list[dict[str, Any]] = []
    total = len(candidates)
    for index, candidate in enumerate(candidates):
        row = reaction_oe_tradeoff(
            build.model,
            candidate["reaction_id"],
            candidate["common_name"],
            candidate["category"],
            build.exchange_reaction_id,
            metabolic,
            target_secretory,
            target_combined,
            mu_points,
            baseline_by_mu,
            complex_subunits,
            factor,
            write_ribosome_translation_constraint,
            write_misfolding_constraints,
            candidate_kind="complex_oe_hypothesis",
            hypothesis_note=COMPLEX_OE_HYPOTHESIS_ASSUMPTION,
            time_limit_seconds=time_limit_seconds,
            timeout_retry_time_limit_seconds=timeout_retry_time_limit_seconds,
        )
        row["target_id"] = target.target_id
        _attach_wildtype_comparison(row, wildtype_best)
        rows.append(row)
        if progress_callback:
            progress_callback(index + 1, total, candidate["reaction_id"])

    return {
        "target_id": target.target_id,
        "mode": mode,
        "reference_growth_rate": reference_growth_rate,
        "mu_points": mu_points,
        "wildtype_max_feasible_mu": wildtype_best["mu"] if wildtype_best else None,
        "wildtype_secretion_at_max_feasible_mu": wildtype_best["secretion_flux"] if wildtype_best else None,
        "gene_count": total,
        "rows": rows,
    }


def _attach_wildtype_comparison(row: dict[str, Any], wildtype_best: dict[str, Any] | None) -> dict[str, Any]:
    wt_mu = wildtype_best["mu"] if wildtype_best else None
    wt_secretion = wildtype_best["secretion_flux"] if wildtype_best else None
    row["wildtype_max_feasible_mu"] = wt_mu
    row["wildtype_secretion_at_max_feasible_mu"] = wt_secretion
    row["growth_retention_ratio"] = (
        row["max_feasible_mu"] / wt_mu if row["max_feasible_mu"] is not None and wt_mu else None
    )
    row["secretion_ratio_vs_wildtype"] = (
        row["secretion_at_max_feasible_mu"] / wt_secretion
        if row["secretion_at_max_feasible_mu"] is not None and wt_secretion
        else None
    )
    return row


def run_genome_wide_tradeoff_screen(
    model: CobraModel,
    target: TargetSpec,
    amino_acids: AminoAcidStoichiometry,
    metabolic: MetabolicEnzymeData,
    secretory: SecretoryEnzymeData,
    combined: CombinedEnzymeData,
    gene_ids: list[str],
    mode: str = "fast",
    reference_growth_rate: float = DEFAULT_REFERENCE_GROWTH_RATE,
    factor: float = DEFAULT_OE_FACTOR,
    write_ribosome_translation_constraint: bool = False,
    write_misfolding_constraints: bool = False,
    progress_callback: Callable[[int, int, str], None] | None = None,
    time_limit_seconds: float = DEFAULT_SOLVER_TIME_LIMIT_SECONDS,
    timeout_retry_time_limit_seconds: float | None = DEFAULT_TIMEOUT_RETRY_TIME_LIMIT_SECONDS,
) -> dict[str, Any]:
    """Run KO+OE growth/secretion tradeoff screening for a list of genes against one target."""
    build = build_supported_target_model(model, target, amino_acids)
    if not build.supported or build.model is None or build.exchange_reaction_id is None:
        raise ValueError(f"Target {target.target_id!r} is not supported by this model (build_status={build.build_status!r}).")

    target_enzymedata = build_target_enzymedata(target, build.model, secretory)
    target_secretory = secretory.with_reaction_coefficients(target_enzymedata.reaction_coefficients)
    target_combined = combined.with_target(target_enzymedata)
    complex_subunits = getattr(secretory, "complex_subunits", None)

    mu_points = mu_points_for_mode(reference_growth_rate, mode)
    baseline_by_mu = wildtype_secretion_by_mu(
        build.model,
        build.exchange_reaction_id,
        metabolic,
        target_secretory,
        target_combined,
        mu_points,
        write_ribosome_translation_constraint,
        write_misfolding_constraints,
        time_limit_seconds,
        timeout_retry_time_limit_seconds,
    )
    wildtype_best = _max_feasible_point(
        [{"mu": mu, "success": entry["success"], "secretion_flux": entry["objective_value"]} for mu, entry in baseline_by_mu.items()]
    )

    rows: list[dict[str, Any]] = []
    total = len(gene_ids)
    for index, gene_id in enumerate(gene_ids):
        ko_plan = plan_gene_knockout(build.model, gene_id)
        oe_plan = plan_gene_overexpression(build.model, gene_id, complex_subunits=complex_subunits)

        ko_row = gene_ko_tradeoff(
            build.model,
            gene_id,
            ko_plan,
            build.exchange_reaction_id,
            metabolic,
            target_secretory,
            target_combined,
            mu_points,
            complex_subunits,
            write_ribosome_translation_constraint,
            write_misfolding_constraints,
            time_limit_seconds,
            timeout_retry_time_limit_seconds,
        )
        oe_row = gene_oe_tradeoff(
            build.model,
            gene_id,
            oe_plan,
            build.exchange_reaction_id,
            metabolic,
            target_secretory,
            target_combined,
            mu_points,
            baseline_by_mu,
            complex_subunits,
            factor,
            write_ribosome_translation_constraint,
            write_misfolding_constraints,
            time_limit_seconds,
            timeout_retry_time_limit_seconds,
        )
        for row in (ko_row, oe_row):
            row["target_id"] = target.target_id
            _attach_wildtype_comparison(row, wildtype_best)
            rows.append(row)

        if progress_callback:
            progress_callback(index + 1, total, gene_id)

    return {
        "target_id": target.target_id,
        "mode": mode,
        "reference_growth_rate": reference_growth_rate,
        "mu_points": mu_points,
        "wildtype_max_feasible_mu": wildtype_best["mu"] if wildtype_best else None,
        "wildtype_secretion_at_max_feasible_mu": wildtype_best["secretion_flux"] if wildtype_best else None,
        "gene_count": total,
        "rows": rows,
    }


__all__ = [
    "FAST_MU_FRACTIONS",
    "PRECISE_MU_FRACTIONS",
    "DEFAULT_REFERENCE_GROWTH_RATE",
    "DEFAULT_OE_FACTOR",
    "DEFAULT_TIMEOUT_RETRY_TIME_LIMIT_SECONDS",
    "COMPLEX_OE_HYPOTHESIS_ASSUMPTION",
    "mu_points_for_mode",
    "wildtype_secretion_by_mu",
    "gene_ko_tradeoff",
    "gene_oe_tradeoff",
    "run_genome_wide_tradeoff_screen",
    "catalog_reaction_candidates",
    "reaction_ko_tradeoff",
    "reaction_oe_tradeoff",
    "run_catalog_reaction_tradeoff_screen",
    "resolve_complex_subunit_oe_hypothesis_candidates",
    "run_complex_subunit_oe_hypothesis_screen",
]
