from __future__ import annotations

from collections.abc import Iterable

from pcsec_pichia.homology.cache_schema import BlastHit, ReciprocalBestHit


def best_hits_by_query(hits: Iterable[BlastHit]) -> dict[str, BlastHit]:
    """Pick best hit by evalue, bitscore, identity, coverage, and subject id."""

    best: dict[str, BlastHit] = {}
    for hit in hits:
        current = best.get(hit.query_id)
        if current is None or _sort_key(hit) < _sort_key(current):
            best[hit.query_id] = hit
    return best


def compute_reciprocal_best_hits(
    forward_hits: Iterable[BlastHit],
    reverse_hits: Iterable[BlastHit],
) -> tuple[ReciprocalBestHit, ...]:
    """Return RBH calls for forward query -> subject candidates."""

    forward_best = best_hits_by_query(forward_hits)
    reverse_best = best_hits_by_query(reverse_hits)
    calls: list[ReciprocalBestHit] = []
    for query_id in sorted(forward_best):
        forward = forward_best[query_id]
        reverse = reverse_best.get(forward.subject_id)
        is_rbh = reverse is not None and reverse.subject_id == query_id
        if reverse is None:
            failure = "no_reverse_hit"
        elif not is_rbh:
            failure = f"reverse_best_is_{reverse.subject_id}"
        else:
            failure = None
        calls.append(
            ReciprocalBestHit(
                query_id=query_id,
                subject_id=forward.subject_id,
                is_rbh=is_rbh,
                forward_hit=forward,
                reverse_hit=reverse,
                failure_reason=failure,
            )
        )
    return tuple(calls)


def _sort_key(hit: BlastHit) -> tuple[float, float, float, float, str]:
    coverage = min(hit.query_coverage, hit.subject_coverage)
    return (hit.evalue, -hit.bitscore, -hit.identity_pct, -coverage, hit.subject_id)
