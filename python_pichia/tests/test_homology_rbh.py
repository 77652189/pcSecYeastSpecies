from __future__ import annotations

from pcsec_pichia.homology.cache_schema import BlastHit
from pcsec_pichia.homology.rbh import best_hits_by_query, compute_reciprocal_best_hits


def _hit(
    query: str,
    subject: str,
    *,
    evalue: float = 1e-20,
    bitscore: float = 100.0,
    identity: float = 50.0,
    qcov: float = 90.0,
    scov: float = 90.0,
) -> BlastHit:
    return BlastHit(
        query_id=query,
        subject_id=subject,
        identity_pct=identity,
        alignment_length=90,
        query_length=100,
        subject_length=100,
        evalue=evalue,
        bitscore=bitscore,
        query_coverage=qcov,
        subject_coverage=scov,
    )


def test_best_hits_by_query_uses_deterministic_ranking() -> None:
    hits = [
        _hit("sce1", "pichia_b", evalue=1e-30, bitscore=200, identity=80),
        _hit("sce1", "pichia_a", evalue=1e-30, bitscore=200, identity=80),
        _hit("sce1", "pichia_c", evalue=1e-20, bitscore=300, identity=90),
    ]

    best = best_hits_by_query(hits)

    assert best["sce1"].subject_id == "pichia_a"


def test_compute_reciprocal_best_hits_identifies_rbh_and_paralog_risk() -> None:
    calls = compute_reciprocal_best_hits(
        [_hit("sce1", "pp1"), _hit("sce2", "pp2")],
        [_hit("pp1", "sce1"), _hit("pp2", "sce3")],
    )

    by_query = {call.query_id: call for call in calls}
    assert by_query["sce1"].is_rbh is True
    assert by_query["sce1"].failure_reason is None
    assert by_query["sce2"].is_rbh is False
    assert by_query["sce2"].failure_reason == "reverse_best_is_sce3"
