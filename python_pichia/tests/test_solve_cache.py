from __future__ import annotations

import dataclasses

from pcsec_pichia.simulation import SecretionSimulationResult
from pcsec_pichia.solve_cache import (
    SecretionSolveCacheKey,
    cache_key_digest,
    cached_secretion_solve,
    cache_path,
    load_cached_secretion_result,
    store_secretion_result,
)


def _key(**overrides) -> SecretionSolveCacheKey:
    base = dict(
        target_id="hLF",
        carbon_source_id="glucose",
        media_type=4,
        growth_rate=0.10,
        write_ribosome_translation_constraint=False,
        write_misfolding_constraints=False,
        solver_method="highs-ds",
    )
    base.update(overrides)
    return SecretionSolveCacheKey(**base)


def _result(success: bool = True, objective: float | None = 0.00328) -> SecretionSimulationResult:
    return SecretionSimulationResult(
        success=success,
        target_id="hLF",
        objective_value=objective,
        growth_rate=0.10,
        secretion_flux=objective if success else None,
        status="0" if success else "2",
        message="ok" if success else "infeasible",
        constraint_counts={"eq_total": 100, "ub_total": 1},
        result_status="draft",
        target_parameter_status="draft_matlab_alignment_pending",
        matlab_alignment_status="pending",
        exchange_reaction_id="Ex_protein_hLF",
        build_status="ok",
        lp_sensitivity={"eq_marginals": (1.0, 2.0), "ub_marginals": (3.0,)},
        key_fluxes={"Ex_glc_D": -1.0},
        open_growth_reaction_ids=("BIOMASS", "BIOMASS_glyc"),
        warnings=("w1", "w2"),
    )


def test_cache_key_digest_is_deterministic_and_field_sensitive() -> None:
    base = _key()
    assert cache_key_digest(base) == cache_key_digest(_key())  # 确定性

    # 改动决定结果的任一字段都必须改变摘要（否则会误命中）。
    mutations = (
        _key(target_id="OPN_ALPHA_FULL_PROJECT"),
        _key(carbon_source_id="methanol"),
        _key(media_type=2),
        _key(growth_rate=0.013),
        _key(write_ribosome_translation_constraint=True),
        _key(write_misfolding_constraints=True),
        _key(solver_method="highs-ipm"),
        _key(target_fingerprint="custom-abc"),
        _key(model_variant_fingerprint="ko:PEP4"),  # KO 变体不得与野生型同键
        _key(schema_version="different"),
    )
    base_digest = cache_key_digest(base)
    for mutated in mutations:
        assert cache_key_digest(mutated) != base_digest


def test_secretion_result_round_trips_with_tuple_fidelity(tmp_path) -> None:
    key = _key()
    original = _result()

    store_secretion_result(key, original, cache_dir=tmp_path)
    loaded = load_cached_secretion_result(key, cache_dir=tmp_path)

    assert loaded == original
    # JSON 往返后 tuple 字段必须仍是 tuple（不是 list）。
    assert isinstance(loaded.open_growth_reaction_ids, tuple)
    assert isinstance(loaded.warnings, tuple)
    assert isinstance(loaded.lp_sensitivity["eq_marginals"], tuple)


def test_cached_secretion_solve_miss_then_hit(tmp_path) -> None:
    key = _key()
    calls: list[int] = []

    def compute() -> SecretionSimulationResult:
        calls.append(1)
        return _result()

    first, from_cache_1 = cached_secretion_solve(key, compute, cache_dir=tmp_path)
    assert from_cache_1 is False
    assert len(calls) == 1
    assert cache_path(key, tmp_path).exists()

    second, from_cache_2 = cached_secretion_solve(key, compute, cache_dir=tmp_path)
    assert from_cache_2 is True
    assert len(calls) == 1  # 命中不得再次求解
    assert second == first


def test_cached_secretion_solve_force_recomputes(tmp_path) -> None:
    key = _key()
    calls: list[int] = []

    def compute() -> SecretionSimulationResult:
        calls.append(1)
        return _result()

    cached_secretion_solve(key, compute, cache_dir=tmp_path)
    assert len(calls) == 1
    _, from_cache = cached_secretion_solve(key, compute, cache_dir=tmp_path, force=True)
    assert from_cache is False
    assert len(calls) == 2  # force 必须重算


def test_unsuccessful_solve_is_not_cached(tmp_path) -> None:
    key = _key()
    calls: list[int] = []

    def compute() -> SecretionSimulationResult:
        calls.append(1)
        return _result(success=False, objective=None)

    _, from_cache_1 = cached_secretion_solve(key, compute, cache_dir=tmp_path)
    assert from_cache_1 is False
    assert not cache_path(key, tmp_path).exists()  # 失败结果不固化

    _, from_cache_2 = cached_secretion_solve(key, compute, cache_dir=tmp_path)
    assert from_cache_2 is False
    assert len(calls) == 2  # 每次都重算


def test_digest_mismatch_is_treated_as_miss(tmp_path) -> None:
    key = _key()
    path = cache_path(key, tmp_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    # 结果 payload 带一个错误的 cache_key_digest → 防御性当作未命中。
    path.write_text('{"cache_key_digest": "deadbeef", "result": {}}', encoding="utf-8")

    assert load_cached_secretion_result(key, cache_dir=tmp_path) is None
