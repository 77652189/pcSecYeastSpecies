from __future__ import annotations

from types import SimpleNamespace

import numpy as np
from scipy import sparse

import pcsec_pichia.probe._prototype as prototype
from pcsec_pichia.probe._prototype import CobraModel


def _tiny_model() -> CobraModel:
    return CobraModel(
        source_file="tiny.mat",
        rxns=["R_IN", "R_OUT"],
        mets=["A_c"],
        genes=[],
        lb=np.array([0.0, 0.0], dtype=float),
        ub=np.array([1000.0, 1000.0], dtype=float),
        b=np.zeros(1, dtype=float),
        s_matrix=sparse.csc_matrix([[1.0, -1.0]]),
        rules=["", ""],
        gr_rules=["", ""],
    )


def test_solve_pcsec_maximize_passes_default_and_explicit_time_limit(monkeypatch):
    captured_time_limits: list[float] = []

    def fake_build_pcsec_constraint_matrices(*args, **kwargs):
        return (
            sparse.csc_matrix([[1.0, -1.0]]),
            np.zeros(1, dtype=float),
            sparse.csc_matrix((0, 2)),
            np.zeros(0, dtype=float),
            {"stoichiometric": 1},
        )

    def fake_linprog(*args, **kwargs):
        captured_time_limits.append(kwargs["options"]["time_limit"])
        return SimpleNamespace(success=False, x=None, status=1, message="Time limit reached")

    monkeypatch.setattr(prototype, "build_pcsec_constraint_matrices", fake_build_pcsec_constraint_matrices)
    monkeypatch.setattr(prototype, "linprog", fake_linprog)

    prototype.solve_pcsec_maximize(
        _tiny_model(),
        "R_OUT",
        metabolic=object(),  # type: ignore[arg-type]
        secretory=object(),  # type: ignore[arg-type]
        combined=object(),  # type: ignore[arg-type]
        mu=0.1,
    )
    prototype.solve_pcsec_maximize(
        _tiny_model(),
        "R_OUT",
        metabolic=object(),  # type: ignore[arg-type]
        secretory=object(),  # type: ignore[arg-type]
        combined=object(),  # type: ignore[arg-type]
        mu=0.1,
        time_limit_seconds=0.25,
    )

    assert captured_time_limits == [prototype.DEFAULT_SOLVER_TIME_LIMIT_SECONDS, 0.25]
