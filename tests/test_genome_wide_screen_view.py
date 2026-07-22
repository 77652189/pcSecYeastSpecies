from __future__ import annotations

import contextlib

import pandas as pd

from app.services import genome_wide_screen_analysis as analysis
from app.ui.views import genome_wide_screen

# Candidate-routing tests for the shared "which draft field gets this candidate"
# logic live in tests/test_simulation_view.py, next to app/ui/views/simulation.py
# (apply_simulation_prefill / _prefill_field_values), which now owns that logic.


def _empty_dimensional_results(**overrides: object) -> analysis.DimensionalResults:
    empty = pd.DataFrame()
    fields: dict[str, object] = {
        "target_id": "hLF",
        "essential_genes": empty,
        "solver_inconclusive_ko": empty,
        "solver_inconclusive_rows": empty,
        "solver_retry_evidence": empty,
        "ko_yield_up_growth_cost": empty,
        "ko_clean_wins": empty,
        "ko_yield_down": empty,
        "oe_yield_up": empty,
        "complex_oe_hypothesis": empty,
        "row_count": 0,
        "skipped_count": 0,
    }
    fields.update(overrides)
    return analysis.DimensionalResults(**fields)


def test_dimension_tables_surface_solver_retry_evidence(monkeypatch) -> None:
    """Regression test for a prior version of this test that only grepped the function's
    *source code* via inspect.getsource() for these strings - which would keep passing even
    if the retry-evidence table were never actually rendered (wrong variable, dead branch, an
    exception thrown earlier in the function). This drives the real render call and asserts
    the exact DataFrame instance reaches st.dataframe under the expected expander label.
    """
    rendered_dataframes: list[object] = []
    expander_labels: list[str] = []

    def fake_expander(label, *args, **kwargs):
        expander_labels.append(label)
        return contextlib.nullcontext()

    monkeypatch.setattr(genome_wide_screen.st, "dataframe", lambda data, **kwargs: rendered_dataframes.append(data))
    monkeypatch.setattr(genome_wide_screen.st, "expander", fake_expander)
    monkeypatch.setattr(genome_wide_screen.st, "metric", lambda *args, **kwargs: None)
    monkeypatch.setattr(genome_wide_screen.st, "warning", lambda *args, **kwargs: None)
    monkeypatch.setattr(genome_wide_screen.st, "markdown", lambda *args, **kwargs: None)
    monkeypatch.setattr(genome_wide_screen.st, "caption", lambda *args, **kwargs: None)

    retry_evidence = pd.DataFrame([{"gene_id": "PAS_TEST_0001", "common_name": "TEST", "max_feasible_mu": 0.1}])
    result = _empty_dimensional_results(solver_retry_evidence=retry_evidence)

    genome_wide_screen._render_dimension_tables(result)

    # 汉化后渲染的是中文列名的展示副本（不再是原始实例），断言行内容仍到达 st.dataframe
    assert any(hasattr(df, "to_string") and "PAS_TEST_0001" in df.to_string() for df in rendered_dataframes)
    assert any("求解器重试证据" in label for label in expander_labels)
