from __future__ import annotations

import contextlib

import pandas as pd

from app.services import genome_wide_screen_analysis as analysis
from app.services.genome_wide_screen_registry import RunInfo
from app.ui.views import genome_wide_screen
from app.ui.views.genome_wide_screen import _verify_prefill_field_values


def test_gene_candidate_routes_to_gene_inputs() -> None:
    values = _verify_prefill_field_values("PAS_chr1-4_0047", "KO", "gene")

    assert values["pichia_draft_ko_genes"] == "PAS_chr1-4_0047"
    assert values["pichia_draft_ko_reactions"] == ""
    assert values["pichia_draft_oe_genes"] == ""
    assert values["pichia_draft_oe_reactions"] == ""


def test_catalog_reaction_candidate_routes_to_reaction_inputs() -> None:
    values = _verify_prefill_field_values("sec_PDI1_ERV2_Ero1p_complex_formation", "OE", "catalog_reaction")

    assert values["pichia_draft_oe_reactions"] == "sec_PDI1_ERV2_Ero1p_complex_formation"
    assert values["pichia_draft_oe_genes"] == ""
    assert values["pichia_draft_ko_genes"] == ""
    assert values["pichia_draft_ko_reactions"] == ""


def test_complex_oe_hypothesis_candidate_routes_to_reaction_inputs_not_gene_inputs() -> None:
    """Regression test: candidate_kind != "catalog_reaction" used to be read as "is a gene",
    which would have sent this reaction id into the gene-ID box (and let it silently fail
    GPR resolution) instead of the reaction-ID box where it actually belongs.
    """
    values = _verify_prefill_field_values("ATPS3m_no_1_fwd", "OE", "complex_oe_hypothesis")

    assert values["pichia_draft_oe_reactions"] == "ATPS3m_no_1_fwd"
    assert values["pichia_draft_oe_genes"] == ""


def test_unknown_future_candidate_kind_fails_safe_to_reaction_routing() -> None:
    values = _verify_prefill_field_values("some_id", "KO", "some_new_kind_nobody_added_yet")

    assert values["pichia_draft_ko_reactions"] == "some_id"
    assert values["pichia_draft_ko_genes"] == ""


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

    assert any(df is retry_evidence for df in rendered_dataframes)
    assert any("求解器重试证据" in label for label in expander_labels)


def test_results_section_shows_message_instead_of_crashing_on_empty_target_ids(tmp_path, monkeypatch) -> None:
    """Regression test: a "done" run whose CSV has a header but zero parseable target_id
    rows (e.g. every candidate task failed, or a catalog/complex_hypothesis-scope run
    matched zero qualifying candidates) used to crash with
    StreamlitAPIException("st.tabs must contain at least one tab label") because
    st.tabs([]) was called with an empty label list instead of being guarded.
    """
    csv_path = tmp_path / "gene_tradeoff_rows.csv"
    csv_path.write_text("target_id\n", encoding="utf-8")

    run = RunInfo(
        run_name="ui_all_tasks_failed",
        status="done",
        done=5,
        total=5,
        targets=("hLF",),
        mode="fast",
        pid=None,
        updated_at="2026-07-15T18:55:23",
        is_stale=False,
        scope="gene",
        csv_path=str(csv_path),
        error_count=5,
    )

    tabs_calls: list[object] = []
    warnings: list[str] = []
    captions: list[str] = []
    monkeypatch.setattr(genome_wide_screen.st, "subheader", lambda *a, **k: None)
    monkeypatch.setattr(genome_wide_screen.st, "selectbox", lambda *a, **k: run)
    monkeypatch.setattr(genome_wide_screen.st, "warning", lambda msg, **k: warnings.append(msg))
    monkeypatch.setattr(genome_wide_screen.st, "caption", lambda msg, **k: captions.append(msg))
    monkeypatch.setattr(genome_wide_screen.st, "tabs", lambda labels: tabs_calls.append(labels))

    genome_wide_screen._render_results_section(None, [run])

    assert tabs_calls == []
    assert any("target_id" in msg for msg in warnings)
    assert any("5" in msg for msg in captions)


def test_results_section_omits_error_hint_when_no_errors_recorded(tmp_path, monkeypatch) -> None:
    """A catalog/complex_hypothesis run that legitimately matched zero candidates (no task
    ever failed) should not be told about task failures it did not have."""
    csv_path = tmp_path / "gene_tradeoff_rows.csv"
    csv_path.write_text("target_id\n", encoding="utf-8")

    run = RunInfo(
        run_name="ui_zero_candidates",
        status="done",
        done=0,
        total=0,
        targets=("hLF",),
        mode="fast",
        pid=None,
        updated_at="2026-07-15T18:55:23",
        is_stale=False,
        scope="catalog",
        csv_path=str(csv_path),
        error_count=0,
    )

    captions: list[str] = []
    monkeypatch.setattr(genome_wide_screen.st, "subheader", lambda *a, **k: None)
    monkeypatch.setattr(genome_wide_screen.st, "selectbox", lambda *a, **k: run)
    monkeypatch.setattr(genome_wide_screen.st, "warning", lambda *a, **k: None)
    monkeypatch.setattr(genome_wide_screen.st, "caption", lambda msg, **k: captions.append(msg))
    monkeypatch.setattr(genome_wide_screen.st, "tabs", lambda labels: (_ for _ in ()).throw(AssertionError("st.tabs should not be called")))

    genome_wide_screen._render_results_section(None, [run])

    assert captions == []
