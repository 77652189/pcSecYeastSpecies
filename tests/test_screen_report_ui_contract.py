from __future__ import annotations

import inspect

from app.ui.views import genome_wide_screen


def test_genome_wide_screen_imports_audited_report_service_not_legacy_llm_generator() -> None:
    source = inspect.getsource(genome_wide_screen)

    assert "screen_report_service" in source
    assert "get_default_generator" not in source
    assert "生成研发建议报告" in source
    assert "build_fact_pack_for_runs" in source


def test_screen_report_button_is_explicit_and_not_called_on_page_load() -> None:
    source = inspect.getsource(genome_wide_screen._render_llm_report_section)

    assert "st.button(\"生成研发建议报告\"" in source
    assert "generate_judged_screen_report" in source
    assert source.index("st.button(\"生成研发建议报告\"") < source.index("generate_judged_screen_report")
    assert "OPENAI_API_KEY" in source
