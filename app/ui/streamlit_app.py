from __future__ import annotations

import importlib

import streamlit as st

st.set_page_config(
    page_title="pcSecYeastSpecies",
    page_icon=None,
    layout="wide",
    initial_sidebar_state="expanded",
)

from app.ui.common import app_css, page_header, sidebar_navigation  # noqa: E402
from app.ui.views.logs import render_logs  # noqa: E402
from app.ui.views.overview import render_overview  # noqa: E402
from app.ui.views.results import render_results_browser  # noqa: E402
import app.ui.views.simulation as simulation_view  # noqa: E402


def _render_simulation_reloaded() -> None:
    import app.services.pichia_secretion_service as pichia_service

    importlib.reload(pichia_service)
    importlib.reload(simulation_view)
    simulation_view.render_simulation()


def main() -> None:
    app_css()
    page_header()
    page = sidebar_navigation()
    if page == "项目总览":
        render_overview()
    elif page == "结果浏览":
        render_results_browser()
    elif page == "仿真验证":
        _render_simulation_reloaded()
    elif page == "运行日志":
        render_logs()
    else:
        render_overview()


if __name__ == "__main__":
    main()
