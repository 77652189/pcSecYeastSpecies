from __future__ import annotations

import importlib
from pathlib import Path

import streamlit as st

APP_ICON = Path(__file__).resolve().parent / "assets" / "pcsecyeast_8502.png"

st.set_page_config(
    page_title="pcSecYeastSpecies",
    page_icon=str(APP_ICON),
    layout="wide",
    initial_sidebar_state="expanded",
)

from app.ui.common import EXPERIMENT_FEEDBACK_PAGE, HOMOLOGY_AUDIT_PAGE, SHADOW_CROSS_CHECK_PAGE, app_css, page_header, sidebar_navigation  # noqa: E402
from app.ui.views.experiment_feedback import render_experiment_feedback  # noqa: E402
from app.ui.views.genome_wide_screen import render_genome_wide_screen  # noqa: E402
from app.ui.views.homology_audit import render_homology_audit  # noqa: E402
from app.ui.views.logs import render_logs  # noqa: E402
from app.ui.views.overview import render_overview  # noqa: E402
from app.ui.views.results import render_results_browser  # noqa: E402
from app.ui.views.shadow_cross_check import render_shadow_cross_check  # noqa: E402
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
    elif page == "全基因组KO/OE筛查":
        render_genome_wide_screen()
    elif page == EXPERIMENT_FEEDBACK_PAGE:
        render_experiment_feedback()
    elif page == HOMOLOGY_AUDIT_PAGE:
        render_homology_audit()
    elif page == SHADOW_CROSS_CHECK_PAGE:
        render_shadow_cross_check()
    elif page == "运行日志":
        render_logs()
    else:
        render_overview()


if __name__ == "__main__":
    main()
