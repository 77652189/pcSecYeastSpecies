from __future__ import annotations

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
from app.ui.views.simulation import render_simulation  # noqa: E402


def main() -> None:
    app_css()
    page_header()
    page = sidebar_navigation()
    if page == "项目总览":
        render_overview()
    elif page == "结果浏览":
        render_results_browser()
    elif page == "仿真验证":
        render_simulation()
    elif page == "运行日志":
        render_logs()
    else:
        render_overview()


if __name__ == "__main__":
    main()
