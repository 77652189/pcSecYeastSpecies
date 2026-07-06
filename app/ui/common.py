from __future__ import annotations

from pathlib import Path

import pandas as pd
import streamlit as st

from app.adapters.powershell import PowerShellAdapter
from app.core.i18n import (
    HEALTH_COLUMN_LABELS,
    RUN_FILE_COLUMN_LABELS,
    SOPLEX_COLUMN_LABELS,
    category_label,
    file_type_label,
    short_path,
    species_label,
    status_label,
)
from app.services.health import HealthService
from app.services.pichia_secretion_service import discover_project_paths
from app.services.results import ResultCatalog, ResultLoader

def _paths():
    return discover_project_paths(Path(__file__))


PATHS = _paths()


@st.cache_data(show_spinner=False)
def cached_datasets() -> list[dict]:
    return [dataset.model_dump() for dataset in ResultCatalog(PATHS).list_datasets()]


@st.cache_data(show_spinner=False)
def cached_loaded_dataset(dataset_id: str) -> dict:
    catalog = ResultCatalog(PATHS)
    loaded = ResultLoader().load_dataset(catalog.get_dataset(dataset_id))
    return loaded.model_dump()


@st.cache_data(show_spinner=False, ttl=60)
def cached_health() -> dict:
    report = HealthService(PATHS, PowerShellAdapter()).check()
    return report.model_dump()


def dataset_frame() -> pd.DataFrame:
    frame = pd.DataFrame(cached_datasets())
    if frame.empty:
        return frame
    frame["path"] = frame["path"].astype(str)
    frame["species_label"] = frame["species"].map(species_label)
    frame["category_label"] = frame["category"].map(category_label)
    frame["size_kb"] = (frame["size_bytes"] / 1024).round(1)
    frame["suffix"] = frame["suffix"].map(file_type_label)
    return frame


def app_css() -> None:
    st.markdown(
        """
        <style>
        .block-container { padding-top: 1.25rem; padding-bottom: 2rem; }
        [data-testid="stMetricValue"] { font-size: 1.35rem; }
        h1, h2, h3 { letter-spacing: 0; }
        .small-note { color: #475569; font-size: 0.95rem; line-height: 1.55; }
        .concept-box {
            border-left: 4px solid #0f766e;
            padding: 0.75rem 1rem;
            background: #f8fafc;
            margin: 0.25rem 0 1rem 0;
        }
        </style>
        """,
        unsafe_allow_html=True,
    )


def page_header() -> None:
    st.title("pcSecYeastSpecies 酵母蛋白分泌模型")
    st.caption("面向生物学专家的跨物种蛋白质组约束模型结果浏览与小规模仿真工具")


NAV_RADIO_KEY = "app_page_nav"
# Streamlit forbids writing to st.session_state[key] once that key's widget has
# rendered in the current script run (raises StreamlitAPIException) - and the nav
# radio below always renders first, before any page-specific code runs. So a page
# that wants to navigate elsewhere (e.g. the genome-wide screen's "verify in
# simulation" button) can't set NAV_RADIO_KEY directly. It sets this key instead;
# request_navigation() below applies it to NAV_RADIO_KEY before the radio is
# instantiated on the *next* run, which is legal.
PENDING_NAV_KEY = "app_page_nav_pending"


def request_navigation(target_page: str) -> None:
    """Queue a jump to target_page for the next rerun. Call this, then st.rerun()."""
    st.session_state[PENDING_NAV_KEY] = target_page


def sidebar_navigation() -> str:
    if PENDING_NAV_KEY in st.session_state:
        st.session_state[NAV_RADIO_KEY] = st.session_state.pop(PENDING_NAV_KEY)
    st.sidebar.title("演示导航")
    page = st.sidebar.radio(
        "选择功能",
        ["项目总览", "结果浏览", "仿真验证", "全基因组KO/OE筛查", "运行日志"],
        index=0,
        key=NAV_RADIO_KEY,
    )
    st.sidebar.divider()
    st.sidebar.markdown(
        """
        **推荐演示顺序**

        1. 项目总览
        2. 结果浏览
        3. 全基因组KO/OE筛查（探索：找出哪些基因值得关注）
        4. 仿真验证（核实：可从筛查结果候选行直接跳转过来并自动填好靶点/基因）
        5. 运行日志
        """
    )
    if page == "结果浏览":
        st.sidebar.caption("当前页的筛选器在下方，可以按物种、结果主题和关键词过滤。")
    return page


def rename_columns(frame: pd.DataFrame, labels: dict) -> pd.DataFrame:
    return frame.rename(columns={key: value for key, value in labels.items() if key in frame.columns})


def compact_path(value: object) -> str:
    if not value:
        return ""
    text = str(value)
    try:
        path = Path(text)
        if path.is_absolute():
            return str(path.relative_to(PATHS.repo_root))
    except (ValueError, OSError):
        pass
    return text


def _download_file_button(column, path: Path | None, label: str, mime: str) -> None:
    if path is not None and path.exists():
        column.download_button(label, path.read_bytes(), file_name=path.name, mime=mime)
    else:
        column.button(label, disabled=True)

