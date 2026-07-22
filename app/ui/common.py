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
    """Custom styling only - deliberately does NOT set .streamlit/config.toml's [theme]
    table. Setting even a bare `base` there disables Streamlit's native Light/Dark/System
    picker (Main menu -> Theme) for every viewer, which is a bigger loss than the custom
    accent/surface colors are worth. So the accent palette lives here as plain CSS custom
    properties with light values by default and a prefers-color-scheme override, layered
    on top of whichever native theme (including the OS-default "System" choice) the viewer
    has picked - it never fights that picker for control.
    """
    st.markdown(
        """
        <style>
        :root {
            --pc-accent: #0F766E;
            --pc-accent-soft: rgba(15, 118, 110, 0.12);
            --pc-surface: #EEF2F1;
            --pc-surface-strong: #FFFFFF;
            --pc-border: #D8DEDC;
            --pc-text-muted: #52605C;
        }
        @media (prefers-color-scheme: dark) {
            :root {
                --pc-accent: #5EEAD4;
                --pc-accent-soft: rgba(94, 234, 212, 0.16);
                --pc-surface: #172321;
                --pc-surface-strong: #1B2624;
                --pc-border: rgba(255, 255, 255, 0.14);
                --pc-text-muted: #9CB0AB;
            }
        }

        html, body, [class*="css"] {
            font-family: "Segoe UI", "Microsoft YaHei", "Source Sans", Arial, sans-serif;
        }
        code, pre, kbd, [data-testid="stCodeBlock"] {
            font-family: "Cascadia Mono", "Consolas", "SFMono-Regular", monospace;
        }

        .block-container { padding-top: 1.25rem; padding-bottom: 3rem; }

        /* Numeric hierarchy: secretion ratios / growth ratios only mean something in
           comparison to each other, so digits must line up column-to-column. */
        [data-testid="stMetricValue"],
        [data-testid="stDataFrame"] div[role="gridcell"],
        [data-testid="stTable"] td {
            font-variant-numeric: tabular-nums;
        }

        /* Headings: give them presence without the generic Streamlit flatness -
           heavier weight scaled down from h1 to h3. */
        h1 { font-weight: 700; letter-spacing: 0; }
        h2 { font-weight: 600; letter-spacing: 0; }
        h3 { font-weight: 600; letter-spacing: 0; }

        .small-note { color: var(--pc-text-muted); font-size: 0.95rem; line-height: 1.55; }
        .concept-box {
            border-left: 4px solid var(--pc-accent);
            padding: 0.75rem 1rem;
            background: var(--pc-surface);
            border-radius: 0 0.375rem 0.375rem 0;
            margin: 0.25rem 0 1rem 0;
        }

        /* Metric readouts read as a flat label+number pair by default; a hairline
           card gives each one a boundary so a row of metrics doesn't blur together. */
        [data-testid="stMetric"] {
            background: var(--pc-surface-strong);
            border: 1px solid var(--pc-border) !important;
            border-radius: 0.5rem !important;
            padding: 0.85rem 1rem;
        }

        /* Buttons: Streamlit's default has no press feedback and a flat hover -
           add a lift on hover and a settle on press so clicks feel registered. */
        .stButton > button, .stDownloadButton > button {
            transition: transform 120ms ease, box-shadow 120ms ease;
        }
        .stButton > button:hover, .stDownloadButton > button:hover {
            transform: translateY(-1px);
            box-shadow: 0 4px 10px -4px var(--pc-accent-soft);
        }
        .stButton > button:active, .stDownloadButton > button:active {
            transform: translateY(0);
            box-shadow: none;
        }

        /* Sidebar nav: the plain radio-dot list reads as "a form field", not
           navigation. Turn each option into a nav row and highlight the active one. */
        [data-testid="stSidebar"] [data-testid="stRadio"] label {
            padding: 0.4rem 0.6rem;
            border-radius: 0.375rem;
            transition: background-color 120ms ease;
        }
        [data-testid="stSidebar"] [data-testid="stRadio"] label:has(input:checked) {
            background: var(--pc-accent-soft);
            font-weight: 600;
        }

        [data-testid="stDataFrame"] { border-radius: 0.5rem; overflow: hidden; }
        </style>
        """,
        unsafe_allow_html=True,
    )


def page_header() -> None:
    st.title("pcSecYeastSpecies 酵母蛋白分泌模型")
    st.caption("面向生物学专家的跨物种蛋白质组约束模型结果浏览与小规模仿真工具")


NAV_RADIO_KEY = "app_page_nav"
EXPERIMENT_FEEDBACK_PAGE = "实验反馈闭环"
HOMOLOGY_AUDIT_PAGE = "基因命名与同源规则审计"
SHADOW_CROSS_CHECK_PAGE = "Shadow LP一致性验证"
# Streamlit forbids writing to st.session_state[key] once that key's widget has
# rendered in the current script run (raises StreamlitAPIException) - and the nav
# radios below always render first, before any page-specific code runs. So a page
# that wants to navigate elsewhere (e.g. the genome-wide screen's "verify in
# simulation" button) can't set NAV_RADIO_KEY directly. It sets this key instead;
# request_navigation() below applies it before the radios are instantiated on the
# *next* run, which is legal.
PENDING_NAV_KEY = "app_page_nav_pending"

# Grouped so the sidebar reads as "总览 -> 生成候选 -> 复核证据 -> 内部工具" instead of
# nine equally-weighted flat entries. Each group is its own radio (Streamlit has no
# built-in grouped-radio widget); _sync_nav_group keeps exactly one selected across
# all of them so it still behaves like a single page switch.
NAV_GROUPS: tuple[tuple[str, tuple[str, ...]], ...] = (
    ("总览与结果", ("项目总览", "结果浏览")),
    ("候选生成与筛查", ("全基因组KO/OE筛查", "仿真验证")),
    ("证据与历史数据", (EXPERIMENT_FEEDBACK_PAGE, HOMOLOGY_AUDIT_PAGE)),
    ("内部工具", (SHADOW_CROSS_CHECK_PAGE, "运行日志")),
)
DEFAULT_PAGE = NAV_GROUPS[0][1][0]


def _group_widget_key(group_index: int) -> str:
    return f"app_page_nav_group_{group_index}"


def request_navigation(target_page: str) -> None:
    """Queue a jump to target_page for the next rerun. Call this, then st.rerun()."""
    st.session_state[PENDING_NAV_KEY] = target_page


def _sync_nav_group(chosen_index: int) -> None:
    """on_change callback: apply the just-clicked group's choice, clear the rest.

    Without this, clicking a page in one group would leave a stale selection
    visibly highlighted in whichever group was active before it.
    """
    choice = st.session_state.get(_group_widget_key(chosen_index))
    if not choice:
        return
    _set_current_page(choice)


def _set_current_page(target_page: str) -> None:
    st.session_state[NAV_RADIO_KEY] = target_page
    for group_index, (_, options) in enumerate(NAV_GROUPS):
        st.session_state[_group_widget_key(group_index)] = target_page if target_page in options else None


def sidebar_navigation() -> str:
    if PENDING_NAV_KEY in st.session_state:
        _set_current_page(st.session_state.pop(PENDING_NAV_KEY))
    elif NAV_RADIO_KEY not in st.session_state:
        _set_current_page(DEFAULT_PAGE)

    st.sidebar.title("功能导航")
    for group_index, (label, options) in enumerate(NAV_GROUPS):
        st.sidebar.caption(f"**{label}**")
        st.sidebar.radio(
            label,
            options,
            index=None,
            key=_group_widget_key(group_index),
            label_visibility="collapsed",
            on_change=_sync_nav_group,
            args=(group_index,),
        )
    page = str(st.session_state.get(NAV_RADIO_KEY, DEFAULT_PAGE))

    st.sidebar.divider()
    st.sidebar.markdown(
        """
        **推荐使用顺序**

        1. 项目总览
        2. 结果浏览
        3. 全基因组KO/OE筛查（探索：找出哪些基因值得关注，可直接跳去核实/同源证据）
        4. 基因命名与同源规则审计（复核：查看离线 BLAST/RBH 证据和规则迁移状态）
        5. Shadow LP一致性验证
        6. 仿真验证（核实：可从筛查候选行直接跳转过来并自动填好靶点/基因）
        7. 运行日志
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

