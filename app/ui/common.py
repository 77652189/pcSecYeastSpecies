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
from app.core.paths import ProjectPaths
from app.services.cds_design import CdsDesignService
from app.services.health import HealthService
from app.services.opn import OpnCandidateCatalog
from app.services.opn_signal_peptides import OpnSignalPeptideCandidateSource
from app.services.results import ResultCatalog, ResultLoader
from app.services.signal_peptide_library import SignalPeptideLibraryService
from app.services.signal_peptide_screening import SignalPeptideScreeningService

def _paths() -> ProjectPaths:
    return ProjectPaths.discover(Path(__file__))


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


@st.cache_data(show_spinner=False)
def cached_opn_candidates() -> list[dict]:
    return [candidate.model_dump() for candidate in OpnCandidateCatalog(PATHS).list_candidates()]


@st.cache_data(show_spinner=False, ttl=30)
def cached_opn_rankings() -> list[dict]:
    return [ranking.model_dump() for ranking in OpnCandidateCatalog(PATHS).rank_candidates()]


@st.cache_data(show_spinner=False, ttl=30)
def cached_opn_construct_designs() -> list[dict]:
    return [design.model_dump() for design in OpnCandidateCatalog(PATHS).construct_designs()]


@st.cache_data(show_spinner=False)
def cached_signal_peptide_library() -> list[dict]:
    return signal_peptide_library_service().library_rows()


@st.cache_data(show_spinner=False, ttl=3600)
def cached_uniprot_signal_peptides(taxon_id: int, size: int, reviewed_only: bool) -> dict:
    result = SignalPeptideScreeningService(PATHS).discover_and_persist_uniprot_candidates(
        taxon_id=taxon_id,
        max_records=size,
        reviewed_only=reviewed_only,
        exclude_existing=True,
    )
    return {
        "rows": result.rows,
        "source_url": result.source_url,
        "errors": result.errors,
        "duplicate_count": result.duplicate_count,
        "duplicate_rows": result.duplicate_rows,
    }


@st.cache_data(show_spinner=False, ttl=3600)
def cached_opn_cds_designs(candidate_ids: tuple[str, ...], per_construct: int, seed: int) -> dict:
    result = CdsDesignService(PATHS).design_opn_shortlist(
        list(candidate_ids),
        cds_candidates_per_construct=per_construct,
        seed=seed,
    )
    return result.model_dump()


def signal_peptide_library_service() -> SignalPeptideLibraryService:
    source = OpnSignalPeptideCandidateSource(OpnCandidateCatalog(PATHS))
    return SignalPeptideLibraryService(source.list_candidates())


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


def sidebar_navigation() -> str:
    st.sidebar.title("演示导航")
    page = st.sidebar.radio(
        "选择功能",
        ["项目总览", "结果浏览", "OPN 信号肽", "仿真验证", "运行日志"],
        index=2,
    )
    st.sidebar.divider()
    st.sidebar.markdown(
        """
        **推荐演示顺序**

        1. OPN 信号肽
        2. 结果浏览
        3. 运行日志
        """
    )
    if page == "OPN 信号肽":
        st.sidebar.info("当前建议：首轮小试做 PAS_chr3_0030、DDDK18，并保留 alpha-factor 作为对照。")
    elif page == "结果浏览":
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

