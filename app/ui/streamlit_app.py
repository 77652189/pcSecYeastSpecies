from __future__ import annotations

from pathlib import Path

import pandas as pd
import streamlit as st

from app.adapters.matlab import MatlabAdapter
from app.adapters.powershell import PowerShellAdapter
from app.core.i18n import (
    DATASET_COLUMN_LABELS,
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
from app.services.logs import RunLogService
from app.services.opn import (
    DEFAULT_OPN_CANDIDATE,
    DEFAULT_OPN_PRODUCTION_RATIO,
    OPN_SHORTLIST,
    OpnCandidateCatalog,
    OpnSimulationService,
    opn_category_label,
)
from app.services.results import PlotBuilder, ResultCatalog, ResultLoader
from app.services.simulation import SimulationService


st.set_page_config(
    page_title="pcSecYeastSpecies",
    page_icon=None,
    layout="wide",
    initial_sidebar_state="expanded",
)


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


@st.cache_data(show_spinner=False, ttl=3600)
def cached_opn_cds_designs(candidate_ids: tuple[str, ...], per_construct: int, seed: int) -> dict:
    result = CdsDesignService(PATHS).design_opn_shortlist(
        list(candidate_ids),
        cds_candidates_per_construct=per_construct,
        seed=seed,
    )
    return result.model_dump()


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


def render_overview() -> None:
    datasets = dataset_frame()
    health = cached_health()
    items = pd.DataFrame(health["items"])

    st.markdown(
        """
        <div class="concept-box">
        这个工具用于查看和验证酵母蛋白分泌模型的计算结果。你可以把它理解成一个“细胞工厂模拟器”：
        选择物种、条件或结果文件后，观察模型预测的生长、代谢通量、蛋白成本和求解状态。
        </div>
        """,
        unsafe_allow_html=True,
    )

    col1, col2, col3, col4 = st.columns(4)
    col1.metric("可浏览结果文件", len(datasets))
    col2.metric("结果主题数量", datasets["category_label"].nunique() if not datasets.empty else 0)
    col3.metric("模型物种数量", 3)
    col4.metric("最近输出文件", len(RunLogService(PATHS).recent_files()))

    with st.expander("这个网站能做什么、输入是什么、输出是什么", expanded=True):
        st.markdown(
            """
            **能做什么：** 比较三种酵母在蛋白分泌、温度变化、碳源变化、目标蛋白生产等场景下的模型预测结果。

            **输入是什么：** 当前版本主要读取项目中已经计算好的 `Results/` 文件；仿真验证页额外允许输入一个生长速率 `mu`。

            **输出是什么：** 表格、趋势图、求解器 SoPlex 输出、线性规划 LP 文件，以及最近一次仿真的目标函数值（objective value）。
            """
        )

    st.subheader("部署状态")
    if not items.empty:
        items["status_label"] = items["status"].map(status_label)
        status_order = {"ok": 0, "warning": 1, "missing": 2, "error": 3}
        items["sort"] = items["status"].map(status_order).fillna(9)
        display = items.sort_values(["sort", "name"])[["name", "status_label", "detail"]]
        st.dataframe(rename_columns(display, HEALTH_COLUMN_LABELS), use_container_width=True, hide_index=True)
    with st.expander("怎么看部署状态"):
        st.markdown(
            """
            - **正常：** 依赖、模型文件或结果目录已找到。
            - **缺失/错误：** 相关依赖没有安装、路径不对，或预检脚本运行失败。
            - Windows 原生 `soplex` 可以缺失；本项目当前通过 Docker 里的 SoPlex 求解。
            """
        )

    st.subheader("模型覆盖范围")
    st.dataframe(
        pd.DataFrame(
            [
                {"物种": "酿酒酵母（S. cerevisiae）", "代码": "SCE", "模型": "pcSecYeast", "当前功能": "结果浏览 + 已验证小规模仿真"},
                {"物种": "毕赤酵母（K. phaffii）", "代码": "PPA", "模型": "pcSecPichia", "当前功能": "结果浏览"},
                {"物种": "马克斯克鲁维酵母（K. marxianus）", "代码": "KMX", "模型": "pcSecKmarx", "当前功能": "结果浏览"},
            ]
        ),
        use_container_width=True,
        hide_index=True,
    )

    if PATHS.phase1_png.exists():
        st.subheader("最近一次图形验证")
        st.image(str(PATHS.phase1_png), width=620)


def render_results_browser() -> None:
    datasets = dataset_frame()
    st.markdown(
        """
        <div class="concept-box">
        在这里浏览项目已经计算好的结果文件。建议先按“结果主题”和“物种”筛选，再打开一个数据集查看表格和图。
        </div>
        """,
        unsafe_allow_html=True,
    )
    if datasets.empty:
        st.warning("没有找到 Results 目录中的 Excel 或 MATLAB 结果文件。")
        return

    with st.sidebar:
        st.subheader("结果筛选")
        species_options = ["全部", *sorted(datasets["species_label"].unique())]
        category_options = ["全部", *sorted(datasets["category_label"].unique())]
        species = st.selectbox("物种", species_options, index=0)
        category = st.selectbox("结果主题", category_options, index=0)
        keyword = st.text_input("关键词", placeholder="例如 Insulin、CSource、SCE")

    filtered = datasets.copy()
    if species != "全部":
        filtered = filtered[filtered["species_label"] == species]
    if category != "全部":
        filtered = filtered[filtered["category_label"] == category]
    if keyword:
        mask = filtered["name"].str.contains(keyword, case=False, na=False) | filtered["id"].str.contains(keyword, case=False, na=False)
        filtered = filtered[mask]

    with st.expander("如何选择数据集", expanded=True):
        st.markdown(
            """
            - **结果主题** 对应原项目 `Results/` 下的文件夹，例如“碳源分析”“目标蛋白生长分析”。
            - **物种** 是从文件名推断的；如果显示“未知物种”，通常说明文件名没有明确物种缩写。
            - **Excel 结果表** 更适合直接阅读；**MATLAB 结果文件** 通常包含矩阵或结构体，页面会先显示变量摘要。
            """
        )

    st.subheader("结果数据集")
    display_columns = ["name", "category_label", "species_label", "suffix", "size_kb", "modified_at", "id"]
    st.dataframe(
        rename_columns(filtered[display_columns], DATASET_COLUMN_LABELS),
        use_container_width=True,
        hide_index=True,
    )
    if filtered.empty:
        st.info("当前筛选条件下没有结果文件。")
        return

    selected_id = st.selectbox(
        "打开数据集",
        filtered["id"].tolist(),
        format_func=lambda value: filtered.loc[filtered["id"] == value, "name"].iloc[0],
    )
    loaded_data = cached_loaded_dataset(selected_id)
    loaded = ResultLoader().load_dataset(ResultCatalog(PATHS).get_dataset(selected_id))

    st.subheader(loaded.info.name)
    st.caption(f"{category_label(loaded.info.category)} / {species_label(loaded.info.species)} / {short_path(selected_id)}")

    if loaded.variable_summary:
        st.markdown("**MATLAB 变量摘要**")
        st.dataframe(pd.DataFrame(loaded.variable_summary), use_container_width=True, hide_index=True)

    table_names = list(loaded.tables)
    if table_names:
        table_name = st.selectbox("表格或变量", table_names)
        table = pd.DataFrame(loaded_data["tables"][table_name])
        st.dataframe(table, use_container_width=True, hide_index=True)

        with st.expander("图怎么看"):
            st.markdown(
                """
                图表是对当前表格中数值列的快速预览。对于碳源分析，散点图会比较实验测得生长速率和模型预测生长速率；
                对于一般矩阵结果，页面会展示前几个数值变量的趋势，帮助你判断数据是否可读、是否需要下载原始文件进一步分析。
                """
            )
        figure = PlotBuilder().build_chart(loaded, table_name)
        if figure is not None:
            st.plotly_chart(figure, use_container_width=True)
        else:
            st.info("这个数据集没有适合自动绘图的数值列，建议下载原始文件进一步查看。")

    file_path = Path(loaded_data["info"]["path"])
    if file_path.exists():
        st.caption("下载的是原始结果文件，不是重新计算生成的报告。")
        st.download_button(
            "下载原始结果文件",
            file_path.read_bytes(),
            file_name=file_path.name,
            mime="application/octet-stream",
        )


def render_simulation() -> None:
    st.markdown(
        """
        <div class="concept-box">
        这里运行一个已验证的小规模示例计算：酿酒酵母（SCE）在葡萄糖（glucose）条件下，给定生长速率 mu，
        生成线性规划 LP 文件并调用 Docker 中的 SoPlex 求解器。
        </div>
        """,
        unsafe_allow_html=True,
    )
    left, right = st.columns([1, 2])
    with left:
        mu = st.number_input("生长速率 mu（h^-1）", min_value=0.01, max_value=0.44, value=0.10, step=0.01, format="%.2f")
        timeout = st.number_input("SoPlex 超时秒数", min_value=120, max_value=3600, value=300, step=60)
        st.caption("点击运行后会调用 MATLAB 生成 LP，再调用 Docker SoPlex 求解，可能需要数分钟。")
        run_clicked = st.button("运行示例仿真", type="primary")

    with st.expander("这些输入和输出是什么意思", expanded=True):
        st.markdown(
            """
            - **生长速率 mu：** 模型中固定的细胞生长速度，单位是 h^-1。
            - **线性规划 LP：** MATLAB 根据模型和条件生成的优化问题文件。
            - **SoPlex：** 用来求解 LP 的优化求解器。
            - **objective value：** 求解器返回的目标函数值，用来判断这次优化结果的数值结果；它不是直接的产量单位。
            """
        )

    if run_clicked:
        service = SimulationService(PATHS, MatlabAdapter(), PowerShellAdapter())
        with st.spinner("正在生成 LP 并调用 Docker SoPlex，请稍等..."):
            result = service.run_sce_glucose_smoke(mu=float(mu), timeout_seconds=int(timeout))
        st.session_state["last_simulation_result"] = result.model_dump()

    result_data = st.session_state.get("last_simulation_result")
    with right:
        if result_data:
            st.success(result_data["message"]) if result_data["success"] else st.error(result_data["message"])
            summary = {
                "生长速率 mu": result_data["mu"],
                "目标函数值（objective value）": result_data["objective_value"],
                "LP 文件": result_data["lp_file"],
                "SoPlex 输出文件": result_data["output_file"],
            }
            st.dataframe(pd.DataFrame([summary]), use_container_width=True, hide_index=True)
            if result_data["command_output"]:
                with st.expander("查看命令输出"):
                    st.code(result_data["command_output"][-12000:], language="text")
        else:
            latest, summary = RunLogService(PATHS).latest_soplex_summary()
            if latest and summary:
                st.info("下面显示最近一次 SoPlex 输出摘要。")
                display = pd.DataFrame(
                    [
                        {
                            "输出文件": str(latest.relative_to(PATHS.repo_root)),
                            "optimal_label": status_label(summary.optimal),
                            "objective": summary.objective_value,
                            "status": summary.status_line,
                        }
                    ]
                )
                st.dataframe(rename_columns(display, SOPLEX_COLUMN_LABELS), use_container_width=True, hide_index=True)


def render_opn_signal_peptides() -> None:
    st.markdown(
        """
        <div class="concept-box">
        这里展示“用于毕赤酵母生产人骨桥蛋白 OPN”的候选信号肽。候选表把每个 leader 接到你提供的成熟 OPN 序列前端，
        然后用 pcSecPichia 生成小规模 LP 验证。当前页面用于演示和筛选优先级，不等同于真实发酵产量证明。
        </div>
        """,
        unsafe_allow_html=True,
    )
    candidates = pd.DataFrame(cached_opn_candidates())
    if candidates.empty:
        st.error("没有找到 OPN 候选表。请先生成 Data/pcSecPichia/TargetProtein_OPN_candidates_meta.csv。")
        return

    candidates["分类"] = candidates["category"].map(opn_category_label)
    candidates["leader长度aa"] = candidates["leader_length"]
    candidates["完整构建长度aa"] = candidates["construct_length"]

    metric1, metric2, metric3, metric4 = st.columns(4)
    metric1.metric("候选数量", len(candidates))
    metric2.metric("成熟 OPN 长度", "298 aa")
    metric3.metric("默认候选", DEFAULT_OPN_CANDIDATE)
    metric4.metric("默认产量约束", f"{DEFAULT_OPN_PRODUCTION_RATIO:.0e}")

    render_opn_recommendation_board(pd.DataFrame(cached_opn_rankings()))

    with st.expander("这个页面怎么看", expanded=True):
        st.markdown(
            """
            - **候选信号肽表**：直接展示每个候选的 leader 序列、长度和加工路线，便于横向比较。
            - **候选详情**：解释为什么选它，以及真实实验中要注意什么。
            - **运行小规模验证**：调用 MATLAB 生成 OPN/Pichia LP，再调用 Docker SoPlex 求解；结果只表示模型链路是否可行。
            """
        )

    display = candidates[
        [
            "candidate_id",
            "分类",
            "leader长度aa",
            "完整构建长度aa",
            "leader_sequence",
            "processing_route",
        ]
    ].rename(
        columns={
            "candidate_id": "候选 ID",
            "leader_sequence": "leader 序列",
            "processing_route": "加工路线",
        }
    )
    st.subheader("候选信号肽表")
    st.dataframe(display, use_container_width=True, hide_index=True)

    render_opn_construct_design_table(pd.DataFrame(cached_opn_construct_designs()))
    render_opn_cds_design_panel()

    default_index = candidates.index[candidates["candidate_id"] == DEFAULT_OPN_CANDIDATE]
    selected_index = int(default_index[0]) if len(default_index) else 0
    selected_id = st.selectbox(
        "选择一个候选查看详情",
        candidates["candidate_id"].tolist(),
        index=selected_index,
    )
    selected = candidates.loc[candidates["candidate_id"] == selected_id].iloc[0]

    left, right = st.columns([1, 1])
    with left:
        st.subheader("候选详情")
        st.markdown(f"**分类：** {opn_category_label(selected['category'])}")
        st.markdown(f"**leader 长度：** {selected['leader_length']} aa")
        st.markdown(f"**信号肽长度：** {len(selected['signal_peptide_sequence'])} aa")
        st.markdown(f"**完整构建长度：** {selected['construct_length']} aa")
        st.markdown(f"**选择理由：** {selected['rationale']}")
        st.markdown(f"**注意事项：** {selected['caution']}")
        with st.expander("查看序列"):
            st.code(selected["leader_sequence"], language="text")
    with right:
        st.subheader("模型会看什么")
        st.markdown(
            """
            模型主要会把 leader/OPN 的序列长度、氨基酸组成、ER 分泌路径、翻译和分泌机器资源约束写入 LP。
            它不会直接预测真实滴度，也不会模拟宿主蛋白酶是否切割 OPN 内部位点。
            """
        )
        st.warning("成熟 OPN 内部检测到 `RR` 和 `KR` 二碱性位点；使用 alpha-factor pro 路线时需要额外关注 Kex2 类异常切割风险。")
        latest, summary = OpnSimulationService(PATHS, MatlabAdapter()).latest_candidate_result(selected_id)
        if latest and summary:
            st.info(
                f"最近一次验证：{status_label(summary.optimal)}。"
                f"目标函数值 {summary.objective_value or '未读取到'}，输出文件 {latest.name}。"
            )

    st.subheader("运行 OPN 小规模验证")
    with st.form("opn_candidate_smoke_form"):
        form_left, form_right = st.columns([1, 1])
        with form_left:
            run_candidate_id = st.selectbox(
                "要运行的候选",
                candidates["candidate_id"].tolist(),
                index=candidates["candidate_id"].tolist().index(selected_id),
                key="opn_run_candidate_id",
            )
            mu = st.number_input("生长速率 mu（h^-1）", min_value=0.01, max_value=0.44, value=0.10, step=0.01, format="%.2f")
        with form_right:
            production_ratio = st.selectbox(
                "OPN 生产通量约束",
                [1e-8, 1e-6],
                index=0,
                format_func=lambda value: f"{value:.0e}",
            )
            timeout = st.number_input("SoPlex 超时秒数", min_value=120, max_value=3600, value=300, step=60)
        st.caption("运行会调用 MATLAB R2020b+ 生成 LP，再调用 Docker 镜像 `pcsec-soplex:24.04` 求解，通常需要 1-3 分钟。")
        run_clicked = st.form_submit_button("运行 OPN 候选验证", type="primary")

    if run_clicked:
        service = OpnSimulationService(PATHS, MatlabAdapter())
        with st.spinner("正在生成 OPN/Pichia LP 并调用 SoPlex，请稍等..."):
            result = service.run_candidate_smoke(
                candidate_id=run_candidate_id,
                mu=float(mu),
                production_ratio=float(production_ratio),
                timeout_seconds=int(timeout),
            )
        st.session_state["last_opn_result"] = result.model_dump()

    result_data = st.session_state.get("last_opn_result")
    if result_data:
        render_opn_result_explanation(result_data)
    else:
        latest, summary = OpnSimulationService(PATHS, MatlabAdapter()).latest_candidate_result(selected_id)
        if latest and summary:
            render_opn_latest_result_explanation(selected_id, latest, summary)

    if PATHS.opn_candidate_csv.exists():
        st.download_button(
            "下载模型候选 CSV",
            PATHS.opn_candidate_csv.read_bytes(),
            file_name=PATHS.opn_candidate_csv.name,
            mime="text/csv",
        )


def render_opn_construct_design_table(designs: pd.DataFrame) -> None:
    st.subheader("实验构建设计输出")
    if designs.empty:
        st.info("暂时没有可导出的 OPN 构建设计。")
        return

    st.markdown(
        """
        这张表是给实验侧看的：把“推荐角色、信号肽序列、成熟 OPN、完整蛋白序列、Kex2 风险、下一步密码子优化”放在一起。
        网页里只展示关键列，下载 CSV 会包含完整序列。
        """
    )
    preview = designs[
        [
            "candidate_id",
            "experimental_role",
            "leader_sequence",
            "signal_peptide_sequence",
            "full_protein_length",
            "contains_alpha_pro_region",
            "kex2_risk",
            "codon_optimization_next",
        ]
    ].rename(
        columns={
            "candidate_id": "候选 ID",
            "experimental_role": "实验角色",
            "leader_sequence": "leader 序列",
            "signal_peptide_sequence": "信号肽序列",
            "full_protein_length": "完整蛋白长度",
            "contains_alpha_pro_region": "含 alpha pro 区",
            "kex2_risk": "Kex2/切割风险",
            "codon_optimization_next": "下一步",
        }
    )
    st.dataframe(preview, use_container_width=True, hide_index=True)

    export = designs.rename(
        columns={
            "candidate_id": "候选ID",
            "experimental_role": "实验角色",
            "recommendation": "推荐级别",
            "leader_sequence": "leader序列",
            "signal_peptide_sequence": "信号肽序列",
            "mature_opn_sequence": "成熟OPN序列",
            "full_protein_sequence": "完整蛋白序列",
            "leader_length": "leader长度",
            "signal_peptide_length": "信号肽长度",
            "mature_opn_length": "成熟OPN长度",
            "full_protein_length": "完整蛋白长度",
            "contains_alpha_pro_region": "是否含alpha_pro区",
            "processing_route": "加工路线",
            "kex2_risk": "Kex2风险",
            "codon_optimization_next": "下一步密码子优化",
            "note": "备注",
        }
    )
    st.download_button(
        "下载实验构建设计 CSV",
        export.to_csv(index=False).encode("utf-8-sig"),
        file_name="OPN_Pichia_signal_peptide_construct_design.csv",
        mime="text/csv",
    )


def render_opn_cds_design_panel() -> None:
    st.subheader("首轮实验构建方案")
    st.markdown(
        """
        建议首轮做 `OPN_PPA_PASCHR3_0030`、`OPN_PPA_DDDK18`，并保留
        `OPN_ALPHA_FULL_PROJECT` 作为 alpha-factor 对照。这里会调用 PichiaCLM 为每个构建生成
        毕赤酵母 CDS 候选，并导出 CSV、XLSX 和 FASTA，供实验同事讨论。
        """
    )
    st.info("注意：这些 CDS 和模型分数是表达设计参考，不代表真实发酵产量；最终仍需要毕赤酵母小试验证。")
    with st.expander("这两个项目现在怎么通信？", expanded=True):
        st.markdown(
            """
            - **当前做法：函数级调用。** pcSecYeastSpecies 的 Streamlit 页面调用 `CdsDesignService`，服务层再调用 `PichiaClmAdapter`，适配器直接 import 本机的 PichiaCLM 核心模型。
            - **不推荐：Streamlit 调 Streamlit。** 两个网页都偏展示层，互相 HTTP 调用会让错误处理、状态管理和未来上线都变复杂。
            - **后期上线：服务级调用。** PichiaCLM 可以单独提供 FastAPI；pcSecYeastSpecies 只把适配器从“本地 import”换成“HTTP 请求”，前端和业务服务不需要重写。
            """
        )

    candidates = pd.DataFrame(cached_opn_candidates())
    default_ids = [candidate_id for candidate_id in OPN_SHORTLIST if candidate_id in set(candidates["candidate_id"])]
    selected_ids = st.multiselect(
        "选择要生成 CDS 的 OPN 构建",
        candidates["candidate_id"].tolist(),
        default=default_ids,
    )
    left, right = st.columns([1, 1])
    with left:
        per_construct = st.number_input("每个构建生成几个 CDS 候选", min_value=1, max_value=10, value=3, step=1)
    with right:
        seed = st.number_input("随机种子", min_value=1, max_value=9999, value=42, step=1)

    if not selected_ids:
        st.info("请至少选择一个 OPN 构建。")
        return
    if not st.button("生成首轮构建方案", type="primary"):
        st.caption("为了避免页面打开时自动加载深度学习模型，点击按钮后才会调用 PichiaCLM。")
        return

    with st.spinner("正在调用 PichiaCLM 生成 CDS，首次加载模型可能需要几十秒..."):
        result = cached_opn_cds_designs(tuple(selected_ids), int(per_construct), int(seed))

    if not result["available"]:
        st.error(result["message"])
        return
    st.success(result["message"])
    records = pd.DataFrame(result["records"])
    if records.empty:
        return

    matched_count = int(records["translation_matches_input"].sum())
    stop_count = int(records["internal_stop_codons"].map(len).sum())
    summary_col1, summary_col2, summary_col3, summary_col4 = st.columns(4)
    summary_col1.metric("CDS 候选记录", len(records))
    summary_col2.metric("回译匹配蛋白", f"{matched_count}/{len(records)}")
    summary_col3.metric("内部终止密码子", stop_count)
    summary_col4.metric("导出文件", "CSV / XLSX / FASTA")

    preview = records[
        [
            "construct_id",
            "experimental_role",
            "recommendation",
            "cds_candidate_rank",
            "recommended_subset",
            "aa_length",
            "cds_length",
            "gc_percent",
            "gc_status",
            "cai_public",
            "quality_status",
            "warnings",
            "restriction_sites",
            "motif_hits",
            "length_multiple_of_three",
            "translation_matches_input",
            "internal_stop_codons",
        ]
    ].rename(
        columns={
            "construct_id": "OPN 构建 ID",
            "experimental_role": "实验角色",
            "recommendation": "推荐级别",
            "cds_candidate_rank": "CDS 候选序号",
            "recommended_subset": "推荐保留",
            "aa_length": "氨基酸长度",
            "cds_length": "CDS 长度 bp",
            "gc_percent": "GC%",
            "gc_status": "GC 状态",
            "cai_public": "CAI（公开毕赤参考）",
            "quality_status": "质控状态",
            "warnings": "提醒数",
            "restriction_sites": "默认/自定义酶切位点数",
            "motif_hits": "不期望 motif 数",
            "length_multiple_of_three": "长度为3倍数",
            "translation_matches_input": "回译匹配蛋白",
            "internal_stop_codons": "内部终止密码子位置",
        }
    )
    st.dataframe(preview, use_container_width=True, hide_index=True)

    selected_record = st.selectbox(
        "查看一条 CDS 序列",
        records.index.tolist(),
        format_func=lambda idx: f"{records.loc[idx, 'construct_id']} / CDS {records.loc[idx, 'cds_candidate_rank']}",
    )
    st.code(records.loc[selected_record, "cds"], language="text")

    st.markdown("**下载给实验同事讨论的文件**")
    col_csv, col_xlsx, col_fasta = st.columns(3)
    download_specs = [
        (col_csv, result.get("csv_file"), "下载 CSV 构建表", "text/csv"),
        (col_xlsx, result.get("xlsx_file"), "下载 XLSX 构建表", "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"),
        (col_fasta, result.get("fasta_file"), "下载 FASTA 序列", "text/plain"),
    ]
    for column, path_value, label, mime in download_specs:
        if not path_value:
            column.warning("文件未生成")
            continue
        path = Path(path_value)
        if path.exists():
            column.download_button(label, path.read_bytes(), file_name=path.name, mime=mime)
        else:
            column.warning(f"未找到 {path.name}")


def render_opn_recommendation_board(rankings: pd.DataFrame) -> None:
    st.subheader("当前建议先做哪几个")
    if rankings.empty:
        st.info("还没有可用于排序的候选结果。请先运行候选验证。")
        return

    shortlist = rankings[rankings["candidate_id"].isin(OPN_SHORTLIST)].sort_values("rank")
    best = shortlist.iloc[0] if not shortlist.empty else rankings.sort_values("rank").iloc[0]
    st.success(
        f"首选：`{best['candidate_id']}`。建议首轮小试做 3 个构建："
        "`OPN_PPA_PASCHR3_0030`、`OPN_PPA_DDDK18`，再加 `OPN_ALPHA_FULL_PROJECT` 作为 alpha-factor 对照。"
    )
    st.markdown(
        """
        **为什么这样选：** pcSec 模型的 objective 差异很小，只靠模型数字不能定最终信号肽。
        所以首轮策略应优先选择 Pichia 来源、避开 alpha pro/Kex2 加工风险的短信号肽，同时保留工业上常用的 alpha-factor 基线做对照。
        """
    )

    display = rankings.sort_values("rank")[
        [
            "rank",
            "candidate_id",
            "experimental_role",
            "recommendation",
            "objective_text",
            "model_rank",
            "objective_delta_percent",
            "risk_level",
            "reason",
        ]
    ].rename(
        columns={
            "rank": "推荐顺序",
            "candidate_id": "候选 ID",
            "experimental_role": "实验角色",
            "recommendation": "建议",
            "objective_text": "模型目标函数值",
            "model_rank": "模型成本排名",
            "objective_delta_percent": "相对最优差距%",
            "risk_level": "主要风险",
            "reason": "理由",
        }
    )
    st.dataframe(display, use_container_width=True, hide_index=True)

    with st.expander("为什么不是只选 objective 最好的"):
        st.markdown(
            """
            当前 LP 的 objective 主要反映模型在固定生长和固定 OPN 生产通量下的资源/底物需求，不是分泌滴度。
            对信号肽而言，真实表达还强烈受切割位点、宿主蛋白酶、糖基化、mRNA/翻译效率和培养工艺影响。
            因此这里把 objective 当作“模型成本参考”，而不是唯一排序依据。
            """
        )


def render_opn_latest_result_explanation(candidate_id: str, output_file: Path, summary) -> None:
    st.subheader("最近一次运行结果怎么读")
    if summary.optimal:
        st.success(f"结论：`{candidate_id}` 最近一次本地验证已经求解成功，可以作为演示用结果。")
    else:
        st.warning(f"结论：`{candidate_id}` 最近一次本地验证没有显示 optimal，需要重新运行或查看输出。")

    col1, col2, col3 = st.columns(3)
    col1.metric("求解状态", status_label(summary.optimal))
    col2.metric("目标函数值", summary.objective_value or "未读取到")
    col3.metric("候选", candidate_id)

    st.markdown(
        f"""
        **它说明了什么：**
        这个候选在最近一次本地 OPN/Pichia 小规模验证中，SoPlex 能够完成求解。

        **它不说明什么：**
        这不是实际发酵滴度，也不是最终推荐排序。它只是说明这个候选进入模型比较是可行的。

        **目标函数值怎么看：**
        `{summary.objective_value or '未读取到'}` 是模型优化目标的数值结果，不是产量单位。它要和其他候选在相同参数下横向比较才有意义。
        """
    )

    st.markdown("**最值得看的输出文件**")
    st.dataframe(
        pd.DataFrame(
            [
                {
                    "文件": "SoPlex 输出文件",
                    "用途": "确认是否 optimal，并查看 objective value。",
                    "路径": compact_path(output_file),
                }
            ]
        ),
        use_container_width=True,
        hide_index=True,
    )


def render_opn_result_explanation(result_data: dict) -> None:
    success = bool(result_data.get("success"))
    candidate_id = result_data.get("candidate_id", "未知候选")
    objective = result_data.get("objective_value") or "未读取到"
    lp_file = compact_path(result_data.get("lp_file"))
    output_file = compact_path(result_data.get("output_file"))

    st.subheader("这次运行结果怎么读")
    if success:
        st.success(f"结论：`{candidate_id}` 在这组模型约束下已经跑通，SoPlex 求解成功。")
    else:
        st.error(f"结论：`{candidate_id}` 这次没有得到可用解，需要查看错误输出或换参数。")

    col1, col2, col3 = st.columns(3)
    col1.metric("求解状态", "成功" if success else "未通过")
    col2.metric("目标函数值", objective)
    col3.metric("候选", str(candidate_id))

    st.markdown(
        f"""
        **它说明了什么：**
        在固定生长速率 `mu={result_data.get('mu')}`、固定 OPN 生产通量 `{result_data.get('production_ratio')}`、
        培养基类型 `{result_data.get('media_type')}` 的条件下，模型能为 `{candidate_id}` 找到一个满足约束的解。

        **它不说明什么：**
        这个结果不是实际发酵产量，也不是说这个信号肽一定最高产。它只是说明“这条候选在当前模型条件下可计算、可比较”。

        **目标函数值怎么看：**
        当前 LP 的目标和葡萄糖交换通量有关。在模型符号约定里，葡萄糖摄取常显示为负数，所以 `-1.074...` 这类值不是产量单位。
        它最适合在同一套参数下和其他候选横向比较，单独看一个数意义有限。
        """
    )

    st.info("演示时可以这样说：这一步是在筛掉模型上不可行的候选，并为后续实验优先级排序提供参考。真正的表达量还需要小试实验验证。")

    st.markdown("**输出文件在哪里**")
    st.dataframe(
        pd.DataFrame(
            [
                {
                    "文件": "LP 输入文件",
                    "用途": "MATLAB 生成的优化问题，给求解器使用，通常不需要人工阅读。",
                    "路径": lp_file,
                },
                {
                    "文件": "SoPlex 输出文件",
                    "用途": "最值得看，里面包含 optimal 和 objective value。",
                    "路径": output_file,
                },
            ]
        ),
        use_container_width=True,
        hide_index=True,
    )

    with st.expander("给研发/实验同事的简短解释"):
        st.markdown(
            f"""
            本次计算验证的是候选信号肽 `{candidate_id}` 接到成熟 OPN 后，在 pcSecPichia 模型中是否能满足
            生长、分泌和资源约束。求解成功说明它可以进入下一轮候选比较；但模型没有模拟蛋白酶切割、糖基化异质性和真实发酵滴度，
            因此不能直接当作表达量结论。
            """
        )

    with st.expander("查看原始命令输出（开发者排错用）"):
        output = result_data.get("command_output", "")
        st.code(output[-12000:] if output else "这次运行没有命令输出。", language="text")


def render_logs() -> None:
    st.markdown(
        """
        <div class="concept-box">
        这里汇总本地运行产生的文件。`.lp` 是优化问题输入，`.lp.out` 是 SoPlex 求解输出，`.png` 是图形验证结果。
        </div>
        """,
        unsafe_allow_html=True,
    )
    service = RunLogService(PATHS)
    recent = pd.DataFrame(service.recent_files())
    st.subheader("最近运行文件")
    with st.expander("哪些文件最值得看", expanded=True):
        st.markdown(
            """
            - **`.lp.out`：** 最重要，包含 SoPlex 是否求解成功以及 objective value。
            - **`.lp`：** MATLAB 生成的线性规划输入文件，通常很大，不建议手工阅读。
            - **`.png`：** 第一阶段图形验证输出，可用于确认绘图链路正常。
            """
        )
    if recent.empty:
        st.info("local_runs 目录暂无文件。")
    else:
        recent["修改时间"] = pd.to_datetime(recent["修改时间"], unit="s").dt.strftime("%Y-%m-%d %H:%M:%S")
        st.dataframe(rename_columns(recent, RUN_FILE_COLUMN_LABELS), use_container_width=True, hide_index=True)

    latest, summary = service.latest_soplex_summary()
    if latest and summary:
        st.subheader("最近 SoPlex 求解摘要")
        display = pd.DataFrame(
            [
                {
                    "文件": str(latest.relative_to(PATHS.repo_root)),
                    "optimal_label": status_label(summary.optimal),
                    "objective": summary.objective_value,
                    "status": summary.status_line,
                }
            ]
        )
        st.dataframe(rename_columns(display, SOPLEX_COLUMN_LABELS), use_container_width=True, hide_index=True)
        with st.expander("查看输出文件末尾日志"):
            text = latest.read_text(encoding="utf-8", errors="replace")
            st.code(text[-12000:], language="text")


def main() -> None:
    app_css()
    page_header()
    page = sidebar_navigation()
    if page == "项目总览":
        render_overview()
    elif page == "结果浏览":
        render_results_browser()
    elif page == "OPN 信号肽":
        render_opn_signal_peptides()
    elif page == "仿真验证":
        render_simulation()
    elif page == "运行日志":
        render_logs()


if __name__ == "__main__":
    main()
