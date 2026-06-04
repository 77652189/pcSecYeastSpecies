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
from app.services.health import HealthService
from app.services.logs import RunLogService
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


def rename_columns(frame: pd.DataFrame, labels: dict) -> pd.DataFrame:
    return frame.rename(columns={key: value for key, value in labels.items() if key in frame.columns})


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
    tab_overview, tab_results, tab_simulation, tab_logs = st.tabs(["项目总览", "结果浏览", "仿真验证", "运行日志"])
    with tab_overview:
        render_overview()
    with tab_results:
        render_results_browser()
    with tab_simulation:
        render_simulation()
    with tab_logs:
        render_logs()


if __name__ == "__main__":
    main()
