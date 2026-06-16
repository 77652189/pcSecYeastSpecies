from __future__ import annotations

from pathlib import Path

import pandas as pd
import streamlit as st

from app.core.i18n import status_label
from app.services.opn import OPN_SHORTLIST
from app.ui.common import compact_path


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


