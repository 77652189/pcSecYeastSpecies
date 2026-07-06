"""Generates a natural-language report + recommendations from the structured
genome-wide screen dimensional results. The LLM client is a pluggable port:
swap in a different provider by implementing ReportGenerator, without
touching the UI or the analysis layer.

Provider is selected via the REPORT_LLM_PROVIDER env var (default "openai").
API keys are read from environment variables, never hardcoded or passed on
the command line.
"""

from __future__ import annotations

import json
import os
from typing import Protocol

REPORT_SYSTEM_PROMPT = (
    "你是毕赤酵母代谢工程方向的研究助理，正在帮生物研发同事解读一次全基因组 KO/OE 分泌-生长权衡筛查的结果。"
    "读者是熟悉分子生物学但不一定读过原始数据的同事和研发组长。请用中文，先给结论，再给证据，再给可执行的下一步建议。"
    "不要编造数据里没有的数字；如果某个维度是空的，如实说明并解释可能原因，不要回避。"
)

DEFAULT_MODEL = "gpt-4o-mini"


class ReportGenerator(Protocol):
    def generate(self, dimensional_summaries: list[dict[str, object]], run_metadata: dict[str, object]) -> str: ...


class OpenAIReportGenerator:
    """Default adapter. Requires OPENAI_API_KEY in the environment."""

    def __init__(self, model: str = DEFAULT_MODEL) -> None:
        self.model = model

    def generate(self, dimensional_summaries: list[dict[str, object]], run_metadata: dict[str, object]) -> str:
        from openai import OpenAI

        api_key = os.environ.get("OPENAI_API_KEY")
        if not api_key:
            raise RuntimeError(
                "OPENAI_API_KEY 未设置。请在环境变量或 .env 文件里配置后再生成报告；数据本身已经算好，随时可以重新生成报告。"
            )
        client = OpenAI(api_key=api_key)
        user_prompt = _build_user_prompt(dimensional_summaries, run_metadata)
        response = client.chat.completions.create(
            model=self.model,
            messages=[
                {"role": "system", "content": REPORT_SYSTEM_PROMPT},
                {"role": "user", "content": user_prompt},
            ],
        )
        return response.choices[0].message.content or ""


def _build_user_prompt(dimensional_summaries: list[dict[str, object]], run_metadata: dict[str, object]) -> str:
    return (
        f"本次筛查的运行信息：\n{json.dumps(run_metadata, ensure_ascii=False, indent=2)}\n\n"
        f"各靶点的维度分析结果（每个维度只截取了前若干条，完整数据见原始CSV）：\n"
        f"{json.dumps(dimensional_summaries, ensure_ascii=False, indent=2)}\n\n"
        "请输出一份结构化报告，包含：\n"
        "1. 概览（跑了多少基因、多少靶点、整体发现了什么）\n"
        "2. 每个维度的关键发现（必需基因、产量升高但生长受损的候选、零代价候选、OE候选、靶点特异性差异——如果某个维度是空的要明确说明）\n"
        "3. 最值得优先跟进的3-5个候选基因，说明推荐理由\n"
        "4. 对下一步实验/建模工作的具体建议"
    )


def get_default_generator() -> ReportGenerator:
    provider = os.environ.get("REPORT_LLM_PROVIDER", "openai").strip().lower()
    if provider == "openai":
        return OpenAIReportGenerator()
    raise ValueError(f"Unknown REPORT_LLM_PROVIDER: {provider!r}. Supported: 'openai'.")


__all__ = [
    "DEFAULT_MODEL",
    "REPORT_SYSTEM_PROMPT",
    "OpenAIReportGenerator",
    "ReportGenerator",
    "get_default_generator",
]
