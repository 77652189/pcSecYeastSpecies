from __future__ import annotations

from pathlib import Path


SPECIES_LABELS = {
    "SCE": "酿酒酵母（S. cerevisiae）",
    "PPA": "毕赤酵母（K. phaffii）",
    "KMX": "马克斯克鲁维酵母（K. marxianus）",
    "Unknown": "未知物种",
}

CATEGORY_LABELS = {
    "CSource": "碳源分析",
    "Crabtree": "Crabtree 效应分析",
    "Enzyme_sensitivity_analysis": "酶参数敏感性分析",
    "Experimental_validation": "实验验证",
    "FSEOF": "代谢工程靶点分析（FSEOF）",
    "Growth_rate_HNG": "人源化糖基化生长分析",
    "Growth_rate_TDM": "温度相关生长分析",
    "Growth_rate_TP": "目标蛋白生长分析",
    "Protein_cost_HNG": "人源化糖基化蛋白成本分析",
    "Protein_cost_TDM": "温度相关蛋白成本分析",
    "Protein_cost_TP": "目标蛋白蛋白成本分析",
    "SHAP_Analysis": "机器学习解释分析（SHAP）",
    "Temperature-sensitive_parameters_analysis": "温度敏感参数分析",
    "Results": "综合结果",
}

STATUS_LABELS = {
    "ok": "正常",
    "warning": "提醒",
    "missing": "缺失",
    "error": "错误",
    True: "求解成功",
    False: "未通过",
}

DATASET_COLUMN_LABELS = {
    "name": "数据集名称",
    "category": "结果主题",
    "category_label": "结果主题",
    "species": "物种代码",
    "species_label": "物种",
    "suffix": "文件类型",
    "size_bytes": "文件大小（字节）",
    "size_kb": "文件大小（KB）",
    "modified_at": "修改时间",
    "id": "文件路径",
    "path": "完整路径",
}

HEALTH_COLUMN_LABELS = {
    "name": "检查项",
    "status": "状态",
    "status_label": "状态",
    "detail": "说明",
}

SOPLEX_COLUMN_LABELS = {
    "optimal": "求解状态",
    "optimal_label": "求解状态",
    "objective": "目标函数值（objective value）",
    "status": "求解器状态",
    "输出文件": "输出文件",
    "文件": "文件",
}

RUN_FILE_COLUMN_LABELS = {
    "文件": "文件",
    "大小KB": "大小（KB）",
    "修改时间": "修改时间",
}


def species_label(code: str | None) -> str:
    return SPECIES_LABELS.get(code or "Unknown", code or "未知物种")


def category_label(category: str | None) -> str:
    return CATEGORY_LABELS.get(category or "", category or "未分类结果")


def status_label(status: str | bool | None) -> str:
    if status in STATUS_LABELS:
        return STATUS_LABELS[status]
    return str(status) if status is not None else "未知"


def file_type_label(suffix: str) -> str:
    mapping = {
        ".xlsx": "Excel 结果表",
        ".mat": "MATLAB 结果文件",
        ".lp": "线性规划 LP 文件",
        ".out": "求解器输出",
        ".png": "图片",
        ".sh": "求解脚本",
    }
    return mapping.get(suffix.lower(), suffix or "未知文件")


def short_path(path: str | Path) -> str:
    return str(path).replace("\\", "/")
