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


# ---------------------------------------------------------------------------
# 仿真结果页（render_pichia_results）本地化字典 — 单一集中来源。
# 结果页里大量表格是直接把引擎 payload 的英文字段名/枚举值倒出来的；这里统一翻译，
# 避免以后改代码时翻译散落、不一致或漏改。新增 payload 字段时只在这里加一行即可。
# ---------------------------------------------------------------------------

# 结果页各 dataframe 里出现的字段名（列名）。缺失的键回退显示原文。
SIMULATION_RESULT_COLUMN_LABELS = {
    # LP 归因：OE 可缓解 / 下界 floor / bound marginals
    "reaction": "反应",
    "reaction_id": "反应",
    "variable_id": "变量 ID",
    "variable_index_0based": "变量下标",
    "bound_type": "边界类型",
    "secretory_process": "分泌资源层",
    "marginal": "影子价格",
    "abs_marginal": "影子价格绝对值",
    "oe_actionable": "OE 可缓解",
    "flux": "通量",
    # 约束块 / 约束级 marginals
    "constraint_type": "约束类型",
    "block": "约束块",
    "row_count": "行数",
    "nonzero_marginal_count": "非零影子价格数",
    "sum_abs_marginal": "影子价格绝对值之和",
    "max_abs_marginal": "最大影子价格绝对值",
    "row_index_0based": "行下标",
    "row_index_1based": "行号",
    # 求解器稳健性 per_method
    "method": "求解算法",
    "result_status": "状态",
    "top_oe_actionable_reaction_id": "首位 OE 可缓解反应",
    "top_dominant_block": "首位主导约束块",
    # 成本斜率
    "mu": "生长速率 μ",
    "cost_key": "成本类型",
    "slope": "斜率",
    "point_count": "点数",
    "status": "状态",
    "target_exchange_ratio": "目标分泌比例",
    "objective_reaction": "目标反应",
    "objective_value": "目标值",
    "target_exchange_reaction": "目标分泌反应",
    "glucose_flux": "葡萄糖通量",
    "glucose_cost": "葡萄糖成本",
    "glucose_cost_status": "葡萄糖成本状态",
    "ribosome_reaction": "核糖体反应",
    "ribosome_flux": "核糖体通量",
    "ribosome_cost": "核糖体成本",
    "medium_compatibility_mode": "培养基兼容模式",
    "message": "求解信息",
    "success": "求解成功",
    # 生长权衡
    "secretion_flux": "分泌通量",
    "secretion_per_biomass": "单位生物量分泌",
    "interpretation": "解读",
}

# 结果页里出现的枚举 / 编码值。缺失的值回退显示原文。集中收编了此前散落在
# simulation_results.py 的形状 / 资源层 / 歧义类型标签。
SIMULATION_RESULT_VALUE_LABELS = {
    # 求解器 / 排序稳健性分类
    "ranking-insensitive-to-solver": "跨求解器稳定（可信）",
    "ranking-sensitive-to-solver": "跨求解器翻转（数值假象，非结论）",
    "ranking-insensitive-to-capacity": "跨容量假设稳定（可信）",
    "ranking-sensitive-to-capacity": "跨容量假设翻转（不可信）",
    "inconclusive": "不足以判定",
    # 边界类型
    "upper": "上限（OE 可放宽）",
    "lower": "下限（最低要求，OE 动不了）",
    # 分泌资源层 / secretory_process
    "ribosome": "翻译（核糖体）",
    "proteasome_degradation": "蛋白降解（蛋白酶体）",
    "disulfide_folding": "二硫键折叠 / DSB",
    "chaperone_folding": "分子伴侣折叠",
    "n_glycan_processing": "N-糖基化",
    "o_glycan_processing": "O-糖基化",
    "glycosylation": "糖基化",
    "erad_misfolding": "错误折叠 / ERAD",
    "misfolding_erad": "错误折叠 / ERAD",
    "er_translocation": "ER 转运",
    "er_to_golgi_transport": "ER→Golgi 转运",
    "golgi_surface_transport": "Golgi→胞外运输",
    "secretory_capacity": "分泌容量",
    "metabolic_or_other": "代谢 / 其它",
    "target_secretory_reaction": "目标蛋白分泌反应",
    "target_exchange": "目标蛋白分泌交换",
    "target_related": "目标蛋白相关",
    "growth": "生长",
    "medium_exchange": "培养基交换",
    "unknown": "未解析",
    # OE 剂量响应形状
    "saturating": "饱和型（适度过表达就够，再加收益递减）",
    "linear": "线性型（还在涨，值得进一步加大表达）",
    "threshold": "阈值型（要超过某个最小倍数才起效）",
    "flat_no_response": "无响应（任何倍数都几乎没提升，别过表达）",
    "non_monotonic_numerical_artifact": "非单调（数值假象，不可作结论）",
    "insufficient_points": "数据点不足",
    # R4 歧义类型
    "near_tie": "近似并列（模型分不清谁更好）",
    "capacity_flip": "跨容量假设翻转",
    "solver_flip": "跨求解器翻转",
    # result_status / 状态码
    "draft_lp_sensitivity": "LP 敏感度（草稿）",
    "draft_lp_sensitivity_unavailable": "LP 敏感度不可用",
    "draft_solver_robustness": "求解器稳健性（草稿）",
    "draft_oe_dose_response": "OE 剂量响应（草稿）",
    "draft_oe_dose_response_insufficient": "OE 剂量响应：数据点不足",
    "draft_oe_dose_response_unavailable": "OE 剂量响应不可用",
    "draft_value_of_information": "价值-of-information（草稿）",
    "draft_ranking_robustness": "排序稳健性（草稿）",
    "draft_cost_slope_analysis": "蛋白成本分析（草稿）",
    "draft_matlab_compatible_cost_slope": "MATLAB 兼容成本斜率（草稿）",
    "draft_cost_slope_unavailable": "成本斜率不可用",
    # 成本斜率 / 通量状态
    "uptake_flux": "摄取通量",
    "non_uptake_flux": "非摄取通量",
    "slope_estimated": "已估算斜率",
    "zero_variance": "无变化（方差为零）",
    "glucose_cost": "葡萄糖成本",
    "ribosome_cost": "核糖体成本",
    # 约束块名
    "stoichiometric": "化学计量（代谢质量平衡）",
    "secretory_coupling": "分泌耦合",
    "metabolic_coupling": "代谢耦合",
    "protein_mass": "蛋白质量预算",
    "proteasome": "蛋白酶体",
    "ribosome_assembly": "核糖体装配",
    "ribosome_translation": "核糖体翻译",
    "misfolding": "错误折叠",
    "mitochondrial": "线粒体",
    # 约束类型
    "eq": "等式约束",
    "ub": "上界约束",
    "lb": "下界约束",
    # 布尔
    True: "是",
    False: "否",
    "True": "是",
    "False": "否",
}


def sim_result_column_label(name: object) -> str:
    """结果页 dataframe 列名 -> 中文（未知列回退原文）。"""
    return SIMULATION_RESULT_COLUMN_LABELS.get(str(name), str(name))


def sim_result_value_label(value: object) -> str:
    """结果页枚举/编码值 -> 中文（未知值回退原文；None -> —）。"""
    if value is None:
        return "—"
    if value in SIMULATION_RESULT_VALUE_LABELS:
        return SIMULATION_RESULT_VALUE_LABELS[value]
    return SIMULATION_RESULT_VALUE_LABELS.get(str(value), str(value))


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
