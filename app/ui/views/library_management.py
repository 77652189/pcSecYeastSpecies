"""序列库与映射管理页。

这里放的都是**偶尔做一次**的管理动作：新增/修改序列库条目、复核基因↔复合体映射。
它们此前挤在每天都走的仿真验证主流程里——研究员每次跑仿真都要从它们中间穿过去。
搬到独立页面后，主流程只剩"选候选 → 跑 → 看结果"。
"""

from __future__ import annotations

import streamlit as st

from app.ui.views.gene_complex_mapping_review import render_gene_complex_mapping_review
from app.ui.views.target_library_manager import render_target_library_manager


def render_library_management() -> None:
    st.subheader("序列库与映射管理")
    st.caption(
        "这里的操作不常做，但做完会立刻影响仿真页的可选项：序列库条目会出现在构建下拉框里，"
        "映射复核结果会出现在候选选择器的「实验时对应基因」列。"
    )

    st.markdown("**① 序列库**：信号肽 / 引导肽 / 成熟蛋白 / 组合模板的新建、修改、删除。")
    render_target_library_manager()

    st.divider()
    st.markdown("**② 基因 ↔ 复合体映射**：确认“过表达某复合体”对应实验室该动哪几个基因。")
    render_gene_complex_mapping_review()


__all__ = ["render_library_management"]
