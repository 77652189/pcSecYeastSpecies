"""整页冒烟测试：用**真 Streamlit 运行时**把每个页面渲染一遍，断言不崩。

为什么必须是真运行时（2026-07-28 立此规矩）：此前所有 UI 验证都是把 `st` 桩掉后调用渲染函数——
那只能证明数据对，证明不了页面能跑。结果是一连串问题全靠用户看截图发现：
`has no attribute`（进程内模块过旧）、重复控件 key、在控件实例化后写 session_state……
这些都只在真运行时才会暴露。桩测试**不能**替代本文件。

本测试不做 LP 求解，只渲染，所以默认运行（不进慢测网关）。
"""

from __future__ import annotations

from pathlib import Path

import pytest

from streamlit.testing.v1 import AppTest


REPO_ROOT = Path(__file__).resolve().parents[1]
APP_PATH = REPO_ROOT / "app" / "ui" / "streamlit_app.py"
PENDING_NAV_KEY = "app_page_nav_pending"

# 侧边栏的全部页面。新增页面必须加进来——否则它就是下一个"用户先发现"的地方。
PAGES = (
    "项目总览",
    "结果浏览",
    "全基因组KO/OE筛查",
    "仿真验证",
    "实验反馈闭环",
    "基因命名与同源规则审计",
    "Shadow LP一致性验证",
    "运行日志",
)


def _run_page(page: str, *, timeout: int = 120) -> AppTest:
    app_test = AppTest.from_file(str(APP_PATH), default_timeout=timeout)
    app_test.session_state[PENDING_NAV_KEY] = page
    return app_test.run()


@pytest.mark.parametrize("page", PAGES)
def test_every_page_renders_without_exception(page: str) -> None:
    app_test = _run_page(page)

    assert not app_test.exception, (
        f"「{page}」渲染抛异常：\n"
        + "\n".join(str(item.value) for item in app_test.exception)
    )


@pytest.mark.parametrize("page", PAGES)
def test_every_page_renders_actual_content(page: str) -> None:
    """不崩还不够——空白页也不崩。每页都得真渲染出东西。"""
    app_test = _run_page(page)

    produced = (
        len(app_test.markdown) + len(app_test.dataframe) + len(app_test.metric)
        + len(app_test.caption) + len(app_test.subheader) + len(app_test.header)
    )
    assert produced > 0, f"「{page}」没有渲染出任何内容"


def test_simulation_page_exposes_the_primary_ko_oe_inputs() -> None:
    """仿真验证是主流程：两个主输入框必须在场（此前它们被辅助面板挤到屏幕外过）。"""
    app_test = _run_page("仿真验证")

    labels = [widget.label for widget in app_test.text_area]
    assert any("敲除基因" in label for label in labels), labels
    assert any("过表达基因" in label for label in labels), labels


def test_genome_screen_exposes_scope_choice_and_no_stale_scope_label() -> None:
    app_test = _run_page("全基因组KO/OE筛查")

    radio_options = [option for widget in app_test.radio for option in widget.options]
    assert any("全基因组" in option for option in radio_options), radio_options
    # 曾把 61 个反应写成"约30个"，少报一半
    assert not any("约30个反应" in option for option in radio_options), radio_options


def test_no_duplicate_widget_keys_across_the_simulation_page() -> None:
    """重复 key 会让 Streamlit 直接报 DuplicateWidgetID——桩测试完全抓不到。
    仿真验证页控件最多（选择器 + 序列库 + 策展复核 + 输入框），最容易撞。"""
    app_test = _run_page("仿真验证")

    assert not app_test.exception, [str(item.value) for item in app_test.exception]
