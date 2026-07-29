"""仿真构建的分步向导：一次只展开一步、有上一步/下一步、运行按钮只在最后一步出现。

用户 2026-07-28 反馈标签页"太不明显"——三个标签地位平等，看不出先后，也没有"做完这步该干嘛"的指引。
这里用真 Streamlit 运行时验证向导行为，而不是只读源码。
"""

from __future__ import annotations

from pathlib import Path

from streamlit.testing.v1 import AppTest


REPO_ROOT = Path(__file__).resolve().parents[1]
APP_PATH = REPO_ROOT / "app" / "ui" / "streamlit_app.py"
STEP_KEY = "pichia_builder_step"


def _simulation_app(step: int | None = None) -> AppTest:
    app_test = AppTest.from_file(str(APP_PATH), default_timeout=120)
    app_test.session_state["app_page_nav_pending"] = "仿真验证"
    if step is not None:
        app_test.session_state[STEP_KEY] = step
    return app_test.run()


def _button_labels(app_test: AppTest) -> list[str]:
    return [button.label for button in app_test.button]


def test_first_step_offers_next_but_no_back_and_no_run() -> None:
    app_test = _simulation_app(1)

    labels = _button_labels(app_test)
    assert any("下一步" in label for label in labels), labels
    assert not any("上一步" in label for label in labels), labels
    assert not any("运行 Python 分泌仿真" in label for label in labels), "第 1 步不该能直接运行"


def test_middle_step_offers_both_directions_and_still_no_run() -> None:
    app_test = _simulation_app(2)

    labels = _button_labels(app_test)
    assert any("上一步" in label for label in labels), labels
    assert any("下一步" in label for label in labels), labels
    assert not any("运行 Python 分泌仿真" in label for label in labels), "第 2 步不该能直接运行"


def test_last_step_offers_run_and_back_but_no_next() -> None:
    app_test = _simulation_app(3)

    labels = _button_labels(app_test)
    assert any("运行 Python 分泌仿真" in label for label in labels), labels
    assert any("上一步" in label for label in labels), labels
    assert not any("下一步" in label for label in labels), labels


def test_only_the_current_step_is_expanded() -> None:
    """向导的核心：一次只展开当前这步，其余收起。"""
    markers = "①②③"
    for step in (1, 2, 3):
        app_test = _simulation_app(step)
        # 只认最外层的三个步骤容器；步骤内部的折叠区不得再用 ①②③ 编号（会撞号）。
        step_expanders = [item for item in app_test.expander if item.proto.label[:1] in markers]
        expanded = [item.proto.label for item in step_expanders if item.proto.expanded]

        assert len(step_expanders) == 3, [item.proto.label for item in step_expanders]
        assert len(expanded) == 1, f"第 {step} 步应只展开一个，实际 {expanded}"
        assert expanded[0].startswith(markers[step - 1]), (step, expanded)


def test_clicking_next_advances_the_step() -> None:
    app_test = _simulation_app(1)

    next_button = next(button for button in app_test.button if "下一步" in button.label)
    next_button.click().run()

    assert app_test.session_state[STEP_KEY] == 2


def test_prefilled_values_survive_navigating_back_and_forth() -> None:
    """从筛查页跳来的预填不能因为切步骤而消失。"""
    app_test = _simulation_app(1)
    app_test.session_state["pichia_draft_ko_genes"] = "PAS_chr2-2_0107"
    app_test.run()

    next(button for button in app_test.button if "下一步" in button.label).click().run()
    next(button for button in app_test.button if "上一步" in button.label).click().run()

    assert app_test.session_state["pichia_draft_ko_genes"] == "PAS_chr2-2_0107"


def test_user_typed_input_survives_navigating_steps() -> None:
    """**真实用户输入**必须保住——这是比预填更常见的路径，也是曾经真丢过的那条。

    实测：用 st.rerun() 切步会丢掉本轮刚敲进控件的值（填完基因点"下一步"，回头框是空的）；
    改用 on_click 回调后才正常。上一版测试只覆盖了"程序预填"，所以没抓到这个 bug。
    """
    app_test = _simulation_app(1)
    typed = next(area for area in app_test.text_area if area.label == "敲除基因（KO gene）")
    typed.set_value("PAS_chr2-2_0107").run()

    next(button for button in app_test.button if "下一步" in button.label).click().run()
    assert app_test.session_state["pichia_draft_ko_genes"] == "PAS_chr2-2_0107", "点下一步就丢了"

    next(button for button in app_test.button if "上一步" in button.label).click().run()
    back_on_step_one = next(area for area in app_test.text_area if area.label == "敲除基因（KO gene）")
    assert back_on_step_one.value == "PAS_chr2-2_0107", "回到第 1 步时输入框应还留着原内容"


def test_step_navigation_uses_callbacks_not_rerun() -> None:
    """锁住修法本身：切步一旦回退成"按钮分支里 st.rerun()"，用户输入又会开始丢。"""
    import inspect

    from app.ui.views import simulation

    source = inspect.getsource(simulation._render_step_nav)

    assert "on_click=_go_to_step" in source
    assert "st.rerun()" not in source


def test_step_indicator_tells_you_where_you_are() -> None:
    app_test = _simulation_app(2)

    captions = " ".join(item.value for item in app_test.caption)
    assert "第 2 / 3 步" in captions, captions


def test_navigation_sits_inside_the_current_step_not_at_page_bottom() -> None:
    """按钮必须紧跟刚做完的内容。放在三个折叠条之后（页面底部）时，做完第 1 步还得
    往下翻过另外两个面板才找到"下一步"，向导的引导感就没了。"""
    for step in (1, 2, 3):
        app_test = _simulation_app(step)
        nav_keys = [
            button.key for button in app_test.button
            if button.key and button.key.startswith(("pichia_builder_prev_", "pichia_builder_next_"))
        ]
        # 只有当前步渲染导航；其余步是收起状态，不该也不必带导航按钮
        assert all(key.endswith(str(step)) for key in nav_keys), (step, nav_keys)


def test_next_step_name_is_announced_so_you_know_what_is_coming() -> None:
    app_test = _simulation_app(1)

    captions = " ".join(item.value for item in app_test.caption)
    assert "下一步：改造候选" in captions, captions
