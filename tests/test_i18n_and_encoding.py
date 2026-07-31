from __future__ import annotations

from pathlib import Path

from app.core.i18n import category_label, species_label, status_label


def test_chinese_label_mappings() -> None:
    assert species_label("SCE") == "酿酒酵母（S. cerevisiae）"
    assert category_label("CSource") == "碳源分析"
    assert status_label("ok") == "正常"
    assert status_label(True) == "求解成功"


def test_user_facing_python_files_have_no_mojibake() -> None:
    checked_roots = [Path("app/ui"), Path("app/services"), Path("app/core")]
    bad_fragments = ["�", "閰", "绋", "鐢", "鍙", "杩", "鏂", "浠跨湡", "鍛戒护"]

    # 防空转：目录被改名/移走时 rglob 静默返回空，offenders == [] 会 vacuously 通过。
    for root in checked_roots:
        assert root.is_dir(), f"扫描目标不存在，测试会静默通过：{root}"

    offenders: list[str] = []
    for root in checked_roots:
        for path in root.rglob("*.py"):
            text = path.read_text(encoding="utf-8", errors="replace")
            if any(fragment in text for fragment in bad_fragments):
                offenders.append(str(path))

    assert offenders == []


def test_streamlit_entrypoint_keeps_all_demo_pages() -> None:
    text = Path("app/ui/streamlit_app.py").read_text(encoding="utf-8")
    common_text = Path("app/ui/common.py").read_text(encoding="utf-8")

    assert "项目总览" in text
    assert "结果浏览" in text
    assert "仿真验证" in text
    assert "基因命名与同源规则审计" in common_text
    assert "HOMOLOGY_AUDIT_PAGE" in text
    assert "运行日志" in text
    assert "render_overview" in text
    assert "render_results_browser" in text
    assert "render_simulation" in text
    assert "render_homology_audit" in text
    assert "render_logs" in text


def test_streamlit_internal_modules_do_not_use_reserved_pages_directory() -> None:
    assert not Path("app/ui/pages").exists()


def test_app_css_does_not_fetch_external_font_assets() -> None:
    text = Path("app/ui/common.py").read_text(encoding="utf-8")

    assert "@import url(" not in text
    assert "fonts.googleapis.com" not in text
    assert "letter-spacing: -" not in text


def _read_ui_text(*paths: str) -> str:
    return "\n".join(Path(path).read_text(encoding="utf-8") for path in paths)
