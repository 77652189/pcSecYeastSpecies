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

    offenders: list[str] = []
    for root in checked_roots:
        for path in root.rglob("*.py"):
            text = path.read_text(encoding="utf-8", errors="replace")
            if any(fragment in text for fragment in bad_fragments):
                offenders.append(str(path))

    assert offenders == []


def test_opn_page_separates_secretion_model_from_cds_design() -> None:
    text = _read_ui_text(
        "app/ui/pages/opn.py",
        "app/ui/pages/opn_cds.py",
    )

    assert "pcSecPichia 分泌模型筛选" in text
    assert "PichiaCLM 生成毕赤酵母 CDS" in text
    assert "PichiaCLM 不参与分泌模型评分" in text


def test_opn_page_explains_method_comparison_without_default_signalp() -> None:
    text = _read_ui_text("app/ui/pages/opn_library.py")

    assert "从 UniProt 建库并比较筛选方法" in text
    assert "USPNet-fast" in text
    assert "自研规则" in text


def test_streamlit_entrypoint_keeps_all_demo_pages() -> None:
    text = Path("app/ui/streamlit_app.py").read_text(encoding="utf-8")

    assert "项目总览" in text
    assert "结果浏览" in text
    assert "OPN 信号肽" in text
    assert "仿真验证" in text
    assert "运行日志" in text
    assert "render_overview" in text
    assert "render_results_browser" in text
    assert "render_opn_signal_peptides" in text
    assert "render_simulation" in text
    assert "render_logs" in text


def _read_ui_text(*paths: str) -> str:
    return "\n".join(Path(path).read_text(encoding="utf-8") for path in paths)
