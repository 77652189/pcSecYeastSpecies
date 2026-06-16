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
    text = Path("app/ui/streamlit_app.py").read_text(encoding="utf-8")

    assert "pcSecPichia 分泌模型筛选" in text
    assert "PichiaCLM 下游 CDS 设计" in text
    assert "PichiaCLM 不参与分泌模型评分" in text


def test_opn_page_explains_method_comparison_without_default_signalp() -> None:
    text = Path("app/ui/streamlit_app.py").read_text(encoding="utf-8")

    assert "从 UniProt 建库并比较筛选方法" in text
    assert "USPNet-fast" in text
    assert "自研规则" in text
