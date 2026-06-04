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
