"""界面显示用的**统一命名入口**：反应 id / 基因 id → 研究员看得懂的名字。

**为什么要收拢**：命名逻辑此前散在好几处——基因名在 `gene_name_annotation`、复合体俗名在候选面板、
反应 id 在图表里各自 `_short_reaction` 截断、模型自带的反应名压根没被读进来。同一个反应在不同页面
显示成不同样子，改一处不会全局生效。这里是唯一的解析入口，各页面只调它。

反应名按**可信度**依次尝试，全部落空才退回 id：

1. **模型自带名称**（`Model/pcSecPichia.mat` 的 `rxnNames`，29026 条）——最权威，代谢反应几乎全覆盖；
2. **策展俗名**（`SECRETION_GENE_CATALOG`）——分泌机器复合体在 .mat 里**没有**名称，它们的可读名
   （PDI1 / KAR2 / OCH1）只存在于策展库；
3. **借基因名**——形如 `PAS_chr1-1_0187_dilution_misfolding` 的反应内嵌基因位点，可借该基因的蛋白注释。

诚实边界：模型作者把一部分反应直接标成 `unclear reaction`，这类会显式标注出来——
研究员不该把模型自己都没标清的反应当靶点。
"""

from __future__ import annotations

import re
from functools import lru_cache

from app import ensure_python_pichia_on_path
from app.ui.common import PATHS


UNCLEAR_REACTION_NAME = "unclear reaction"
UNCLEAR_REACTION_LABEL = "模型未标注用途"
_GENE_TOKEN = re.compile(r"PAS_chr[\w\-]*?\d(?:_\d+)?")
# 图表 y 轴 / 按钮标签的长度上限：模型里有长达上百字符的反应名，不截断会撑爆图表。
MAX_REACTION_LABEL_CHARS = 44


@lru_cache(maxsize=1)
def _model_reaction_names() -> dict[str, str]:
    """模型自带的反应名（一次性加载并缓存）。加载失败时退化为空字典，不影响页面。"""
    ensure_python_pichia_on_path()
    try:
        from pcsec_pichia.probe import load_pcsec_pichia_model

        model = load_pcsec_pichia_model(PATHS.repo_root)
    except Exception:  # noqa: BLE001 - 命名只是显示增强，读不到不该拖垮任何页面
        return {}
    return {
        reaction_id: str(name).strip()
        for reaction_id, name in zip(model.rxns, model.rxn_names)
        if str(name or "").strip() and str(name).strip() != reaction_id
    }


@lru_cache(maxsize=1)
def _curated_reaction_names() -> dict[str, str]:
    """策展库里复合体反应的俗名（PDI1 / KAR2 / OCH1……）。分泌机器在 .mat 里没有名称。"""
    ensure_python_pichia_on_path()
    try:
        from pcsec_pichia.services.gene_catalog import SECRETION_GENE_CATALOG
    except Exception:  # noqa: BLE001
        return {}
    names: dict[str, str] = {}
    for entry in SECRETION_GENE_CATALOG:
        common_name = str(getattr(entry, "common_name", "") or "").strip()
        if not common_name:
            continue
        for attribute in ("oe_reaction_id", "ko_reaction_id"):
            reaction_id = str(getattr(entry, attribute, "") or "").strip()
            if not reaction_id:
                continue
            existing = names.get(reaction_id)
            # 多条策展条目指向同一复合体（PDI1/ERO1/ERV2）——合并成一个可读标签，不覆盖。
            if existing and common_name not in existing.split("/"):
                names[reaction_id] = f"{existing}/{common_name}"
            elif not existing:
                names[reaction_id] = common_name
    return names


def _gene_borrowed_name(reaction_id: str) -> str:
    """反应 id 里内嵌了基因位点时，借该基因的蛋白注释。"""
    match = _GENE_TOKEN.search(reaction_id)
    if not match:
        return ""
    try:
        from app.services.gene_name_annotation import load_standard_name_lookup

        record = load_standard_name_lookup().get(match.group(0))
    except Exception:  # noqa: BLE001
        return ""
    if record is None:
        return ""
    return str(record.standard_symbol or record.display_name or "").strip()


def reaction_display_name(reaction_id: str) -> str:
    """反应的可读名称；查不到返回空串。三个来源按可信度依次尝试。"""
    resolved = str(reaction_id or "").strip()
    if not resolved:
        return ""
    name = _model_reaction_names().get(resolved)
    if name:
        return UNCLEAR_REACTION_LABEL if name.lower() == UNCLEAR_REACTION_NAME else name
    curated = _curated_reaction_names().get(resolved)
    if curated:
        return curated
    return _gene_borrowed_name(resolved)


def reaction_display_label(reaction_id: str, *, max_chars: int = MAX_REACTION_LABEL_CHARS) -> str:
    """给图表/表格用的标签：有名字就「名称（id）」，没有就只给 id。

    **id 始终在场**——它才是模型实际算的对象，研究员要靠它跟仿真输入对上。
    """
    resolved = str(reaction_id or "").strip()
    name = reaction_display_name(resolved)
    if not name:
        return _shorten(resolved, max_chars)
    return f"{_shorten(name, max_chars)}（{resolved}）"


def _shorten(text: str, max_chars: int) -> str:
    if max_chars <= 0 or len(text) <= max_chars:
        return text
    return text[: max_chars - 1].rstrip(" ,;:-") + "…"


__all__ = [
    "MAX_REACTION_LABEL_CHARS",
    "UNCLEAR_REACTION_LABEL",
    "reaction_display_label",
    "reaction_display_name",
]
