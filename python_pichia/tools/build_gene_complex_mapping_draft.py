"""从同源比对结果自动起草「基因 ↔ 复合体」映射，供人工复核（ADR-007 / E3）。

把策展工作从"从零查 78 个复合体"降成"审一批勾选题"。草稿一律标最保守的
`pending_review` + `auxiliary` + `unknown`，**不会自己生效**——必须有人复核后才可能进入可执行路径。

输出默认落 `local_runs/`（运行产物目录）。复核完成后由**人**显式把成品放到
`Data/pcSecPichia/gene_complex_mapping.json`——策展映射是长期科学资产，
按数据治理不能由脚本自动写入受保护目录。

用法：
    python python_pichia/tools/build_gene_complex_mapping_draft.py
    python python_pichia/tools/build_gene_complex_mapping_draft.py --out 某个路径.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT = REPO_ROOT / "local_runs" / "gene_complex_mapping_draft" / "gene_complex_mapping.draft.json"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT, help="草稿输出路径")
    args = parser.parse_args(argv)

    sys.path.insert(0, str(REPO_ROOT))
    sys.path.insert(0, str(REPO_ROOT / "python_pichia" / "src"))

    from app.services.pichia_gene_catalog_service import load_hlf_opn_candidate_genes
    from pcsec_pichia.services.gene_complex_mapping import (
        build_draft_mappings_from_candidates,
        serialize_gene_complex_mappings,
    )

    candidates = load_hlf_opn_candidate_genes(target_context=None, include_shared=True)
    drafts = build_draft_mappings_from_candidates(candidates)
    payload = serialize_gene_complex_mappings(drafts)

    output_path = Path(args.out)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")

    complexes = {row.complex_reaction_id for row in drafts}
    genes = {row.pichia_gene_id for row in drafts}
    print(f"起草 {len(drafts)} 条映射：{len(genes)} 个基因 × {len(complexes)} 个复合体反应")
    print(f"输出：{output_path}")
    print("全部标为 pending_review（不会自己生效）。复核后由人放到 Data/pcSecPichia/gene_complex_mapping.json。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
