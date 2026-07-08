# SCE/Pichia BLAST 同源可行性验证

日期：2026-07-08
分支：`codex/sce-pichia-homology-feasibility`

## 验证结论

可行。当前仓库已经具备离线 BLAST/RBH 同源证据链的主体基础设施，可以用本地 S. cerevisiae 与 Pichia 全蛋白组序列，对 `SECRETION_GENE_CATALOG` 中的酿酒酵母常用名生成结构化 crosswalk、name audit 和 rule-transfer audit。

这条路径比 UniProt/KEGG 命名匹配更可靠的原因在本次试跑中得到验证：它先用序列双向最佳命中建立候选，再把低 identity、缺 reciprocal hit、未进入当前 Pichia GEM gene_index 的情况降级到人工复核或不可迁移，而不是仅凭名称相似性直接接入模型解释。

## 本次试跑

命令：

```powershell
python scripts\build_pichia_homology_cache.py `
  --catalog-only `
  --blast-bin "local_runs\blast_homolog_feasibility\bin\ncbi-blast-2.17.0+\bin" `
  --output-dir "local_runs\sce_pichia_homology_feasibility_20260708"
```

输入与边界：

- SCE 蛋白序列：`Data/pcSecYeast/Protein_Sequence.mat`
- Pichia 蛋白序列：`Data/pcSecPichia/Protein_Sequence.mat`
- 查询来源：`SECRETION_GENE_CATALOG` 派生出的 SCE symbols/aliases
- BLAST/RBH：本地 NCBI BLAST+ 2.17.0，未联网
- 阈值：`min_identity=30.0`、`min_coverage=50.0`、`max_evalue=1e-10`
- 产物位置：`local_runs/sce_pichia_homology_feasibility_20260708/`，由 `.gitignore` 排除

关键产物：

- `sce_to_pichia_homology_cache.jsonl`
- `sce_to_pichia_homology_cache.tsv`
- `sce_to_pichia_name_audit.jsonl`
- `sce_to_pichia_rule_transfer_audit.jsonl`
- `homology_audit_summary.json`

## 结果摘要

`homology_audit_summary.json`：

| 指标 | 数值 |
| --- | ---: |
| BLAST 状态 | `completed` |
| homology rows | 106 |
| RBH rows | 75 |
| model-ready RBH high confidence | 14 |
| RBH but not in current Pichia GEM gene_index | 53 |
| low identity review required | 7 |
| no reciprocal hit | 4 |
| coverage review required | 1 |
| unresolved query symbol | 27 |

rule-transfer audit 同步给出：

| 状态 | 数值 |
| --- | ---: |
| `rule_transfer_ready` | 14 |
| `rule_transfer_supported_not_model_operable` | 53 |
| `rule_transfer_low_confidence` | 8 |
| `rule_transfer_not_supported` | 4 |
| `rule_transfer_unresolved` | 27 |

## 关键基因 spot check

| SCE symbol | SCE ORF | Pichia candidate | RBH | in model gene_index | identity | review_status | 解释 |
| --- | --- | --- | --- | --- | ---: | --- | --- |
| `KAR2` | `YJL034W` | `PAS_chr2-1_0140` | true | true | 75.532 | `model_ready_rbh_high_confidence` | 可作为高可信同源候选进入后续模型复核 |
| `VPS1` | `YKR001C` | `PAS_chr1-4_0644` | true | true | 71.023 | `model_ready_rbh_high_confidence` | 可作为高可信同源候选进入后续模型复核 |
| `CDC48` | `YDL126C` | `PAS_FragD_0026` | true | true | 83.751 | `model_ready_rbh_high_confidence` | 可作为高可信同源候选进入后续模型复核 |
| `PDI1` | `YCL043C` | `PAS_chr4_0844` | true | false | 45.233 | `rbh_not_in_model` | 序列支持同源，但当前 GEM 不可直接执行该 gene_id |
| `DOA10` | `YIL030C` | `PAS_chr3_0538` | true | false | 29.022 | `low_identity_review_required` | RBH 存在，但低于 identity 阈值且不在模型中，必须人工复核 |
| `HRD1` | `YOL013C` | `PAS_chr4_0156` | true | false | 24.411 | `low_identity_review_required` | RBH 存在，但低于 identity 阈值且不在模型中，必须人工复核 |

## 可行性判断

这次验证说明：

- 本地全蛋白组序列源存在且可解析，不需要依赖不稳定的命名注释作为第一证据。
- builder 已经支持 SCE->Pichia 正向 BLAST、Pichia->SCE 反向 BLAST、RBH 计算、阈值配置、Pichia GEM gene_index join。
- 结果不会把强同源直接升级为 phenotype evidence；它只产生 `review_status` 和 `rule_transfer_status`。
- 对 PDI1、DOA10、HRD1 这类“生物学上值得关注但模型不可直接操作或 identity 偏低”的候选，流程会保守降级，而不是自动放进 KO/OE 可执行集合。
- Streamlit/runtime 应继续只读 cache；BLAST 只能作为离线 builder 或人工验证步骤运行。

## 下一步建议

推荐下一步是保留当前 offline cache builder 边界，补一层“小范围 homology validation report”入口或命令，用固定的 10-20 个 SCE symbols 生成可复核 TSV/Markdown，作为后续把 SCE 文献日志接入候选库前的审查门。

备选路径：

- 将本次 `local_runs/sce_pichia_homology_feasibility_20260708/` 产物提升为一次审计样本，但不要作为稳定生产 cache。
- 追加离线 external-name reference cache，用 UniProt/SGD/NCBI 只做名称交叉复核，不覆盖 RBH 事实。
- 针对 `rbh_not_in_model` 的 53 个候选，单独做模型 gene/reaction/GPR 可操作性复核，避免把同源证据误当成仿真可执行证据。
