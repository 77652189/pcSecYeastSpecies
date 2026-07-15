# pcSecPichia Handoff

状态：active
最后更新：2026-07-15

## 当前执行位置

```yaml
current_phase: phase_2_gene_level_oe
current_round: round_6a_external_capacity_candidates
round_status: in_progress
current_checkpoint: a0b_quantitative_source
checkpoint_status: in_progress_absolute_capacity_gap
```

## 当前事实

- A0a 结构收束已完成；模块和调用边界以 [ADR-001](adr/001-external-capacity-candidate-promotion.md#模块和调用边界) 为准。
- `external_candidates.py` 仅保留兼容重导出；source、schema、IO、evaluation、promotion 和 audit orchestration 已分离。
- CLI 与 `app/services` 均只通过公开 core API 调用候选 audit/review/promotion；既有公共导出、cache 格式、hash 和 UI/service 契约保持兼容。
- candidate review 使用同一 manifest/data snapshot；正式 promotion 使用审核 hash 和跨进程资产锁，不能静默覆盖并发更新。
- UniProt 只确认 `PAS_chr2-1_0308 / G6PDH2` identity，不提供 baseline capacity。
- 正式容量 registry 仍为空，正式 acceptance 仍为 `passed=false / reviewed_baseline_capacity`。
- A0b 已接入 PRIDE `PXD055501` MaxQuant 正式 parser/cache：项目版本 `2025-01-30`、许可 `CC0-1.0`，原始 proteinGroups SHA-256 为 `15b814790186146a1137353eceb83332493ccb826a47db0d7f781e1fe9084a26`。
- `F2QTE5 / ZWF1` 与当前模型 `C4R099 / PAS_chr2-1_0308` 序列 `504/504` 一致；T0 iBAQ 原始值为 `12868000 / 10476000 / 21552000`，在线获取和离线回放均可审计。
- iBAQ 仅是相对强度；来源条件为 glucose chemostat、`mu=0.075 h^-1`，与正式 `mu=0.1` 不匹配，并缺 absolute abundance、biomass normalization 和配对 kcat。因此候选为 `review_required`，`promotion_ready=false`，不能形成 `model_flux` capacity。
- A0b 已增加 ecPichia `Supplementary 8.yml` 正式文件导入和离线回放 adapter；YAML SHA-256 为 `317ab62f77c95feb2758f9ad7ed5efe18ff8430c747fbb880c03bb4d6b943d34`，上游 ZIP SHA-256 为 `bea45233dc4feb81295315c4e73ca2ca4c886f648822dda27347be8892a3620c`。
- 该 YAML 可重复提取 G6PDH2 的 gene `PAS_chr2-1_0308`、enzyme `C4R099`、MW `57689 g/mol`、kcat `8000 s^-1`、reaction coefficient `-0.00200309027777778`、reported concentration `0.752073171936811` 和 protein pool `-219.25`。
- ecPichia 证据仍有 supplement 表头单位与 GECKO 语义待协调、LFQ 到绝对丰度 provenance 缺失、kcat 仅标记 `brenda`、培养条件不完整、许可不可复用确认和 formation-flux 换算缺口，因此只进入 source assessment，`promotion_ready=false`。

## 下一 Checkpoint：继续 A0b source acquisition

继续寻找可把 G6PDH2 定量证据闭合为 absolute、biomass-normalized baseline capacity 的来源，或补齐 condition-matched abundance + kcat 换算链。来源优先级、许可、provenance、条件匹配和 promotion 门禁见：

1. [ADR-001：决策与来源优先级](adr/001-external-capacity-candidate-promotion.md#决策)
2. [数据与结果政策](data_and_results_policy.md)
3. [执行计划：Round 6A](pichia_next_plan.md#round-6a外部-baseline-capacity-候选与审核提升)

现有 `PXD055501` candidate 和 ecPichia raw values 只可作为待复核证据，不得提升为正式容量。下一步优先获取直接 Komagataella G6PDH2 absolute abundance、condition-matched direct kcat，或补齐 ecPichia LFQ→mg/gDCW provenance 与当前 formation flux 换算；不得把 UniProt identity、BLAST/RBH、人工 smoke、通用上界、optimal flux 或 fixture 当作容量。Round 6A 保持 `in_progress`。

## 验证方式

- external-candidate focused tests 和全部 `test_oe_capacity*.py`
- service/UI contract 与 `tests/test_docs_active_boundary.py`
- `python -m compileall -q python_pichia/src`
- CLI/import 边界、`git diff --check`、`local_runs/` ignore 和保护目录检查

## 停止线

当前已完成 A0b 的 PRIDE relative candidate 与 ecPichia formal source-assessment 路径，但 absolute capacity 缺口未闭合；未执行正式 promotion，未开始 Round 6B 或 Phase 3，也不生成 Phase 3 提示词。
