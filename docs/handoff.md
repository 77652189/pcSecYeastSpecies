# pcSecPichia Handoff

状态：active
最后更新：2026-07-14

## 当前执行位置

```yaml
current_phase: phase_2_gene_level_oe
current_round: round_6a_external_capacity_candidates
round_status: in_progress
current_checkpoint: a0b_quantitative_source
checkpoint_status: ready
```

## 当前事实

- A0a 结构收束已完成；模块和调用边界以 [ADR-001](adr/001-external-capacity-candidate-promotion.md#模块和调用边界) 为准。
- `external_candidates.py` 仅保留兼容重导出；source、schema、IO、evaluation、promotion 和 audit orchestration 已分离。
- CLI 与 `app/services` 均只通过公开 core API 调用候选 audit/review/promotion；既有公共导出、cache 格式、hash 和 UI/service 契约保持兼容。
- candidate review 使用同一 manifest/data snapshot；正式 promotion 使用审核 hash 和跨进程资产锁，不能静默覆盖并发更新。
- UniProt 只确认 `PAS_chr2-1_0308 / G6PDH2` identity，不提供 baseline capacity。
- 正式容量 registry 仍为空，正式 acceptance 仍为 `passed=false / reviewed_baseline_capacity`。
- 当前没有真实 Pichia 定量蛋白组、iPichia/ecPichia capacity 或 abundance+kcat 候选。

## 下一 Checkpoint：A0b

目标是接入至少一个可审计的真实定量外部来源，并生成 G6PDH2 数值候选。来源优先级、许可、provenance、条件匹配和 promotion 门禁见：

1. [ADR-001：决策与来源优先级](adr/001-external-capacity-candidate-promotion.md#决策)
2. [数据与结果政策](data_and_results_policy.md)
3. [执行计划：Round 6A](pichia_next_plan.md#round-6a外部-baseline-capacity-候选与审核提升)

A0b 不得把 UniProt identity、人工导入 smoke、source inventory 或 synthetic fixture 当作真实定量来源。没有满足许可、hash、condition、单位换算和 current-model binding 的候选时，Round 6A 保持 `in_progress`。

## 验证方式

- external-candidate focused tests 和全部 `test_oe_capacity*.py`
- service/UI contract 与 `tests/test_docs_active_boundary.py`
- `python -m compileall -q python_pichia/src`
- CLI/import 边界、`git diff --check`、`local_runs/` ignore 和保护目录检查

## 停止线

当前只准备 A0b；尚未接入真实定量来源、未执行正式 promotion、未开始 Round 6B 或 Phase 3，也不生成 Phase 3 提示词。
