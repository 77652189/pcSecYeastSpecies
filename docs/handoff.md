# pcSecPichia Handoff

状态：active
最后更新：2026-07-15

## 当前执行位置

```yaml
current_phase: phase_2_gene_level_oe
current_round: round_6a_external_capacity_candidates
round_status: in_progress
current_checkpoint: a0c_ecpichia_provenance_closure
checkpoint_status: architecture_decision_required
```

## 当前结论

- A0c 已通过正式 source/cache/parser/evaluation/audit/CLI 路径完成 ecPichia G6PDH2 provenance closure，并支持无网络 replay。
- `Supplementary 8.yml` 的 `57689 / (8000 × 3600)` 可严格复算 GECKO coefficient；这只证明模型如何使用该 kcat，不证明 abundance 或当前模型 baseline capacity 成立。
- `Supplementary 11_V2.docx` 与 YAML 在 gene、enzyme、MW 和 concentration 上冲突；表格单位为 `g/L` 且 G6PDH2 concentration 为 `NaN`，YAML 则未声明单位。
- `kcat=8000 s^-1` 可追至 Thermotoga maritima、thio-NADP+、80°C 的 BRENDA 记录，不能直接适用于 Komagataella `PAS_chr2-1_0308 / C4R099`。
- hLF/OPN 已建立正式 `glucose_mu_0.1` current-model crosswalk，精确到 `G6PDH2_no_1_fwd` 和 `G6PDH2_no_1_fwd_complex_formation`；但 catalytic flux 到 formation/dilution `model_flux` 的换算仍无直接证据。
- 正式产物结论为 `architecture_decision_required`：`candidate_count=0`、`promotion_ready_count=0`、`nominal_capacity=null`、`promotion_preview_available=false`。
- 正式容量 registry 未修改，formal acceptance 仍因 `reviewed_baseline_capacity` 缺口为 `passed=false`；未进入 Round 6B 或 Phase 3。

## 下一项决策

先决定是否拆分 Phase 2 验收等级：

1. 保持现有硬门禁，只接受经审核的绝对 baseline capacity；Round 6A 继续阻塞，等待可闭合的同条件绝对 abundance/direct capacity 与 formation conversion 证据。
2. 新增独立的“相对、未校准 gene-level OE 场景”产品等级，同时保留绝对容量校准为单独未通过门禁；该等级不得解释为绝对产量或正式 capacity。

在决策完成前，不再增加只重复相对强度、名称或缺单位模型字段的 source adapter，不进入 Round 6B。

## 必读材料

1. [ADR-001：证据闭合与停止决策](adr/001-external-capacity-candidate-promotion.md#证据闭合与停止决策)
2. [执行计划：A0c 与状态迁移](pichia_next_plan.md#a0c-ecpichia-provenance-closure)
3. [数据与结果治理策略](data_and_results_policy.md#oe-capacity-参数与映射)
4. ignored 运行产物：`local_runs/oe_capacity/round6a/a0c_ecpichia_provenance/formal_run/g6pdh2_ecpichia_provenance_gap.json`

## 验证与停止线

- 重新运行 A0c focused tests、全部 `test_oe_capacity*.py`、service/UI contract、文档边界测试和 CLI offline replay。
- 运行 `compileall`、`git diff --check`、`local_runs/` ignore、密钥、依赖和保护目录检查。
- 不得用通用上界、baseline optimal flux、固定 `1.0`、PRIDE iBAQ、未确认单位的 ecPichia 值或 fixture 伪造容量。
- 不得进入 Round 6B、Phase 3，也不得生成 Phase 3 Round 0 提示词。
