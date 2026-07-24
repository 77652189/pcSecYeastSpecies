# pcSecPichia Handoff

状态：active
最后更新：2026-07-23

## 当前目标

当前 slice：**碳源条件标定 + 跨条件稳健性**（方向5 有界升级）。详见 [项目级执行计划](EXECUTION_PLAN.md) 阶段① 与 [ADR-006](adr/006-carbon-source-condition-calibration.md)。

```yaml
current_slice: direction_5_carbon_source_calibration_and_condition_robustness
slice_status: in_progress
current_program: mvp_directions_plus_relative_signal_deepening
absolute_capacity_status: unavailable_waiting_for_qualified_evidence
carbon_source_calibration_status: glucose_corrected_others_draft_pending_internal_calibration
wetlab_fermentation_data_status: in_hand_local_private_store_out_of_repo
rnaseq_expression_constraint_status: contract_ready_adr005_waiting_for_data
direction_4_combination_status: strategically_deferred
```

一句话背景：模型是"折叠/分泌容量这一片"的相对候选生成器，不预测绝对产量；绝对容量数据永久缺失。五方向完整状态、分层架构见 [需求与架构文档](pichia_current_architecture_and_requirements.md)，历史逐次记录见 git。

## 下一步

阶段① 任务（清单见执行计划）：
- **A 碳源标定 — 已完成**：蛋白含量条件化（甲醇 0.40）、蛋白成本/生长约束认 formulation 选定的生长反应、核实 `*_meoh` 生物量组成、五条件内部验证、三档 formulation 状态（非葡萄糖=`internally_calibrated`）+ UI 标注；glucose 逐字不变（护栏 10 passed）。
- **B 跨条件稳健性 + 数据接线 — 进行中**：已完成 B1 求解缓存 / B2 条件矩阵 / B4 跨条件稳健性面板 / B5 μ 验证（默认 μ=0.10 = hLF 甘油生长相实测）。**B3 噪声门控暂不建**（真实工艺矩阵 0 翻转、hLF 短名单跨碳源全稳健、无表观敏感可甄别，未来出现敏感再建）。
- **剩余（下一步）**：B5 在手发酵数据本地接线（gitignored 私有路径 + 提交护栏 + titer 锚点 + UPR×折叠）、B4 湿实验一致性标注（依赖前者）。

## 范围边界（硬约束）

- 无界完整跨条件排名产品仍不做；只产相对信号，绝对容量恒 unavailable。
- **glucose 的 corrected_reference 结果不得改动**（回归锁定）。
- **保密湿实验数据只存仓库外本地私有区**（`CursorProject/pcSec_wetlab_private/`）、提交产物只含机制层抽象、不上传云端/GitHub。
- 方向4 组合搜索、目标蛋白降解通路建模、换默认 solver：明确不做。

## 必读材料

1. [项目级执行计划](EXECUTION_PLAN.md)：当前阶段、任务、范围边界。
2. [需求与架构文档](pichia_current_architecture_and_requirements.md)：五方向状态、分层架构、碳源标定状态、数据与产物治理。
3. [ADR-006](adr/006-carbon-source-condition-calibration.md)：碳源标定 + 数据契约。
4. [ADR-002](adr/002-relative-oe-and-absolute-capacity-layers.md) / [ADR-004](adr/004-relative-signal-deepening-under-permanent-data-gap.md) / [ADR-005](adr/005-rnaseq-expression-constrained-enzyme-capacity.md)。

## 验证方式

- 文档锚点：`tests/test_docs_active_boundary.py`、`test_data_results_boundaries.py`、`test_slow_test_gates.py`。
- 阶段① 落地后：碳源标定内部验证 + **glucose 回归结果逐字不变** + 短名单/契约测试。
- 最近基线（改文档前）：根 `tests/` 334 passed；`python_pichia/tests/` 全量隔离 555 passed / 21 skipped / 0 failed（见 git 历史）。
