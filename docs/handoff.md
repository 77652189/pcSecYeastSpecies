# pcSecPichia Handoff

状态：active
最后更新：2026-07-23

## 当前目标

当前 slice：**改造后 per-strain 瓶颈 → 下一步 OE 候选**——复用 per-solve 瓶颈归因（R1）+ 有界剂量响应（R2）接进仿真验证，让已改造的菌株也能"再找瓶颈、排下一候选"。详见 [执行计划](EXECUTION_PLAN.md) 阶段② + [ADR-004](adr/004-relative-signal-deepening-under-permanent-data-gap.md)。方向5（**碳源条件标定 + 跨条件稳健性**）主体已完成（阶段①：A 全 + B1/B2/B4/B5-μ/私有护栏，短名单跨碳源稳健 0 翻转），剩 titer 锚点 / 方向1 摄入等数据门控尾巴。

```yaml
current_slice: per_strain_next_oe_candidate_readout
slice_status: in_progress
current_program: mvp_directions_plus_relative_signal_deepening
absolute_capacity_status: unavailable_waiting_for_qualified_evidence
carbon_source_calibration_status: glucose_corrected_others_internally_calibrated
wetlab_fermentation_data_status: in_hand_local_private_store_out_of_repo
rnaseq_expression_constraint_status: contract_ready_adr005_waiting_for_data
direction_4_combination_status: strategically_deferred
```

一句话背景：模型是"折叠/分泌容量这一片"的相对候选生成器，不预测绝对产量；绝对容量数据永久缺失。五方向完整状态、分层架构见 [需求与架构文档](pichia_current_architecture_and_requirements.md)，历史逐次记录见 git。

## 下一步

当前 slice（#1 迭代候选，清单见执行计划阶段②）：
- **C1 service**（纯逻辑 + 单测）：从 solve 结果 `oe_actionable_bottlenecks`（top-N binding 复合体）派生候选 → 有界 OE 剂量响应（复用 R2）→ 按效应 + 形状排序。
- **C2 编排 + C3 UI**：两趟 pipeline（solve 拿瓶颈 → 对瓶颈复合体 dose-response），opt-in；仿真验证页"下一步 OE 候选（针对当前改造菌株）"section（复合体级 + 分泌层 + 诚实 caveat）。
- **C4**：测试 + app 验证。

方向5 收尾（数据门控、暂缓）：B5 titer 锚点（待验证数据）+ 方向1 本地摄入（经护栏读私有区）、B4 湿实验一致性标注。B3 噪声门控暂不建（0 翻转、无表观敏感可甄别）。

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
