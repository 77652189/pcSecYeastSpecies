# pcSecPichia Handoff

状态：active
最后更新：2026-07-23

## 当前目标

当前 slice：**改造后候选系统（OE + KO · 分层复用）**——迭代1（下一步 OE 候选：R1 瓶颈 + R2 剂量响应接进仿真验证）已完成（本地 `d2b4cd0` 未 push）；迭代2 把它扩成 OE+KO 完整短名单 + **分层复用**（改造只重算受影响层、其余复用野生型缓存）+ **全基因组后台层**，整合进仿真验证页。详见 [执行计划](EXECUTION_PLAN.md) 阶段② 迭代2 + [ADR-004](adr/004-relative-signal-deepening-under-permanent-data-gap.md)。方向5（碳源条件标定 + 跨条件稳健性）主体已完成，剩 titer 锚点 / 方向1 摄入等数据门控尾巴。

```yaml
current_slice: modified_strain_ko_oe_layered_shortlist
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

阶段② 迭代1（C1–C4，下一步 OE 候选）✅ **已完成**（本地 `d2b4cd0` 未 push）：`strain_modifications` 叠 KO/OE + `solve_secretion_capacity`/`run_oe_dose_response_sweep` opt-in 参数（默认 None → glucose 逐字不变）；引擎 `next_oe_candidates` + 服务 `per_strain_oe_candidate_run` + 仿真验证页 UI；586/342 全绿、app 实跑通过。

迭代2（改造后候选系统 · 分层复用 + 全基因组后台，用户 2026-07-24 拍板；清单见执行计划阶段② 迭代2）：
- **D1 KO 候选引擎**（当前起点）：`run_knockout_screen` 加 opt-in `strain_modifications`（改造后基线跑 KO 扰动，与 sweep 同款 additive）→ 引擎产改造后 KO 候选。
- **D2/D3 L1**：改造后瓶颈（C2）对比野生型 R1 缓存 → 受影响层判定 + 给野生型短名单打标（复用/失效）；即时复用 + 按需重算受影响层。
- **D4 L3**：改造后全基因组 KO/OE 离线工具 + 缓存（菌株指纹 key，复用 B1 fingerprint 思路）。
- **D5 UI** 整合进仿真验证页；**D6** 测试/验证/文档。
- 核心前提（诚实）：复用是**近似**（同层/紧邻才需重算），须显式标注失效范围；**KO 无免费派生、必求解**（OE 有影子价格捷径）。

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
