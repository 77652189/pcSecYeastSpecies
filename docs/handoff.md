# pcSecPichia Handoff

状态：active

## 当前目标

当前 slice：**可用性 + 可达性（阶段④）**——让已经建好的生成能力真正能被研究员用上。见 [执行计划](EXECUTION_PLAN.md) 阶段④ 与 [ADR-007](adr/007-secretory-machinery-gene-complex-reachability.md)。

> **2026-07-28 重要更正（别再照抄旧结论）**：此前 handoff 与执行计划都写"剩余有意义工作全部数据门控"——**该判断已不成立**。易用性走查查实：策展 69 条分泌机器候选里只有 13 条能按基因跑，其余 **56 条全部带着模型里真实存在的复合体反应**（含 `sec_Pdi1p_complex_formation`＝hLF 头号 OE 杠杆 +8.15%），而界面把它们显示成"仅复核/模型外"的死胡同、且全基因组基因筛查**覆盖这些反应 0 个**。即：**真正的杠杆区此前在界面上基本不可达**，这是一条实质性的、**不需要任何湿实验数据**的工作线。

阶段④ 的已完成项（E0–E3）与逐条 commit 见 [执行计划](EXECUTION_PLAN.md) 阶段④——本文不复述。

下一步：**策展数据到位**——把映射写进 `Data/pcSecPichia/gene_complex_mapping.json` 即生效、无需改代码（该文件按设计**尚未创建**，等策展拍板）。**待拍板：策展范围与责任人**（哪些复合体先做、谁判断映射可信；生物学判断，非软件问题）。阶段④ 软件侧至此收口。

前一 slice（改造后候选系统 · 分层复用）已收口。结论与逐 D 细节见 [执行计划](EXECUTION_PLAN.md)「已完成能力」。

```yaml
current_slice: usability_and_secretory_machinery_reachability
slice_status: in_progress
previous_slice: modified_strain_ko_oe_layered_shortlist
previous_slice_status: done
current_program: mvp_directions_plus_relative_signal_deepening
project_stage: generation_side_complete_usability_reachability_in_progress
reachability_gap_status: secretory_complexes_runnable_but_were_unreachable_in_ui_adr007
absolute_capacity_status: unavailable_waiting_for_qualified_evidence
carbon_source_calibration_status: glucose_corrected_others_internally_calibrated
wetlab_fermentation_data_status: in_hand_is_condition_validation_only_no_modification_titer_screen
experiment_feedback_validation_status: machinery_ready_data_gated_no_modification_titer_dataset
rnaseq_expression_constraint_status: contract_ready_adr005_waiting_for_data
direction_4_combination_status: strategically_deferred
```

一句话背景：模型是“折叠 / 分泌容量这一片”的相对候选生成器，不预测绝对产量；绝对容量数据永久缺失。五方向定义与能力边界见 [需求](requirements.md)，分层与产物治理见 [架构](architecture.md)，各方向当前进展见 [执行计划](EXECUTION_PLAN.md)，历史逐次记录见 git。

## 下一步

**仍然数据门控**（注意：不再是"剩余全部"，见上更正）：#2 一致性验证、RNA-seq 表达约束、绝对容量，以及方向 5（碳源条件标定 + 跨条件稳健性）的收尾项。逐条门控原因与现状见 [执行计划](EXECUTION_PLAN.md)「仍然数据门控的部分」。

## 范围边界（硬约束）

- 无界完整跨条件排名产品仍不做；只产相对信号，绝对容量恒 unavailable。
- **glucose 的 corrected_reference 结果不得改动**（回归锁定）。
- **保密湿实验数据只存仓库外本地私有区**（`CursorProject/pcSec_wetlab_private/`）、提交产物只含机制层抽象、不上传云端 / GitHub。
- 方向 4 组合搜索、目标蛋白降解通路建模、换默认 solver：明确不做。

## 必读材料

1. [项目级执行计划](EXECUTION_PLAN.md)：当前状态、待数据工作、范围边界。
2. [需求与能力边界](requirements.md)：五方向定义、能力边界、证据分层、碳源标定契约。
3. [架构与边界](architecture.md)：分层与所有权、能力清单、数据与产物治理。
4. [ADR-006](adr/006-carbon-source-condition-calibration.md)：碳源标定 + 数据契约。
5. [ADR-007](adr/007-secretory-machinery-gene-complex-reachability.md)：分泌机器可达性层（当前 slice 的决策依据）。
6. [ADR-002](adr/002-relative-oe-and-absolute-capacity-layers.md) / [ADR-004](adr/004-relative-signal-deepening-under-permanent-data-gap.md) / [ADR-005](adr/005-rnaseq-expression-constrained-enzyme-capacity.md)。

## 验证方式

- 文档锚点：`tests/test_docs_active_boundary.py`、`test_data_results_boundaries.py`、`test_slow_test_gates.py`。
- 改造求解 / 短名单：glucose 回归结果逐字不变 + 短名单 / 契约测试；碳源标定内部验证。
- 全量回归基线见 git 历史（根 `tests/` 与 `python_pichia/tests/` 隔离全绿）。
