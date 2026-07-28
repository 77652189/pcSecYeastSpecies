# pcSecPichia Handoff

状态：active
最后更新：2026-07-28

## 当前目标

当前 slice：**改造后候选系统（OE + KO · 分层复用）——迭代2 D1–D6 全落地并已全部 push**。约束档开的新鲜全基因组基线也跑完 + 验证（6h / 1025-0err / cached `697fb401`）。**关键结论（生成侧到此为止）**：分层复用机制正确 + 端到端全验证，但**真复用增益不出现**——top-N 短名单被代谢候选占满、代谢按约定恒保守 `已失效`，可复用的分泌专属层候选排不进 top-N（换约束档 / 改造都不变）；用户定不改短名单口径。**#2 实验反馈验证已诊断为数据门控**：在手私有数据是 PH / 温度**条件验证**（时间序列），**没有**“改造→titer”可链数据集，#2 落不了地（机器现成、只差数据）。→ **产品生成侧已完善；剩余有意义工作（#2 验证 / RNA-seq / 绝对容量）全部数据门控，当前收在“等数据”**。详见 [执行计划](EXECUTION_PLAN.md) + [ADR-004](adr/004-relative-signal-deepening-under-permanent-data-gap.md)。

```yaml
current_slice: modified_strain_ko_oe_layered_shortlist
slice_status: done
current_program: mvp_directions_plus_relative_signal_deepening
project_stage: generation_side_complete_holding_for_wetlab_data
absolute_capacity_status: unavailable_waiting_for_qualified_evidence
carbon_source_calibration_status: glucose_corrected_others_internally_calibrated
wetlab_fermentation_data_status: in_hand_is_condition_validation_only_no_modification_titer_screen
experiment_feedback_validation_status: machinery_ready_data_gated_no_modification_titer_dataset
rnaseq_expression_constraint_status: contract_ready_adr005_waiting_for_data
direction_4_combination_status: strategically_deferred
```

一句话背景：模型是“折叠 / 分泌容量这一片”的相对候选生成器，不预测绝对产量；绝对容量数据永久缺失。五方向完整状态、分层架构见 [需求与架构文档](pichia_current_architecture_and_requirements.md)，历史逐次记录见 git。

## 下一步

迭代候选（改造后 OE + KO · 分层复用）**迭代 1（C1–C4）+ 迭代 2（D1–D6）已全部完成并 push**：叠 KO/OE 重解 → 改造后瓶颈 → 分层短名单（L1 打标即时复用 / L2 按需重算），整合进仿真验证页。**诚实结论**：机制正确、端到端全验证，但**真复用增益只在瓶颈集中（约束档开、folding-limited）时兑现**——约束档关的基线下瓶颈弥散、候选全保守标 `已失效`（正确但无增益）。逐 D 细节 + 分类两套词表经 `PROCESS_LABELS` 桥归并的须知见 [执行计划](EXECUTION_PLAN.md) 已完成能力 + git。

生成侧完成后又做了一轮**易用性 / 性能打磨（2026-07-28，已 push）**：筛查 / 仿真页去 dev 术语（研究员向）、结果段 `st.cache_data`、剂量响应曲线 + 明细表 + fact pack 改按需渲染（修首屏卡顿 /“无法滑动”）、移除与 catalog 冗余的 complex_hypothesis。

**剩余全部数据门控**（收在“等数据”）：#2 实验反馈一致性验证（缺“改造→titer”可链集）、RNA-seq 表达约束（[ADR-005](adr/005-rnaseq-expression-constrained-enzyme-capacity.md)）、绝对容量恒 unavailable。方向 5（碳源条件标定 + 跨条件稳健性）收尾项——B5 titer 锚点（待验证数据）、方向 1 本地摄入（经护栏读私有区）、B4 湿实验一致性标注——同为数据门控、暂缓；B3 噪声门控暂不建（0 翻转、无表观敏感可甄别）。**后置（非数据门控）**：D4b 改造后全量 L3 兜底。

## 范围边界（硬约束）

- 无界完整跨条件排名产品仍不做；只产相对信号，绝对容量恒 unavailable。
- **glucose 的 corrected_reference 结果不得改动**（回归锁定）。
- **保密湿实验数据只存仓库外本地私有区**（`CursorProject/pcSec_wetlab_private/`）、提交产物只含机制层抽象、不上传云端 / GitHub。
- 方向 4 组合搜索、目标蛋白降解通路建模、换默认 solver：明确不做。

## 必读材料

1. [项目级执行计划](EXECUTION_PLAN.md)：当前状态、待数据工作、范围边界。
2. [需求与架构文档](pichia_current_architecture_and_requirements.md)：五方向状态、分层架构、碳源标定状态、数据与产物治理。
3. [ADR-006](adr/006-carbon-source-condition-calibration.md)：碳源标定 + 数据契约。
4. [ADR-002](adr/002-relative-oe-and-absolute-capacity-layers.md) / [ADR-004](adr/004-relative-signal-deepening-under-permanent-data-gap.md) / [ADR-005](adr/005-rnaseq-expression-constrained-enzyme-capacity.md)。

## 验证方式

- 文档锚点：`tests/test_docs_active_boundary.py`、`test_data_results_boundaries.py`、`test_slow_test_gates.py`。
- 改造求解 / 短名单：glucose 回归结果逐字不变 + 短名单 / 契约测试；碳源标定内部验证。
- 全量回归基线见 git 历史（根 `tests/` 与 `python_pichia/tests/` 隔离全绿）。
