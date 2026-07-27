# pcSecPichia Handoff

状态：active
最后更新：2026-07-27

## 当前目标

当前 slice：**改造后候选系统（OE + KO · 分层复用）——迭代2 D1–D6 全落地并全 push**（origin/main=`46b5a52`，5 笔）。约束档开的新鲜全基因组基线也跑完+验证（6h / 1025-0err / cached `697fb401`）。**关键结论（生成侧到此为止）**：分层复用机制正确+端到端全验证，但**真复用增益不出现**——top-N 短名单被代谢候选占满、代谢按约定恒保守 `已失效`，可复用的分泌专属层候选排不进 top-N（换约束档/改造都不变）；用户定不改短名单口径。**#2 实验反馈验证已诊断为数据门控**：在手私有数据是 PH/温度**条件验证**（时间序列），**没有**"改造→titer"可链数据集，#2 落不了地（机器现成、只差数据）。→ **产品生成侧已完善；剩余有意义工作（#2 验证 / RNA-seq / 绝对容量）全部数据门控，当前收在"等数据"**。详见 [执行计划](EXECUTION_PLAN.md) + [ADR-004](adr/004-relative-signal-deepening-under-permanent-data-gap.md)。

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

一句话背景：模型是"折叠/分泌容量这一片"的相对候选生成器，不预测绝对产量；绝对容量数据永久缺失。五方向完整状态、分层架构见 [需求与架构文档](pichia_current_architecture_and_requirements.md)，历史逐次记录见 git。

## 下一步

阶段② 迭代1（C1–C4，下一步 OE 候选）✅ **已完成**（本地 `d2b4cd0` 未 push）：`strain_modifications` 叠 KO/OE + `solve_secretion_capacity`/`run_oe_dose_response_sweep` opt-in 参数（默认 None → glucose 逐字不变）；引擎 `next_oe_candidates` + 服务 `per_strain_oe_candidate_run` + 仿真验证页 UI；586/342 全绿、app 实跑通过。

迭代2（改造后候选系统 · 分层复用 + 全基因组后台，用户 2026-07-24 拍板；清单见执行计划阶段② 迭代2）：
- **D1 KO 候选引擎** ✅（`8e502e7`）：`run_knockout_screen` opt-in `strain_modifications` → 改造后基线跑 KO 扰动。
- **D4 复用地基** ✅（本地，待 push）：野生型**口径指纹**基线缓存——引擎 `strain_baseline_cache`（口径+菌株指纹、复用 `solve_cache`/`model_variant_fingerprint`；schema post-`87f99ac`、旧 unknown 分类永不命中）+ 服务 `strain_baseline_service`（跟随 run 口径读/报未构建/CSV ingest）+ 并行筛查工具加口径参数&完成后蒸馏进缓存。8+4 单测、2 基因真跑闭环、2050 行真 CSV 蒸馏都验过。
- **D2/D3/D5** ✅（`00b6f70`+`3d236ac`）：D2 打标（`per_strain_layer_reuse`）+ D3 L1 编排 + L2 重算（`per_strain_shortlist_run` + 引擎 `recompute_modified_strain_candidate_effects`）+ D5 面板（`simulation_results._render_modified_strain_shortlist`，复用 C3 stash）。L2 引擎真跑 smoke 验过。
- **D6** ✅（`34f1a5c`+`46b5a52`）：路2 live demo（ingest overnight CSV 真跑 L1→L2 机制通）+ 两处实跑修（`top_n=0` 跳 L1 冗余 sweep、效应显示 `:.3g`）。
- **收尾状态（2026-07-27，全 push origin=`46b5a52`）**：① 约束档开新鲜基线已跑（`hlf_folding_baseline`、6h/1025-0err、cached `697fb401`）并验证——**真复用增益仍不出现**（top-N 全代谢→恒保守 `已失效`，见上"关键结论"）；② **#2 已诊断 = 数据门控**：在手私有数据是 PH/温度条件验证（时间序列）、**无"改造→titer"可链集**，#2 落不了地——**机器现成、收在"等数据"**（数据到位当天可落）；③ **D4b** 改造后全量 L3（兜底、后置）；④ **项目当前收在"生成侧完成、等数据"**：剩余有意义工作（#2 验证 / RNA-seq / 绝对容量）全部数据门控。
- 核心前提（诚实）：复用是**近似**（同层/紧邻才需重算），须显式标注失效范围；**KO 无免费派生、必求解**（OE 有影子价格捷径）；层复用只对**分泌专属层**干净、代谢桶保守（分类两套词表经 `PROCESS_LABELS` 桥归并粗模块）。
- **D2 起点须知（D4 时坐实的两处更正）**：① 分类**两套词表不同**——LP 瓶颈用 `classify_secretory_process`（英文键）、筛查短名单用 `gene_perturbation_map`（中文展示标签），**不能直接比对**；D2 要把候选的 `affected_reactions` 经 `classify_secretory_process` 重分类成同词表再比（基线已保留 `affected_reactions`）。② 分类只有 5 个粗桶、绝大多数落 `metabolic_or_other`——层级复用对**分泌专属层**（折叠/糖基化/ER 转运/ERAD/Golgi）才干净有效，代谢桶要保守。

方向5（碳源条件标定 + 跨条件稳健性）收尾（数据门控、暂缓）：B5 titer 锚点（待验证数据）+ 方向1 本地摄入（经护栏读私有区）、B4 湿实验一致性标注。B3 噪声门控暂不建（0 翻转、无表观敏感可甄别）。

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
