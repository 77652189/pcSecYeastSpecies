# pcSecPichia 文档索引

状态：active  
最后更新：2026-07-07

当前文档只保留对现状判断和下一步执行有用的入口。历史计划、阶段验证和长排查记录放在 `docs/archive/`。

## 当前入口

| 文档 | 用途 |
| --- | --- |
| [当前架构与需求](pichia_current_architecture_and_requirements.md) | 当前工作目标、系统分层、能力边界和 COBRApy / Shadow LP 状态 |
| [下一步计划](pichia_next_plan.md) | 接下来应做什么、暂不做什么、每项的验证命令 |
| [BLAST/RBH 同源映射架构](pichia_homology_crosswalk_architecture.md) | 酿酒酵母到 Pichia 的离线同源证据层设计 |
| [数据与结果治理策略](data_and_results_policy.md) | 哪些目录只读、哪些产物进 `local_runs/`、何时需要 checkpoint |

## 已归档内容

| 归档文档 | 原因 |
| --- | --- |
| [COBRApy Phase 0 baseline](archive/cobrapy_phase0_baseline_assessment_2026-07-06.md) | 阶段性验证已完成，当前结论已合并到架构文档 |
| [COBRApy Phase 3 installed validation](archive/cobrapy_phase3_installed_shadow_validation_2026-07-06.md) | 阶段性验证已完成，当前结论已合并到架构文档 |
| [KO/OE genome screen design](archive/pichia_ko_oe_genome_screen_design_2026-07-02.md) | 长设计和历史排查记录已不适合作为当前入口 |

## 当前推荐阅读顺序

1. 先读 [当前架构与需求](pichia_current_architecture_and_requirements.md)。
2. 再读 [下一步计划](pichia_next_plan.md)。
3. 若要做 BLAST/RBH，读 [同源映射架构](pichia_homology_crosswalk_architecture.md)。
4. 执行前确认 [数据与结果治理策略](data_and_results_policy.md)。
