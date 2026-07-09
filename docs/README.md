# pcSecPichia 文档索引

状态：active  
最后更新：2026-07-09

当前文档只保留对现状判断和下一步执行有用的入口。历史计划、阶段验证和长排查记录放在 `docs/archive/`。

## 当前入口

| 文档 | 用途 |
| --- | --- |
| [当前架构与需求](pichia_current_architecture_and_requirements.md) | 当前工作目标、系统分层、能力边界和 COBRApy / Shadow LP 状态 |
| [下一步计划](pichia_next_plan.md) | 接下来应做什么、暂不做什么、每项的验证命令 |
| [BLAST/RBH 同源映射架构](pichia_homology_crosswalk_architecture.md) | 酿酒酵母到 Pichia 的本地同源证据层、BLAST/RBH cache、name audit、rule-transfer audit 和 Streamlit 审计边界 |
| [在线外部数据库证据层架构](pichia_online_external_reference_architecture.md) | UniProt / NCBI / SGD 的受控联网 fetcher、external reference cache、命名校对和函数契约 |
| [数据与结果治理策略](data_and_results_policy.md) | 哪些目录只读、哪些产物进 `local_runs/`、何时需要 checkpoint |

## 已归档内容

历史计划、阶段验证和长排查记录保留在本地 `docs/archive/`，不纳入版本控制（包含项目内部设计细节，不适合公开分发）。

## 当前推荐阅读顺序

1. 先读 [当前架构与需求](pichia_current_architecture_and_requirements.md)。
2. 再读 [下一步计划](pichia_next_plan.md)。
3. 若要做 BLAST/RBH 或 Streamlit 同源审计，读 [同源映射架构](pichia_homology_crosswalk_architecture.md)。
4. 若要接入 UniProt / NCBI / SGD 或扩展新的外部来源，读 [在线外部数据库证据层架构](pichia_online_external_reference_architecture.md)。
5. 执行前确认 [数据与结果治理策略](data_and_results_policy.md)。
