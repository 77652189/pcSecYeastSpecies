# pcSecPichia 文档索引

状态：active  
最后更新：2026-07-14

项目文档围绕一个原始研发目标组织：通过 KO/OE 和分泌路径改造，提高 hLF、OPN 等目标蛋白在 Pichia 中的分泌表现。

当前系统已经具备 KO/OE 候选生成、证据复核、实验反馈回放、gene-enzyme-reaction capacity 求解，以及外部 baseline capacity 候选的人工导入、审核和显式提升链路，但不能承诺真实发酵产量、mg/L 绝对值或实验成功率。Round 6A 仍在进行：已接入 PRIDE 相对 iBAQ 与 ecPichia supplement 的可审计 source adapter，但绝对 abundance、条件匹配和 `model_flux` 换算链仍未闭合。

## 当前入口

| 文档 | 用途 |
| --- | --- |
| [当前架构与能力边界](pichia_current_architecture_and_requirements.md) | 当前能做什么、不能做什么、核心模块和证据边界 |
| [下一阶段执行计划](pichia_next_plan.md) | 从实验反馈闭环到 gene-level OE、组合筛查和前瞻验证的固定顺序 |
| [数据与结果治理策略](data_and_results_policy.md) | 模型、外部证据、实验数据、LLM 输入和运行产物的保存规则 |
| [当前 handoff](handoff.md) | 当前工作树、真实完成度、下一 checkpoint 和验证方式 |
| [ADR-001：外部容量候选与正式资产分层](adr/001-external-capacity-candidate-promotion.md) | 外部参数如何成为正式容量锚点及其适用范围 |

## 阅读顺序

1. 判断当前能力或讨论产品方向：先读“当前架构与能力边界”。
2. 准备下一 checkpoint：读取“当前 handoff”，只执行其中记录的下一项。
3. 需要理解后续阶段时再读“下一阶段执行计划”，不得重新从 Round 1 开始。
4. 写入数据、模型、缓存或报告前：检查“数据与结果治理策略”。

## 归档边界

`docs/archive/` 保存已经完成的阶段计划、可行性验证、历史设计决策和长排查记录。归档内容保留原始证据，但不再作为当前执行入口。

已归档的主要主题包括：

- COBRApy import、GEM QA 与 Shadow LP 产品化阶段计划。
- BLAST/RBH 同源 crosswalk 详细架构和 SCE/Pichia 可行性试跑。
- 在线外部数据库证据层 Round 1-9 设计与函数契约。
- 旧 KO/OE 全基因组筛查设计、碳源计划和 Python 迁移记录。

归档目录按当前策略不进入公开版本控制；active 文档必须独立表达现状和下一步，不依赖归档文件才能执行。
