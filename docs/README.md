# pcSecPichia 文档索引

状态：active  
最后更新：2026-07-14

项目文档围绕一个原始研发目标组织：通过 KO/OE 和分泌路径改造，提高 hLF、OPN 等目标蛋白在 Pichia 中的分泌表现。

当前系统已经具备 KO/OE 候选生成、证据复核、实验反馈回放和 gene-enzyme-reaction capacity 求解链路，但不能承诺真实发酵产量、mg/L 绝对值或实验成功率。当前重点是从外部 Pichia 蛋白组、酶约束模型和动力学数据库建立可审核的 baseline capacity 候选，经人工提升后完成 hLF/OPN 正式验收。

## 当前入口

| 文档 | 用途 |
| --- | --- |
| [当前架构与能力边界](pichia_current_architecture_and_requirements.md) | 当前能做什么、不能做什么、核心模块和证据边界 |
| [下一阶段执行计划](pichia_next_plan.md) | 从实验反馈闭环到 gene-level OE、组合筛查和前瞻验证的固定顺序 |
| [数据与结果治理策略](data_and_results_policy.md) | 模型、外部证据、实验数据、LLM 输入和运行产物的保存规则 |
| [ADR-001：外部容量候选与正式资产分层](adr/001-external-capacity-candidate-promotion.md) | 外部参数如何成为正式容量锚点及其适用范围 |

## 阅读顺序

1. 判断当前能力或讨论产品方向：先读“当前架构与能力边界”。
2. 准备下一轮开发：再读“下一阶段执行计划”，从当前未完成轮次继续，不重新从 Round 1 开始。
3. 写入数据、模型、缓存或报告前：检查“数据与结果治理策略”。

## 归档边界

`docs/archive/` 保存已经完成的阶段计划、可行性验证、历史设计决策和长排查记录。归档内容保留原始证据，但不再作为当前执行入口。

已归档的主要主题包括：

- COBRApy import、GEM QA 与 Shadow LP 产品化阶段计划。
- BLAST/RBH 同源 crosswalk 详细架构和 SCE/Pichia 可行性试跑。
- 在线外部数据库证据层 Round 1-9 设计与函数契约。
- 旧 KO/OE 全基因组筛查设计、碳源计划和 Python 迁移记录。

归档目录按当前策略不进入公开版本控制；active 文档必须独立表达现状和下一步，不依赖归档文件才能执行。
