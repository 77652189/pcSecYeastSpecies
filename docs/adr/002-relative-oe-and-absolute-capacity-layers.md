# ADR-002：相对 OE 决策层与绝对容量研究层分离

状态：accepted
日期：2026-07-15
关系：补充 ADR-001，不取代其正式容量 promotion 门禁

## 背景

Phase 2 已实现 gene-enzyme-reaction mapping、剂量、约束、求解、报告和 UI，但 A0c 证明现有外部证据无法形成可审核的 baseline capacity。若继续把绝对容量作为所有 OE 能力的唯一入口，已完成的相对比较能力无法服务候选选择；若用 proxy、默认值或求解结果填充容量，则会破坏科学边界。

项目原始目标是降低 KO/OE 候选选择成本，不是预测绝对产量。因此需要把研发决策能力与绝对容量研究能力分开验收。

## 决策

系统保留两个互不替代的 OE 层级：

1. **相对、未校准的 OE 决策层**：允许使用 reaction proxy 和具备明确 mapping、剂量及不确定性的相对 gene-capacity 场景进行候选比较、风险解释和实验优先级排序。
2. **绝对 gene-capacity 研究层**：只有存在符合 ADR-001 的审核 baseline capacity 时才可执行；否则状态必须为 unavailable/not executable。

相对层不得被命名或展示成绝对 capacity。绝对层不可用时，不得静默降级到相对层；调用方必须显式选择或展示执行模式和校准状态。

实验反馈只校准排序、方向一致性和风险判断。它不能自动生成 baseline capacity，也不能直接修改代谢矩阵、GPR、curated phenotype 或正式科学资产。

## 不变量

- reaction proxy、相对 gene-capacity 和绝对 gene-capacity 是三个可区分状态。
- 没有审核 baseline capacity 时，绝对层保持 unavailable。
- 不使用通用上界、baseline optimal flux、固定 `1.0`、fixture、相对蛋白组强度或未审核外部模型值伪造绝对容量。
- 相对层输出必须保留来源、mapping、剂量、uncertainty、warning 和不可解释范围。
- 相对层不能输出 mg/L、真实表达倍数、实验成功概率或跨条件保证。
- 新实验数据通过独立反馈层进入校准，不覆盖原始预测和科学资产。

## 验收边界

相对决策层完成至少要求：

- 用户可区分 reaction proxy、relative uncalibrated 和 absolute unavailable。
- 旧 proxy 数值与 feature-off 行为保持回归。
- 缺 baseline capacity 时不会调用绝对容量求解或产生 nominal capacity。
- hLF/OPN 输出分别保留 target、context、模型版本和风险说明。
- service/UI 不重新实现模式判断或科学降级逻辑。

绝对容量研究层继续使用 ADR-001 的 source、evaluation、review 和 promotion 门禁，不因相对层可用而视为通过。

## 影响

该决策使项目可以继续利用已完成的模型和证据能力支持实验候选选择，同时诚实保留绝对容量缺口。代价是产品和报告必须长期维护清晰的模式、校准状态和禁用承诺，不能用单一“OE 结果”概括不同科学含义。
