# pcSecPichia ADR 索引

状态：active

| ADR | 状态 | 决策 |
| --- | --- | --- |
| [ADR-001](001-external-capacity-candidate-promotion.md) | accepted | 外部容量候选如何经过审核提升为正式科学资产 |
| [ADR-002](002-relative-oe-and-absolute-capacity-layers.md) | accepted | 相对 OE 决策层与绝对容量研究层保持独立验收 |
| [ADR-003](003-fermentation-feedback-minimal-fields.md) | accepted | 实验反馈发酵模板以现场真实 14 字段结构为基准，只补最小必要缺口 |
| [ADR-004](004-relative-signal-deepening-under-permanent-data-gap.md) | accepted | 承认绝对 OE/容量数据永久缺失，授权相对层四项免数据信号（影子价格瓶颈归因、OE 剂量响应、排序对容量假设的稳健性、价值-of-information） |
| [ADR-005](005-rnaseq-expression-constrained-enzyme-capacity.md) | accepted | RNA-seq 表达约束的菌株特异建模数据契约（transcript→酶丰度上界，经 curated 基因→复合体映射，相对/opt-in，绝对恒 unavailable，实现待数据） |
| [ADR-006](006-carbon-source-condition-calibration.md) | accepted | 碳源条件标定与三档状态（corrected_reference / internally_calibrated / draft_boundary）+ 升 corrected 的数据契约 |
| [ADR-007](007-secretory-machinery-gene-complex-reachability.md) | accepted | 分泌机器「基因 ↔ 复合体」映射作为可达性层：界面不把基因/反应分裂交给用户，映射用于扩大筛查覆盖 + 把复合体过表达翻译成实验基因；可达性≠准确度，绝对容量恒 unavailable；明确不给代谢基因补外部 GPR |

新增 ADR 取代旧决策时，必须在新旧文档中记录替代关系；未声明替代时视为互补决策。
