# pcSecPichia 项目级执行计划

状态：active  
最后更新：2026-07-28

> 本文件只列**当前状态与待启动 / 待数据的工作**。已完成阶段不在此展开——其形成的能力见 [需求与架构文档](pichia_current_architecture_and_requirements.md)，逐次历史见 git。五方向的完整定义与状态、分层架构、数据与产物治理同样以需求与架构文档为准；技术计划不得绕过本文件扩大范围。

## 项目目标

围绕 hLF、OPN 等目标蛋白，持续生成、解释、校准 KO/OE/分泌通路候选，降低实验候选选择成本。**不是绝对产量预测器**；绝对 OE/容量数据经三次独立调研确认公开来源永久缺失，主力放在相对决策层（[ADR-002](adr/002-relative-oe-and-absolute-capacity-layers.md) / [ADR-004](adr/004-relative-signal-deepening-under-permanent-data-gap.md)）。湿实验数据不定期到来，经独立校准层回填，不直接改代谢矩阵 / GPR / 正式科学资产。

## 当前状态：生成侧完成，收在“等数据”

到 2026-07-28，模型的**候选生成侧已完整**：单基因 KO/OE 筛查、复合体级分泌机器干预、相对信号深化 R1–R4、碳源条件标定（内部档）、改造后候选系统（叠 KO/OE 重解 + 分层复用短名单）都已落地、测试并 push。

**剩余有意义的工作全部数据门控**——不是没建、是没数据（机器现成，数据到位当天可落）：

| 待启动工作 | 门控 | 现状 |
| --- | --- | --- |
| #2 实验反馈一致性验证（模型 KO/OE 预测 ↔ 实验改造→titer 对齐记分卡） | 缺“改造→titer”可链数据集 | 机器现成（`direction_consistent` / `rank_correlation` 信号已在）；在手私有数据只有 PH / 温度条件验证时间序列，链不上 |
| RNA-seq 表达约束建模（阶段③） | 缺生产菌株 RNA-seq | 契约与方法已定，见 [ADR-005](adr/005-rnaseq-expression-constrained-enzyme-capacity.md)；数据到位后 transcript→酶丰度上界→经 curated 基因→复合体映射触达分泌层，相对 / opt-in |
| 碳源升 `corrected_reference` + B5 titer 锚点 + 方向 1 本地摄入 | 缺各碳源定量 / validated titer | hLF μ 已核（默认 0.10 与甘油生长相一致），甲醇 / OPN 定量仍缺；数据契约见 [ADR-006](adr/006-carbon-source-condition-calibration.md) |
| 绝对 gene-level 容量层 | 缺 reviewed baseline capacity（公开源永久缺失） | 恒 `unavailable`，按 [ADR-001](adr/001-external-capacity-candidate-promotion.md) 门禁 |

下一实际切片只从 [handoff](handoff.md) 继续；数据到来前**不再堆生成侧新功能**。

## 已完成能力（细节见架构文档，不在此展开）

- **相对信号深化 R1–R4**（横切层，[ADR-004](adr/004-relative-signal-deepening-under-permanent-data-gap.md)）：影子价格瓶颈归因（含 `bound_type` 过滤，下界不报成 OE 瓶颈）、OE 剂量响应形状、排序对容量假设的稳健性、价值-of-information。均只产相对信号，绝对层恒 unavailable。
- **碳源条件标定 + 短名单跨条件稳健性**（原方向 5 有界升级，[ADR-006](adr/006-carbon-source-condition-calibration.md)）：五碳源 bound / 生长反应 / 蛋白含量条件化到 `internally_calibrated`（glucose 逐字不变）；求解结果内容寻址缓存；短名单跨条件矩阵 + 稳健性标注。**结论**：hLF OE 短名单跨碳源（甘油 / 葡萄糖）top-15 全 `saturating`、排序逐字不变（分泌机器瓶颈与碳源无关，合 folding-limited）；真实工艺条件集 0 表观敏感，噪声门控（B3）留待未来出现敏感时再建。
- **改造后候选系统**（迭代 1 + 迭代 2，[ADR-004](adr/004-relative-signal-deepening-under-permanent-data-gap.md)）：野生型短名单绑全基因筛查，改造过的菌株要“再找瓶颈、排下一候选”。走真·改造后重解（叠 KO/OE 再 solve，瓶颈随改造转移）。迭代 1 = 下一步 OE 候选；迭代 2 = OE + KO 分层复用（受影响层重算、其余复用同口径野生型基线缓存）+ 全基因组后台。**诚实结论**：机制正确、端到端全验证，但**真复用增益只在瓶颈集中（约束档开、folding-limited）时兑现**——约束档关的基线下瓶颈弥散于所有层，候选被保守全标“已失效”（正确但无增益）。KO 无免费派生、必求解；层复用只对分泌专属层干净、代谢桶保守。
- **易用性 / 性能打磨**（2026-07-28）：筛查 / 仿真页去 dev 术语（研究员向）、结果段 `st.cache_data`、剂量响应曲线 + 明细表 + fact pack 改按需渲染（修首屏卡顿 / “无法滑动”），移除与 catalog 剂量响应冗余的 complex_hypothesis。

## 后置（非数据门控、低优先）

- **D4b 改造后全量 L3（兜底）**：分层复用不可信时，以改造后为基线全基因组重跑 + 按菌株指纹缓存（需把 `strain_modifications` 接进 `genome_wide_tradeoff` 路径，hour-scale）。仅作 escape hatch，L1 / L2 分层复用是默认。

## 明确不做

- **方向 4 组合 / 多基因搜索**（含 GA / SA / MILP）：真实组合改造在模型范围之外，模型内搜遗传组合低价值；留到有稳定性标注的可信排序后，仅在 top 短名单做有界两两上位性（O(k²)、小 k、仍相对）。
- **目标蛋白降解通路（PEP4 / PRB1 / YPS 家族）建模**：基因身份低置信度待复核、无合理动力学路径；等真实湿实验结果，不是“以后再看”。
- 伪造绝对容量（通用上界 / 最优 flux / 固定 1.0 / fixture）、改写受保护科学资产、引入新默认 solver。
- 固定湿实验 pilot / 预留湿实验预算；EVO2 / GPU / 云端推理或占位接口。

## 硬边界

无界完整跨条件排名产品**仍不做**；只产相对信号、不产绝对容量；**不动 glucose 的 `corrected_reference` 基准**（回归锁定）；不换默认 solver；**保密湿实验数据只存仓库外本地私有区**（`CursorProject/pcSec_wetlab_private/`），提交产物只含机制层抽象、不上传云端 / GitHub。

## 文档治理

- **本执行计划**：只列当前状态与待启动 / 待数据工作，剔除已完成阶段细节。
- **[需求与架构文档](pichia_current_architecture_and_requirements.md)**：五方向完整定义与状态、分层架构、能力边界、数据与产物治理。
- **ADR**：长期高影响决策（分层 ADR-002、相对信号深化 ADR-004、RNA-seq 数据契约 ADR-005、碳源标定数据契约 ADR-006 等）。
- **handoff**：当前目标、下一步、必读材料、验证方式。
