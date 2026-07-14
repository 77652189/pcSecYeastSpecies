# pcSecPichia 当前架构与能力边界

状态：active  
最后更新：2026-07-14

## 原始研发目标

项目要解决的问题是：围绕 hLF、OPN 等目标蛋白，通过 KO、OE 和分泌路径改造，找到更可能提高 Pichia 分泌表现的工程方案，并降低研发同事阅读模型结果和选择实验候选的成本。

当前项目的正确定位是“KO/OE 候选生成、模型验证和证据排序系统”，不是实验产量预测器。

## 当前结论

当前系统可以回答：

- 哪些单基因 KO 在模型 GPR 中可执行。
- 哪些 OE 方向在反应容量层面可能改善目标蛋白分泌。
- 候选对生长、分泌和模型可行性的相对影响。
- 候选是否存在 essentiality、growth risk、proxy、同源映射或外部证据风险。
- hLF 和 OPN 是否出现不同的候选优先级。
- 哪些候选值得进入人工复核和小规模实验。

当前系统不能可靠回答：

- 某个改造会提高多少 mg/L。
- 某个候选的真实实验成功概率。
- reaction-level OE proxy 等同于真实 gene-level OE 后的表达变化。
- 模型外基因、调控网络和全部分泌通路基因的定量作用。
- 多基因组合的上位性和长期培养稳定性。

## 研发工作流

```text
目标蛋白与培养条件
  -> pcSecPichia / Shadow LP 模型准备
  -> 单基因 KO 与 OE proxy 筛查
  -> 分泌、生长、可行性和资源代价计算
  -> gene/GPR/同源/外部数据库/表型证据合并
  -> recommendation tier 与风险分层
  -> fact pack + 程序 validator + LLM writer + Judge
  -> hLF / OPN 实验候选清单
  -> 实验反馈与下一轮排序校准
```

实验反馈与排序校准的结构化闭环已经完成；gene-level OE capacity 的契约、约束、求解、报告、Streamlit 和正式门禁也已实现。当前核心缺口是缺少独立于求解结果、带来源和适用条件的 baseline capacity 锚点，因此真实 hLF/OPN gene-capacity 仍保持不可执行。

## 已有能力

### 模型与筛查

- 加载 pcSecPichia 参考模型和目标蛋白输入。
- 构建内置目标和 custom target 的分泌路径。
- 设置培养基、碳源、固定生长率和目标蛋白负担。
- 执行 genome-wide / catalog-level 单基因 KO 筛查。
- 生成 OE reaction-capacity proxy，并保留“不是 gene-level OE”的 warning。
- 输出 secretion ratio、growth retention、最大可行生长率、protein cost 和 solver 状态。
- 对 essential KO、growth risk、不可解析和模型不可执行候选进行降级。

### 基因、GPR 与证据

- 模型 gene_index、蛋白序列和 UniProt/NCBI/KEGG 等命名注释可交叉查询。
- BLAST/RBH cache 可提供 SCE 到 Pichia 的同源候选、name audit 和 rule-transfer audit。
- 外部数据库 fetcher 可受控生成 UniProt、NCBI、SGD 和外部 GEM/GPR 证据缓存。
- 外部 GPR 只作为候选证据；映射到当前 reaction/gene 后才可成为模型可执行规则。
- phenotype evidence、数据库注释、模型 GPR、OE proxy 和同源证据保持分层。

### 求解与质量验证

- 当前大模型主路径使用结构化约束和 SciPy HiGHS / reference solve。
- Shadow LP 已承载 pcSec 约束并与内置 hLF/OPN reference 结果高精度对齐。
- COBRApy 用于外部 GEM import、基础 QA 和小模型语义验证，不作为 full pcSec 默认后端。
- Streamlit 可手动运行和读取 Shadow LP cross-check 报告。

### 结果解释

- KO/OE preview、screen rows、recommendation、历史结果和 Streamlit 均可透传标准命名与证据字段。
- LLM 只读取程序生成的 fact pack。
- 程序 validator 检查 schema、evidence_id、target、gene/reaction、数值、KO/OE 分组和最低信息覆盖。
- Judge 读取 fact pack、writer 输出和 validator 结果，不合格时反馈重写。
- 最终报告按 hLF / OPN 分区，引用表由程序生成。

### 实验反馈

- CSV/XLSX/JSONL 实验记录可导入统一 schema，并保留原始值、单位、重复、失败和排除状态。
- prediction linkage 显式区分 matched、ambiguous、missing 和 context mismatch。
- hLF/OPN 回放可生成校准记录、方向一致性、Top-K 描述指标和证据不足提示。
- 实验导入不会自动修改 recommendation tier、curated evidence 或模型约束。

## 核心证据边界

| 层级 | 能说明什么 | 不能说明什么 |
| --- | --- | --- |
| 数据库注释 | 名称、功能和来源一致性 | 真实 KO/OE 表型 |
| BLAST/RBH | 序列层同源候选 | 当前 GEM 可执行性或功能完全等价 |
| 模型 GPR | KO 可以在当前模型中执行 | 实验一定可行 |
| OE reaction proxy | 关联反应容量变化方向 | 启动子、拷贝数和 gene-level 表达量 |
| gene-level enzyme capacity | 在明确映射和参数场景下比较某基因 OE 的模型内容量与资源代价 | 真实表达倍数、绝对产量或跨条件保证 |
| curated phenotype | 特定 intervention/context 的已有表型方向 | 通用 mg/L 或跨条件保证 |
| 实验反馈 | 当前宿主、目标和条件下的观测结果 | 未测试条件的自动外推 |

`experiment_calibrated` 只表示存在严格匹配的高置信表型或内部实验支持，不表示能够预测绝对产量。

## 当前主要缺口

1. gene-level OE capacity 链路已经实现，但缺少可审核的 baseline capacity 锚点，真实 hLF/OPN 候选尚不能通过正式验收。
2. Phase 1 尚未接入获批真实实验数据，当前只有脱敏回放和数据契约证据。
3. 分泌、折叠、UPR、ERAD、糖基化和囊泡运输的机制约束仍不完整。
4. 大量模型外蛋白有序列和注释，但没有可执行 GPR，不能直接进入 KO/OE 求解。
5. 筛查以单基因为主，缺少组合改造和上位性分析。
6. 缺少跨碳源、生长率、氧供和参数扰动的鲁棒性排名。
7. 尚未用前瞻实验评价 top-K 命中率、富集倍数和排序校准程度。

## 架构边界

```text
app/ui
  -> app/services
    -> python_pichia
       -> loading / media / targets
       -> constraints / simulation / shadow_lp
       -> screens / analysis / reports
       -> services / homology / external_refs / experimental_feedback
       -> oe_capacity
```

- `python_pichia/`：核心科学逻辑、数据契约、约束、求解、筛查、证据和报告事实层。
- `app/services/`：facade、任务触发、缓存路径和错误汇总，不实现科学判断。
- `app/ui/`：Streamlit 展示和用户操作，不直接修改模型语义。
- `Code/`、`Model/`、`Enzymedata/`、`Results/`：legacy/reference 科学资产，默认只读。
- `local_runs/`：运行产物、缓存、报告和验证证据，默认 ignored。

### Phase 2 所有权边界

- `oe_capacity` 负责 gene-enzyme-reaction 映射、OE 剂量、参数不确定性、约束计划和 proxy 对照结果。
- `screens/gene_interventions.py` 继续拥有现有 GPR 解析与兼容 proxy 规划；Phase 2 调用它，不复制另一套 GPR 解释器。
- `core/pichia_enzymes.py` 和当前模型中的酶变量/稀释反应是本地执行语义的首要来源。
- `external_refs` 只提供外部 GEM、GPR、kcat、丰度和同源映射候选；未经当前模型 reaction/gene/enzyme 复核不得成为可执行约束。
- `shadow_lp`、SciPy HiGHS 和现有 reference solve 仍是求解路径；Phase 2 不引入新的默认 solver。
- gene-level capacity 与旧 reaction proxy 必须作为两个明确模式并行存在，不能通过改名掩盖降级。

### 外部容量候选与正式资产边界

研发组无法提供内部 baseline capacity 时，系统允许从外部来源建立候选，但外部候选不能直接成为可执行约束。来源优先级为：

1. 同宿主、同菌株且培养条件接近的 Pichia 定量蛋白组。
2. 经版本和许可审计的 iPichia/ecPichia 等酶约束模型。
3. Pichia 文献以及 BRENDA、SABIO-RK 等动力学来源与蛋白组的组合换算。
4. 通过 BLAST/RBH 确认的 S. cerevisiae 同源参数转移，仅作为低置信区间。

外部参数分为四种适用范围：

- `target_specific`：同宿主、同条件且针对当前目标蛋白获得的参数。
- `host_condition`：宿主基础容量，可在宿主、菌株、培养基、碳源和生长状态匹配时供 hLF/OPN 复用。
- `external_model_calibrated`：来自外部酶约束模型的校准参数，不等同于实测丰度。
- `homolog_transferred`：由其他物种同源转移，只能生成宽不确定性区间。

匹配优先级固定为 `target_specific > host_condition > external_model_calibrated > homolog_transferred`。候选必须完成单位换算、当前模型 gene/enzyme/formation 映射、条件匹配、来源版本/hash/license、冲突检查和人工审核后，才可提升到 `Enzymedata/oe_capacity_baseline_capacity.json`。`1000` 通用上界、baseline optimal flux、fixture 和未审核同源参数永远不能作为正式锚点。

长期决策见 [ADR-001](adr/001-external-capacity-candidate-promotion.md)。

## 项目成功标准

项目是否更接近原始目标，不以新增多少字段或页面衡量，而以这些结果衡量：

- top-K 候选相对随机或人工经验基线是否有更高实验命中率。
- hLF 与 OPN 的候选排序是否能在重复实验中复现。
- 模型预测方向、实验方向和不确定性是否得到校准。
- 生长风险、不可执行候选和 proxy 结论是否被正确拦截。
- 每条实验建议是否能追溯到模型结果和证据来源。

下一步固定从[下一阶段执行计划](pichia_next_plan.md)中当前未完成轮次继续。
