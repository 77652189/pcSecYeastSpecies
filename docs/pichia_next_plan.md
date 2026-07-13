# pcSecPichia 下一阶段执行计划

状态：active  
最后更新：2026-07-13

## 当前执行位置

```yaml
current_phase: phase_2_gene_level_oe
current_round: round_0_architecture
round_status: ready
```

执行会话必须从这里记录的阶段和轮次继续。完成一轮后，只更新到下一轮；不得在新会话中自动重置为 Phase 1 Round 0。

Phase 1 Round 0-5 全部完成并验收后，将状态更新为：

```yaml
current_phase: phase_2_gene_level_oe
current_round: round_0_architecture
round_status: ready
```

进入 Phase 2 前必须先完成 Phase 1 的端到端验收和 checkpoint，不能在同一轮中越过阶段边界。

## 总目标

把现有 KO/OE 候选生成系统逐步升级为能够利用实验反馈、提高候选命中率的研发系统。

执行顺序固定为：

1. 实验反馈闭环。
2. gene-level OE 与酶容量语义。
3. 分泌资源和蛋白稳态约束。
4. 组合改造与条件鲁棒性。
5. 前瞻实验验证和排序校准。

不得跳过前一阶段的数据契约和验收，直接选择后面最容易实现的页面、字段或报告功能。

## 已完成基线

以下能力已经完成，不再重复立项：

- 单基因 KO、OE reaction proxy、growth tradeoff 和 recommendation tier。
- 基因命名标准化、BLAST/RBH、外部数据库证据和外部 GPR 候选审计。
- Shadow LP 约束求解、reference cross-check、COBRApy import 和 GEM basic QA。
- Streamlit 筛查、仿真、同源审计、cross-check 和 LLM judged report。
- 数据产物写入 `local_runs/`，保护 MATLAB/reference 资产目录。

历史实现轮次和详细函数契约已经归档；后续执行只能补当前阶段缺口。

## 固定执行循环

每个 Round 必须按以下顺序完成：

1. 读取 `docs/README.md`、当前架构、当前计划和数据治理策略。
2. 确认本轮输入、输出、数据边界和验收标准。
3. 先实现核心 `python_pichia` 数据契约或科学逻辑。
4. 增加 focused tests 和 fixture。
5. 再接 `app/services` facade；只有用户工作流需要时才接 Streamlit。
6. 运行 smoke，检查产物进入 ignored `local_runs/`。
7. review、自动修复、复验，最多 3 轮。
8. 更新当前状态；不得重写已完成阶段的大段设计。
9. 检查 `Code/Model/Enzymedata/Results` 和依赖声明没有非预期 diff。

## 子 agent 协作策略

Phase 1 可以使用内置子 agent 提高速度，但所有子任务仍属于同一个目标和同一个主分支。不要再创建多个用户会话并手动切换共享工作区分支。

### 主 agent 职责

- 持有唯一 active goal、架构决策、当前 Round 顺序和验收标准。
- 处理当前关键路径和会阻塞其他工作的接口定义。
- 冻结子任务之间的 API 和文件所有权后再并行派发。
- 审阅、整合所有子 agent 结果，解决跨模块冲突。
- 独占 stage、commit、push、文档状态更新和目标完成判定。

### 子 agent 规则

- 优先使用 1-3 个子 agent；没有真正并行工作时不为“形式上并行”创建 agent。
- explorer 只做具体、可回答的代码发现或影响分析，不修改文件。
- worker 必须获得明确的文件所有权、输入契约、输出契约和 focused test。
- 两个 worker 不得编辑同一个文件、同一 fixture 或同一测试模块。
- worker 在自己的 forked workspace 中工作，不切换主工作区分支，不 commit、不 push。
- worker 必须说明修改文件、验证命令、已知限制和未解决问题。
- 主 agent 不重复实现已经委派的任务；等待期间处理不重叠的关键路径工作。
- 子 agent 产物整合后必须由主 agent 重新 review 和运行测试，不能直接视为完成。
- smoke 输出使用独立的 tmp_path 或唯一 `local_runs/experiment_feedback/<round>/<agent>/`，避免相互覆盖。

### 适合并行的任务

- 多个互不依赖的现有代码路径审计。
- 核心接口冻结后的独立实现与独立测试/fixture。
- service 和 UI 在 facade contract 已冻结后的并行开发。
- hLF、OPN 独立 smoke 和只读安全审计。

### 不适合并行的任务

- Phase 1 canonical schema 和 public API 的最终决定。
- 会影响所有后续模块的 ID、单位、control matching 和 calibration eligibility。
- 同一个文件或同一个测试 fixture 的修改。
- stage、commit、push、迁移稳定科学资产或修改保护目录。

### 各 Round 推荐分工

| Round | 主 agent 关键路径 | 可委派子任务 |
| --- | --- | --- |
| 0 | 冻结 schema、public API、status 和文件边界 | explorer 审计现有 schema/cache 模式；worker 在 API 冻结后补独立 contract tests/fixtures |
| 1 | 确认 bundle/manifest 语义并整合 | worker A 实现 IO；worker B 实现独立坏输入 fixtures/tests，文件互斥 |
| 2 | 冻结 PredictionIndex 和 linkage status | explorer 追踪 screen/report 字段；worker 实现 linkage tests 或 prediction adapter |
| 3 | 决定 control matching 和 calibration 公式 | worker 分别实现统计测试与报告序列化；不得并行设计两套公式 |
| 4 | 冻结 service contract 和页面工作流 | worker A 负责 service+service tests；worker B 负责 UI+UI tests，不能修改对方文件 |
| 5 | 汇总验收和最终决策 | worker A 跑 hLF smoke；worker B 跑 OPN smoke；explorer/worker C 做只读数据泄漏和保护目录审计 |

每轮最多进行一次主并行波次。若子任务结果暴露接口问题，先由主 agent 收束接口，再决定是否进行第二波；不能无限创建 agent 掩盖架构不清。

## Phase 1：实验反馈闭环

这是当前阶段，优先级最高。

目标：让每个 hLF/OPN KO/OE 实验能够与原始预测一一关联，并形成可审计的校准数据，而不是只把实验结果写进自由文本。

### Phase 1 预先约定

以下决策在进入开发前固定，执行会话不得自行改成另一套数据模型。

#### 所有权与模块边界

- 核心模块放在 `python_pichia/src/pcsec_pichia/experimental_feedback/`。
- 建议拆分为 `schema.py`、`io.py`、`quality.py`、`linkage.py`、`calibration.py`。
- `app/services/pichia_experiment_feedback_service.py` 只做 facade、路径和错误汇总。
- Streamlit 页面只调用 service，不实现单位换算、匹配、质量判断或校准公式。
- 原始实验记录不进入 LLM；只有程序生成的脱敏 calibration fact pack 才可在后续显式调用 LLM。
- 核心 schema 默认使用 frozen dataclass、显式 `validate()` 和结构化序列化，保持与现有 `python_pichia` 模式一致；Round 0 不新增 schema 依赖。

首批公开 API 固定为：

```python
load_experiment_bundle(path) -> ExperimentBundle
validate_experiment_bundle(bundle) -> ExperimentValidationResult
write_experiment_feedback_cache(bundle, output_dir) -> ExperimentFeedbackOutputs
build_prediction_index(screen_runs) -> PredictionIndex
link_experiments_to_predictions(bundle, prediction_index) -> PredictionLinkageResult
build_calibration_summary(validated_bundle, linkage_result, config) -> CalibrationSummary
```

具体实现可以增加私有 helper，但 service/UI 不得绕过这些入口重新解析实验表格或计算校准指标。

#### Canonical entity model

内部 canonical schema 使用四类结构化实体，不采用一个不断加列的超宽表作为长期契约：

1. `ExperimentRecord`：一次独立生物培养/一个 biological replicate，包含 target、host、batch、condition、timepoint context 和质量状态。
2. `InterventionRecord`：关联到 experiment 的一个改造组件；支持 control、KO、OE，未来组合改造通过多个 component 表示。
3. `MeasurementRecord`：一次具体 assay 结果；technical replicate 是独立 measurement，不覆盖原值。
4. `PredictionLinkRecord`：experiment/intervention 与 `prediction_run_id`、`evidence_id`、gene/reaction 的匹配及匹配状态。

一个 experiment 可以有多个 intervention component 和多个 measurement。不能把组合改造压成无法解析的自由文本。

#### ID 与版本

- 所有记录包含 `schema_version=1`。
- `experiment_id` 由实验团队或导入模板提供并保持稳定；不得根据会变化的测量值重新计算。
- `intervention_id`、`measurement_id` 在同一 experiment 下唯一。
- import manifest 保留 source file、hash、imported_at、schema version 和 warning。
- 重复 ID、同 ID 不同内容和跨 target 复用必须显式失败或进入 conflict，不静默覆盖。

#### Target、host 与 intervention

- 内置 target 使用 canonical `hLF`、`OPN`；custom target 保留稳定 target_id 和原始名称。
- host 至少记录 species、strain 和 parent strain；空 strain 不能进入正式校准。
- intervention type 固定为 `control`、`KO`、`OE`；组合实验包含多个 intervention component。
- KO 至少记录 gene_id 和构建方法。
- OE 至少记录 gene_id、construct_id、promoter/induction mode；copy number 未知时为空并标记，不猜测。
- common name 只能作为注释，不能代替内部 gene_id 或 PredictionLink。

#### Condition 与对照匹配

- condition 至少包含 medium、carbon source、culture mode、temperature、pH、oxygen/agitation 和 sampling time。
- 未记录的条件使用 missing/unknown，不使用默认实验条件静默补齐。
- fold change、secretion improvement 和 growth retention 都是派生值，不作为唯一原始输入。
- 派生比较必须找到同 target、host、batch、condition、timepoint、assay method 和 unit 的 control/parent；无法匹配时标记 `control_match_missing`。
- 不允许用不同批次、不同 assay 或不同时间点的 control 自动计算 fold change。

#### Measurement 与单位

- 原始值、原始单位和 assay method 永久保留。
- canonical unit 首批支持：titer `mg/L`、biomass `gDCW/L`、specific productivity `mg/gDCW/h`、growth rate `1/h`、time `h`、viability `%`。
- OD600 作为独立 measurement 保存；没有实验室换算系数时不得自动转换为 `gDCW/L`。
- intracellular 与 extracellular measurement 必须记录 compartment，不能混合计算 secretion improvement。
- below LOD、below LOQ、above range、missing、assay failed 和 excluded 必须保留状态与原因，不能转换为 0。
- 单位转换只允许通过显式 conversion registry，并保留原始值；不在 parser 中散落硬编码换算。

#### 重复、质量与排除

- biological replicate 使用不同 `experiment_id` 或明确的 biological replicate ID。
- technical replicate 使用不同 `measurement_id`，由统计层聚合。
- 原始重复值始终保留；均值、标准差和置信区间是派生产物。
- exclusion 不能删除记录，必须保留 `excluded=true`、reason 和 reviewer。
- operator 只保存脱敏 ID；不在 tracked fixture、报告或 LLM fact pack 中保存真实姓名。

#### Prediction linkage 与校准资格

- 原始导入允许暂时没有 prediction link，但这类记录不能进入校准统计。
- 自动链接必须同时匹配 target、gene_id、intervention type 和明确的 prediction run；只匹配 common name 时进入 ambiguous。
- 可进入校准的记录必须：schema 有效、target/host/condition 完整、measurement 有效、control 可匹配、PredictionLink 唯一且 intervention/context 一致。
- 失败实验、阴性结果和 assay failure 都保留；只有 assay failure 不参与效果方向统计。
- calibration 输出作为并行字段和报告，不自动修改 `recommendation_tier`、phenotype fixture 或模型约束。
- 任何 curated promotion 都需要人工 review 和独立 checkpoint。

#### 首批统计边界

- 首批只做描述性统计、方向一致性、top-K hit rate、相对基线富集、rank correlation 和 evidence-tier 分层命中率。
- 不在样本量不足时训练复杂机器学习模型。
- 不把 hLF 数据直接校准 OPN，也不把一个培养条件外推到另一条件。
- 阈值、top-K 和方向判定由配置记录在 calibration manifest，不散落在 UI 或脚本中。

#### 输入输出格式

- canonical cache 使用结构化 JSONL + manifest；CSV/XLSX 只作为导入适配层，不是内部 truth format。
- 首批必须提供脱敏的 hLF、OPN、control、KO、OE、重复、坏单位、缺 control 和 ambiguous link fixtures。
- 默认输出目录：`local_runs/experiment_feedback/<run_id>/`，包含 validated records、conflicts、linkage report、calibration summary 和 manifest。
- 不新增数据库和后台服务；先使用可审计文件缓存，等数据规模证明有必要后再评估持久化数据库。

#### 用户可见完成结果

Phase 1 完成时，研发同事应能在 Streamlit 中：

- 导入标准模板或选择已有实验记录。
- 查看 schema、单位、重复、control 和 prediction linkage 问题。
- 修正或导出冲突，不丢失原始数据。
- 查看 hLF/OPN 分开的 prediction-vs-experiment 摘要。
- 查看哪些候选命中、失败、不可评价，以及下一轮排序依据。

页面不得把“已导入实验”自动展示成“模型已校准”；只有通过校准资格检查的记录才进入统计。

#### 真实数据接入前仍需研发组确认的信息

以下信息不阻塞 Round 0-4 使用脱敏 fixture 开发，但在 Round 5 接入真实数据前必须确认：

- 研发组当前使用的实验记录模板和仪器导出样例。
- hLF/OPN 的实际 assay method、LOD/LOQ 和常用单位。
- biological replicate、fermentation batch 和 technical replicate 的实验室定义。
- 每类实验使用的 control/parent strain 和对照匹配规则。
- 是否存在经实验室确认的 OD600 -> gDCW/L 换算系数。
- 哪些字段涉及人员、客户或项目敏感信息，以及允许的脱敏方式。
- 哪批历史 hLF/OPN 数据获准用于首次回放，哪些只能本地查看。

### Round 0：实验数据契约和验收冻结

交付内容：

- 将上述预先约定转换为明确的 dataclass/Pydantic 契约草案和字段字典。
- 定义 module API、错误类型、status 枚举、import bundle 和 manifest contract。
- 增加脱敏 fixture 与 contract tests，测试先失败以证明缺口存在。
- 不接 Streamlit、不训练模型、不读取真实实验数据。

验收：schema 能区分 KO/OE、hLF/OPN、组合 component、条件、重复、compartment、原始/canonical unit 和 PredictionLink；不允许只有“产量提高/降低”的自由文本。本轮必须产生可执行测试和 fixture，不能只修改文档。

### Round 1：实验记录与缓存 IO

- 在 `python_pichia` 定义 experiment record、measurement、quality flag 和 manifest。
- 支持 CSV/JSONL 导入、稳定 ID、重复记录检查和单位校验。
- 原始导入进入 `local_runs/experiment_feedback/inbox/`。
- 不在本轮训练模型，不改 recommendation tier。

验收：fixture round-trip、坏单位、重复 ID、缺条件和 target 混淆均有测试。

### Round 2：预测与实验关联

- 将 experiment record 关联到原始 screen run、evidence_id、gene/intervention 和目标蛋白。
- 输出 matched、ambiguous、missing prediction、context mismatch。
- 禁止仅凭 gene common name 自动关联。

验收：每条进入校准的数据必须能回溯到预测行和实验来源。

### Round 3：校准指标和排序反馈

- 计算方向一致性、top-K hit rate、相对基线富集、rank correlation 和按 evidence tier 分层的命中率。
- 小数据阶段只生成统计和校准报告，不训练复杂黑盒模型。
- 不把一次实验自动提升为跨 target/context 的普遍规律。

验收：指标可由 fixture 手算复核；缺失值和失败实验不会被静默丢弃。

### Round 4：service 与 Streamlit 回填入口

- `app/services` 只调用核心实验反馈 API。
- Streamlit 支持选择预测、录入或导入实验、查看匹配状态和校准摘要。
- 页面不直接修改模型、tier 或稳定科学资产。

验收：页面入口、导航、旧引用、session state、cache key 和错误提示完整检查。

### Round 5：hLF/OPN 回放验收

- 使用脱敏 fixture 或已批准实验数据完成 hLF、OPN 各一组回放。
- 生成 prediction-vs-experiment 报告和下一轮候选排序建议。
- 只有人工确认后，才可把数据提升到稳定 curated experiment cache。

结束条件：项目能够回答“模型推荐了什么、实验观察到什么、排序是否变得更可信”。

## Phase 2：gene-level OE 与酶容量

进入条件：Phase 1 数据契约稳定，至少能评价 OE proxy 与实验方向的一致性。

目标：把 reaction-level proxy 升级为 gene-enzyme-reaction capacity，表达 OE 剂量和蛋白资源成本。

主要内容：

- 评估并人工复核 iPichia/ecPichia 等外部模型资产。
- 建立 gene-enzyme-reaction mapping，不原样覆盖当前 GPR。
- 表达 promoter/copy-number 档位、酶容量、复合体亚基和总蛋白成本。
- 对缺 kcat/丰度数据使用区间和 uncertainty，不硬编码成单一真值。
- 保留旧 OE proxy 作为兼容和对照路径。

验收：同一 OE 在 proxy 和 gene-level capacity 路径中的差异可解释；没有 locus/GPR 的候选不得伪装成 gene-level 可执行。

## Phase 3：分泌资源与蛋白稳态约束

目标：覆盖传统 GEM 难以表达的分泌瓶颈。

优先机制：

- ER 转运、折叠、二硫键和伴侣蛋白。
- UPR、ERAD、蛋白酶体和错误折叠。
- 糖基化、囊泡运输和胞吐。
- 目标蛋白特异的翻译、修饰和降解成本。

实现要求：先定义资源池和约束语义，再映射基因；不能把文献基因名单直接变成模型反应。

验收：每个新增约束都有单位、来源、开关、基线回归和 hLF/OPN 差异解释。

## Phase 4：组合改造与条件鲁棒性

目标：从单基因列表进入可执行的组合设计。

- 只从经过 Phase 1-3 过滤的 top candidates 搜索 KO+KO、KO+OE、OE+OE。
- 使用成熟 MILP/搜索库或清晰的可替换搜索后端，不写候选特例硬编码。
- 在碳源、生长率、氧供、target 和参数扰动下评价方向稳定性。
- 显式记录上位性、不可行组合和生长代价。

验收：组合收益不能由单基因分数简单相加；结果必须包含跨条件稳定性和失败原因。

## Phase 5：前瞻验证

目标：用未参与校准的新实验判断系统是否真正提高研发效率。

- 在实验前冻结候选排名和证据快照。
- 同时设置模型 top-K、研发经验基线和必要的阴性/风险对照。
- 评价 top-K hit rate、富集倍数、排序校准和重复实验一致性。
- 失败结果必须回流，不只保存成功案例。

项目达到这一阶段后，才能声称“帮助提高 KO/OE 实验命中率”；仍不能声称保证绝对产量。

## 暂不做

- 不把 COBRApy/optlang 切成 full pcSec 默认后端。
- 不把外部 GPR、BLAST/RBH 或数据库注释自动写入当前模型。
- 不在缺少实验数据时训练复杂机器学习模型。
- 不把 reaction proxy 改名为 gene-level OE。
- 不让 LLM 直接读取任意目录、修改 tier 或生成未经 validator/Judge 的最终结论。
- 不迁移或覆盖 legacy `Results/`。

## 当前下一步

从 Phase 1 Round 0 开始：冻结实验数据契约、输入样例、质量规则和验收标准。完成后按 Round 1-5 顺序推进，不返回已完成的 COBRApy、BLAST 或在线证据阶段重新开发。
