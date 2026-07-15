# pcSecPichia 下一阶段执行计划

状态：active  
最后更新：2026-07-15

## 当前执行位置

```yaml
current_program: mvp_directions_1_to_3
current_slice: direction_2_oe_product_tiering_closure
slice_status: ready
```

执行会话必须从这里和 `docs/handoff.md` 记录的切片继续。历史 Phase/Round 只用于理解既有实现，不得覆盖项目级执行计划的当前授权顺序。

方向 1 的软件闭环和真实发酵模板适配已经通过验收。获批真实数据到来后执行独立回填 checkpoint；当前按项目级顺序进入方向 2 产品层级与门禁收口。

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

1. 读取 `docs/README.md`、`docs/EXECUTION_PLAN.md`、当前架构、当前计划和数据治理策略；先确认本轮仍在项目级授权范围内。
2. 确认本轮输入、输出、数据边界和验收标准。
3. 先实现核心 `python_pichia` 数据契约或科学逻辑。
4. 增加 focused tests 和 fixture。
5. 再接 `app/services` facade；只有用户工作流需要时才接 Streamlit。
6. 运行 smoke，检查产物进入 ignored `local_runs/`。
7. review、自动修复、复验，最多 3 轮。
8. 更新当前状态；不得重写已完成阶段的大段设计。
9. 检查 `Code/Model/Enzymedata/Results` 和依赖声明没有非预期 diff。

## 子 agent 协作策略

当前阶段可以使用内置子 agent 提高速度，但所有子任务仍属于同一个目标和同一个主分支。不要再创建多个用户会话并手动切换共享工作区分支。

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
- smoke 输出使用独立 tmp_path 或当前阶段约定的唯一 `local_runs/<phase>/<round>/<agent>/`，避免相互覆盖。

### 适合并行的任务

- 多个互不依赖的现有代码路径审计。
- 核心接口冻结后的独立实现与独立测试/fixture。
- service 和 UI 在 facade contract 已冻结后的并行开发。
- hLF、OPN 独立 smoke 和只读安全审计。

### 不适合并行的任务

- 当前阶段 canonical schema、public API 和数学语义的最终决定。
- 会影响后续模块的 ID、单位、执行状态、参数优先级和约束含义。
- 同一个文件或同一个测试 fixture 的修改。
- stage、commit、push、迁移稳定科学资产或修改保护目录。

当前 Round 的具体分工以当前 Phase 章节为准。每轮最多进行一次主并行波次；若子任务暴露接口问题，先由主 agent 收束接口，再决定是否进行第二波，不能无限创建 agent 掩盖架构不清。

## Phase 1：实验反馈闭环

本阶段已完成并通过端到端验收；后续只在真实数据回放 checkpoint 中补充获批实验数据，不再重开实现轮次。

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
- CSV/XLSX 首批共享 `record_type + payload_json` envelope；XLSX 优先读取 `records` 工作表，不存在时读取活动工作表。
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

Phase 1 和方向 1 模板适配已完成并验收。Phase 2 的 mapping、剂量、约束、求解、报告、UI 和外部候选审计已存在，但正式容量资产为空。当前不再扩大同类来源搜索，而是按 ADR-002 收口 reaction proxy、relative uncalibrated、absolute unavailable 和 not executable 的产品状态与门禁。

### 当前执行锁

- 当前只执行 `direction_2_oe_product_tiering_closure`，必须修改实际核心状态/门禁及其测试。
- 只读审计、继续搜索同类绝对容量来源、重复既有 smoke 或只修改文档均不构成交付。
- 不删除或重写既有 `oe_capacity` 能力；在当前实现上收口产品语义并保持回归。
- 验收前不得进入 Phase 3、组合搜索或完整跨条件排名。

### 目标与边界

目标是把 reaction-level proxy 升级为可审计的 gene-enzyme-reaction capacity，表达基因到酶和反应的映射、OE 剂量、复合体/同工酶语义、参数不确定性和蛋白资源成本。

- 保留旧 reaction proxy 作为独立兼容和对照模式。
- 复用当前 GPR parser、enzyme data、protein pool、SciPy HiGHS、Shadow LP 和 reference solve。
- 不新增默认 solver，不重写 Shadow LP，不修改稳定模型资产。
- 外部 GEM/GPR、BLAST/RBH、数据库和同源参数只作为候选证据，不自动覆盖当前模型。
- Phase 2 不实现 UPR、ERAD、糖基化或囊泡资源池，这些属于 Phase 3。
- 输出仍是模型内相对比较，不预测 mg/L、真实表达倍数或实验成功率。

### 模块布局

```text
python_pichia/src/pcsec_pichia/oe_capacity/
  __init__.py
  schema.py
  mapping.py
  parameters.py
  constraints.py
  simulation.py
  reports.py
```

- `schema.py`：frozen dataclass、枚举、status 和 validation。
- `mapping.py`：复用当前模型 GPR、gene index 和 enzyme data，生成 mapping catalog。
- `parameters.py`：dose、kcat、MW、经审核的绝对 baseline capacity、复合体化学计量和 uncertainty scenario。
- `constraints.py`：生成 backend-neutral 结构化容量修改，不调用 solver。
- `simulation.py`：调用现有模型准备和求解路径，执行 baseline/proxy/gene-capacity 对照。
- `reports.py`：输出 rows、coverage、manifest、parameter trace 和差异报告。
- `screens/gene_interventions.py` 继续拥有 GPR 解释和旧 proxy；不得复制第二套 parser。
- `external_refs` 继续拥有联网抓取与外部 cache；正式 screen 求解不联网。

### 核心契约

`ParameterEstimate`：参数名、nominal/lower/upper、unit、source type/ref/version、confidence、是否同源转移和 warnings。区间必须满足 `lower <= nominal <= upper`。

`OEDoseSpec`：dose id、`explicit_multiplier/promoter_copy_mapping/categorical_only`、multiplier、promoter、copy number、induction、mapping source、uncertainty 和 warnings。没有审核 dose mapping 的类别输入不得自动变成 multiplier。

`GeneEnzymeReactionMapping`：model fingerprint、gene/enzyme/reaction id、GPR、GPR role、complex/subunit、enzyme variable、formation/dilution reaction、source、confidence、execution status 和 warnings。

`GeneCapacitySpec`：mapping、kcat、MW、baseline enzyme amount、complex stoichiometry、dose、parameter scenario 和 resource cost mode。

`OECapacityPlan`：gene/target/context、dose、execution mode/status、executable specs、explain-only mappings、proxy reactions、constraint changes、uncertainty scenarios、missing information 和 warnings。

`OECapacityComparisonResult`：baseline/proxy/gene-capacity solver status、secretion objective、growth retention、最大可行生长率、protein resource cost、scenario 结果、差异、traceability 和 skipped reason。

execution mode 固定为 `gene_capacity`、`reaction_proxy`、`comparison`、`not_executable`；不得静默降级或改名。

execution status 至少包括：

- `gene_level_executable`
- `partial_mapping`
- `isoenzyme_ambiguous`
- `complex_limited`
- `external_evidence_only`
- `categorical_dose_only`
- `proxy_only`
- `unresolved`

### 语义规则

- gene、reaction、enzyme entity 都映射到当前模型后才能成为 gene-level executable。
- 同工酶 OE 只能改变明确属于该基因的酶容量项，不能放宽整个反应上界。
- 单个复合体亚基 OE 默认不能提高完整复合体容量；缺限制亚基/组装证据时为 `complex_limited`。
- mixed GPR 必须保留每条 mapping 的作用范围。
- 外部 GPR、同源关系和外部模型参数不能单独提升 execution status。
- 优先使用当前 pcSecPichia GPR 和本地 enzyme data；外部 Pichia 模型、数据库、文献、同源转移和 smoke fixture 依次降级。
- 缺参数使用显式 low/nominal/high scenario，不静默填单值。
- 如果模型存在 enzyme amount 或 formation/dilution variable，gene-level OE 修改对应容量项并计入资源成本。
- baseline capacity 必须来自带模型、上下文、来源版本和审核信息的容量资产；不得使用最优 formation flux、通用 reaction 上界或固定相对值补齐。
- 每个 uncertainty scenario 必须保存使用相同参数的 capacity baseline 与 perturbation；`1.0x` 是该成对模型的恒等变换。
- 仅把代谢反应 bound 乘以 factor 的路径仍是 reaction proxy。
- `1.0x` gene-capacity 必须与 baseline 在容差内一致；feature-off 和 proxy-only 必须保持旧结果回归。
- 放宽容量不保证 secretion 一定提高，必须保留无变化、降低、不可行和资源成本上升结果。

### 公共 API

```python
build_gene_enzyme_reaction_catalog(model, metabolic, combined, external_evidence=None) -> GeneCapacityCatalog
validate_gene_capacity_catalog(catalog) -> GeneCapacityValidationResult
build_oe_dose_spec(payload, dose_mapping=None) -> OEDoseSpec
build_gene_capacity_specs(gene_id, catalog, dose, parameter_policy) -> tuple[GeneCapacitySpec, ...]
plan_gene_level_overexpression(model, gene_id, target_id, context_id, dose, catalog, parameter_policy) -> OECapacityPlan
build_oe_capacity_constraints(prepared_model, plan) -> OECapacityConstraintBundle
run_gene_level_oe_comparison(prepared_model, plan, solver_options=None) -> OECapacityComparisonResult
run_gene_level_oe_screen(prepared_model, requests, screen_config) -> OECapacityScreenResult
write_oe_capacity_outputs(result, output_dir) -> OECapacityOutputs
```

service/UI 不得绕过这些 API 自行解释 GPR、换算剂量或修改约束。

### Round 0-6

1. **Round 0 架构与可执行契约**：审计 GPR/enzyme/protein/solver 链路；创建 `oe_capacity` 包、schema、validation、fixtures、contract tests。禁止 service/UI 和真实 screen。
2. **Round 1 mapping catalog**：从当前模型建立 gene-enzyme-reaction catalog 和 coverage；外部 evidence 只增加 traceability。
3. **Round 2 dose/parameter/uncertainty**：实现 dose、参数优先级、冲突、low/nominal/high、isoenzyme/complex/missing-parameter 状态。
4. **Round 3 constraint 与单候选求解**：构建真实容量约束；运行 baseline、1x、proxy、gene-capacity scenarios；验证资源成本和回归。
5. **Round 4 screen 与对照报告**：实现小批量/catalog screen；并列输出 proxy/gene-capacity；写入 `local_runs/oe_capacity/`。
6. **Round 5 service/Streamlit**：增加薄 facade 和页面工作流；展示 mapping、dose、参数、uncertainty、资源成本和不可执行原因。
7. **Round 6 hLF/OPN 验收**：分别运行 executable 与边界候选，生成 coverage、数值回归、性能、保护目录和数据泄漏报告。

### Round 6A：外部 baseline capacity 候选与审核提升

本轮是 Round 6 验收阻断后的数据获取补充轮，不重新实现 Round 0-5。

固定执行顺序：

1. 扩展容量候选契约，区分 `target_specific`、`host_condition`、`external_model_calibrated` 和 `homolog_transferred`，并定义条件匹配与来源置信度。
2. 复用现有 `external_refs` 联网与 cache 边界，接入 Pichia 定量蛋白组、iPichia/ecPichia 和动力学来源；无法稳定自动获取的来源生成明确的人工导入入口，不硬编码网页结果。
3. 建立 gene/enzyme/formation 当前模型 crosswalk、单位转换、`abundance × kcat -> model_flux capacity` 换算、low/nominal/high 区间和冲突检测。
4. 生成候选清单和 promotion manifest；默认只写 ignored `local_runs/oe_capacity/`，不得自动修改正式资产。
5. 在 Streamlit 增加来源、映射、条件、换算和冲突复核界面；只有用户明确批准后，才通过核心 promotion API 原子更新正式容量资产。
6. 优先为 `PAS_chr2-1_0308 / G6PDH2` 寻找可审核候选，分别评估 hLF 与 OPN 当前 glucose、`mu=0.1` context；宿主条件参数应复用同一 provenance，不复制成虚假的 target-specific 实验。
7. 重跑 focused tests、联网/离线 cache smoke、无网络回放、Streamlit 浏览器验证和安全检查。

本轮验收标准：

- 至少一个真实外部来源完成 fetch 或可重复人工导入、缓存、版本/hash/license 和解析验证。
- 外部 gene/protein ID 必须映射到当前模型 gene/enzyme/formation；BLAST/RBH 只能辅助映射，不能单独生成容量。
- 单位转换链和原始值可追溯，错误单位、缺失 biomass conversion、条件不匹配和冲突不会静默进入正式资产。
- `host_condition` 与 `target_specific` 的匹配行为有测试，hLF/OPN 不会因复制记录获得虚假 target-specific 证据。
- 同源转移保持低置信区间，不能单独把正式 Phase 2 门禁提升为通过。
- promotion 必须由显式批准触发，写入前后均校验 schema、模型指纹和 source hash，并保留旧资产备份或原子替换证据。
- 没有合格候选时仍保持 `reviewed_baseline_capacity`，不得回退到通用上界、baseline flux 或 fixture。

A0b 历史结果：已通过正式 adapter/cache/parser 接入 PRIDE `PXD055501`，为 G6PDH2 保留 T0 iBAQ 原始值、CC0 许可、项目版本、原始文件 hash、培养条件、`504/504` sequence crosswalk 和不可转换的单位链；同时增加 ecPichia `Supplementary 8.yml` 的正式文件导入和离线回放，保留 gene/enzyme/reaction、MW、kcat、GECKO coefficient、reported concentration、protein pool、artifact hash 和上游 ZIP hash。PRIDE 只有相对 iBAQ；ecPichia 仍缺可审核的 LFQ→absolute abundance provenance、确定单位、完整条件、可复用许可、直接 kcat provenance 和当前 formation-flux 换算。因此两者均 `promotion_ready=false`。该结果已由 ADR-002 收口，不再驱动同类来源扩张。

#### A0c：ecPichia provenance closure

A0c 是 Round 6A 的限界闭合轮，不是第三次宽泛 source acquisition。主候选固定为 ecPichia G6PDH2；PRIDE 仅作为独立表达证据和条件对照。本轮不接入只重复提供相对强度、名称或模型字段的新来源。

固定执行顺序：

1. 审计 ecPichia 论文方法、补充材料、正式仓库和 GECKO 版本语义，确认 concentration、protein pool、coefficient 的单位、biomass basis 和生成方式。
2. 用 MW、kcat 和时间单位独立复算 GECKO coefficient；任何不一致必须进入 conflict，不能选择性忽略。
3. 追踪 G6PDH2 kcat 的直接来源，至少保留 source record、enzyme/reaction/substrate、organism、条件、版本和转移关系；只有 `brenda` 标签不构成完整 provenance。
4. 建立外部 enzyme/reaction 到当前 `PAS_chr2-1_0308 / C4R099` 及 formation/dilution handle 的可审计 crosswalk。
5. 形成逐步单位换算 trace，将 concentration/abundance 与 kcat 转为当前模型 `model_flux` capacity；无法闭合任一维度时不得产生 nominal capacity。
6. 对宿主、菌株、培养基、碳源和 `mu` 差异生成 applicability、low/nominal/high uncertainty 和 warning，不复制为虚假的 target-specific 参数。
7. 运行联网 fetch、冻结 cache、离线 replay、candidate evaluation 和 promotion preview；promotion 仍需明确审核批准。
8. 运行 review/fix/verify，最终写出 evidence-chain report 或 gap report，并进入一个且仅一个完成出口。

A0c 验收：

- 关键外部陈述均能回溯到正式文件、论文、数据库记录或规范，并保存 URL/version/hash/license/retrieved_at。
- coefficient 可由公开字段和明确公式复算，单位链测试覆盖正确值、错误单位、缺字段和条件冲突。
- current-model crosswalk 精确到 gene、enzyme、reaction、formation/dilution handle 和模型指纹。
- 在线路径与无网络 replay 生成相同的标准化 source assessment。
- 没有完整链时不会生成 nominal `model_flux` capacity，不会修改正式 registry。
- 本轮不得以继续搜索结束，必须进入以下一个出口：
  - `evidence_chain_closed`：生成 promotion preview；明确批准并成功提升后进入 Round 6B `ready`。
  - `architecture_decision_required`：保存 gap report，停止增加同类 source adapter，下一轮先决定是否拆分绝对容量与相对 OE 场景的验收等级。

A0c 历史结果：已完成 ecPichia `Supplementary 8.yml` 与 `Supplementary 11_V2.docx` 的 hash 固定、正式解析、在线导入和离线 replay，并建立 hLF/OPN 的 `glucose_mu_0.1` current-model crosswalk。GECKO coefficient 可严格复算，但表格与 YAML 的 gene/enzyme/MW/concentration 冲突；concentration 的 `g/L`、未声明 YAML 单位和缺失 biomass conversion 无法闭合；本地研究记录将 `kcat=8000 s^-1` 指向 Thermotoga/thio-NADP+/80°C，但该记录尚未绑定可 replay 的 source artifact，且不能作为 Komagataella direct kcat；catalytic flux 到 formation/dilution `model_flux` 也缺直接转换证据。因此出口为 `architecture_decision_required`，不生成 nominal capacity 或 promotion preview，正式 registry 不变。ADR-002 随后将绝对容量研究层与相对决策层分离。

### Round 6B：hLF/OPN 绝对容量重验收（条件路径）

只有未来获得并批准合格绝对容量资产时，才重新启用本路径。正式 runner 新建四个固定 smoke：hLF/OPN 各一个 `PAS_chr2-1_0308` executable case 和一个 `PAS_chr1-4_0458` boundary case，并验证 target、context、gene、资产 hash、scenario/proxy 证据和相关回归。该条件路径是否通过不再阻塞相对 OE 产品收口或方向 3 Round 0。

绝对容量研究层的条件状态迁移为：

- Round 6A 功能完成但没有已批准候选：`round_6a_external_capacity_candidates / awaiting_candidate_review`。
- A0c 证据链无法闭合：`round_6a_external_capacity_candidates / architecture_decision_required`，不得继续堆叠同类来源。
- 正式容量资产已批准：`round_6b_hlf_opn_acceptance / ready`，立即运行正式验收。
- Round 6B 未通过：`round_6b_hlf_opn_acceptance / awaiting_acceptance`，保留真实缺口，不回退伪参数。
- Round 6B `passed=true`：仅表示绝对容量研究层通过，不替代当前项目级方向状态。

每个 Round 都必须有生产实现、focused tests 和可验证结果。不得跳到后续 Round，也不得反复选择最简单字段。

子 agent 分工建议：Round 0 并行现有链路审计与独立 contract fixtures；Round 1 并行 catalog/coverage tests；Round 2 并行 parameter 实现与坏输入测试；Round 3 并行 constraint builder 与 tiny-model 回归；Round 4 并行 screen 与 report/fact-pack tests；Round 5 在 service API 冻结后并行 service/UI；Round 6 并行 hLF、OPN smoke 和只读安全审计。schema、数学语义、参数优先级、整合、commit 和阶段判断始终由主 agent负责。

### 输出与验收

screen rows 并行保留旧字段，并新增 execution mode/status、dose、expression multiplier、mapping、parameter source/confidence、uncertainty scenarios、gene-capacity objective、proxy objective、差异、resource cost、missing information 和 warnings。

当前方向 2 产品收口必须满足：

- single gene、isoenzyme、complex、mixed、missing parameter 和 external-only 都有测试。
- reaction proxy、relative uncalibrated、absolute unavailable 和 not executable 可区分，且由核心层统一判定。
- relative gene-capacity 是结构化相对场景，不是 proxy 改名，也不声称绝对 capacity。
- 旧 proxy 可逐候选比较且兼容模式数值不变；feature-off 回归通过。
- 缺审核 baseline capacity 时不调用绝对容量求解、不生成 nominal capacity，并明确返回 unavailable。
- hLF/OPN 分别完成产品状态 smoke，保留 target、context、model version 和风险说明。
- service/UI 只透传和展示核心判断。
- 输出进入 ignored `local_runs/oe_capacity/`。
- `Code/Model/Enzymedata/Results` 和依赖声明没有非预期 diff。
- 状态推进到 `direction_3_secretory_resource_round_0 / ready` 后停止。

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

执行 `direction_2_oe_product_tiering_closure`。具体成功条件和预算边界读取 `docs/EXECUTION_PLAN.md`，产品分层读取架构文档与 ADR-002，当前切片和验证方式读取 `docs/handoff.md`；验收完成前不得进入方向 3。
