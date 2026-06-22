# pcSecPichia Python 引擎详细重构计划

日期：2026-06-22

关联文档：[pcSecPichia Python 移植与研发用途决策记录](pichia_python_migration_strategy.md)

## 1. 重构目标

本轮重构不是把整个 `pcSecYeastSpecies` 翻译成 Python，而是把研发真正需要的毕赤酵母目标蛋白分泌优化能力逐步 Python 化。

最终目标：

> 在不依赖 MATLAB 的情况下，使用 Python 对毕赤酵母 pcSecPichia 模型进行目标蛋白分泌仿真，并筛选可能提高 OPN、hLF 等目标蛋白产量的敲除/过表达候选。

第一阶段目标更小：

> 建立可并行存在的 Python pcSecPichia engine，使现有 Streamlit / service 层可以在 MATLAB engine 和 Python engine 之间切换。

## 2. 非目标

首轮不做：

- 三物种全量迁移；
- 全部 MATLAB 文件逐行翻译；
- 全部论文图复现；
- 完整 `buildModel_pcSecPichia.m` 重写；
- 温度依赖 ETC 全量迁移；
- humanized glycosylation 全量迁移；
- 真实发酵产量承诺。

这些可以作为后续版本增强，而不是当前重构验收条件。

## 3. 架构设计原则

本次重构必须遵循 **高内聚、低耦合、可扩展** 的设计原则。

### 3.1 高内聚

每个模块只负责一个清晰领域：

- target input 层只负责目标蛋白、leader、PTM 参数的输入和校验；
- model layer 只负责 pcSecPichia 模型结构、reaction/gene/enzyme 数据；
- engine 层只负责计算；
- perturbation 层只负责 KO/OE 扰动定义和筛选；
- visualization 层只负责把模型关系组织成网络图数据；
- service 层只负责任务编排；
- UI 层只负责展示和交互。

不允许 UI 直接操作 `.mat`、LP 矩阵、MATLAB 命令或 solver 参数。

### 3.2 低耦合

核心边界应通过协议接口和稳定数据结构连接：

```text
UI -> Service -> Engine Protocol -> MATLAB Engine / Python Engine
UI -> Service -> Visualization Builder -> Plotly/Cytoscape Adapter
Target Input Provider -> Target Schema -> Engine
```

要求：

- Streamlit 页面不依赖具体 engine；
- Python engine 不依赖 Streamlit；
- MATLAB engine 和 Python engine 返回同一种 result schema；
- SigScout、PichiaCLM、后续实验数据库只能通过 provider/export adapter 接入；
- 可视化图不能直接耦合 solver 内部结构，应接收标准化 graph model。

### 3.3 可扩展

系统必须支持后续加入其它蛋白质。

OPN 和 hLF 只是内置示例，不应写死到 engine 中。新增蛋白时应只需要：

1. 增加 target protein 配置；
2. 增加 leader/signal peptide 候选；
3. 补充 PTM 参数和证据来源；
4. 通过同一套 target schema 进入模型。

如果目标蛋白 PTM 信息不完整，系统应进入“草案模式”，不能把结果标记为正式推荐。

## 4. 原项目保留策略

重构不是删除 MATLAB 项目，也不是一次性替换原论文工作流。

原项目应完整保留：

- `Code/pcSecPichia/`
- `Code/pcSecYeast/`
- `Code/pcSecKmarx/`
- `Code/Figures/`
- `Model/`
- `Enzymedata/`
- `Results/`

新的 Python engine 作为并行运行层加入：

```text
原 MATLAB 代码：参考实现、论文复现、baseline 对齐
新 Python engine：部署、网页调用、批量筛选、研发决策支持
```

在 Python engine 完全通过对齐测试之前，MATLAB engine 必须保留为回退路径。

原则：

- 不删除原始 MATLAB 构建和仿真脚本；
- 不移动原始模型和结果数据；
- Python 输出写入新的 `local_runs/` 或专门结果目录；
- 每个 Python 迁移阶段都要能和 MATLAB baseline 对比。

## 5. 当前系统问题

### 5.1 MATLAB 依赖过重

当前仿真入口通过 Python 调用 MATLAB：

- `app/adapters/matlab.py`
- `app/services/simulation.py`
- `app/services/opn.py`
- `local_opn_pichia_glc.m`
- `local_pichia_ref_glc.m`

这导致：

- 部署需要 MATLAB 商业授权；
- 生物学专家难以理解运行链路；
- 服务器化和云端部署成本高；
- 后续做批量 KO / OE 筛选时调度复杂。

### 5.2 目标蛋白输入仍偏 OPN

当前已有 `app/core/opn_inputs.py`，但它仍然以 OPN 命名。

后续需要支持：

- OPN；
- hLF；
- 其它目标蛋白；
- 多个 leader / signal peptide 候选；
- PTM 参数；
- 实验角色和证据说明。

因此需要抽象成通用 target protein input，而不是继续扩大 OPN 专用接口。

### 5.3 模型运行和结果解释耦合

当前服务层既要调用 MATLAB，又要找 LP 文件、解析 SoPlex 输出、组织 UI 消息。

重构后需要把这些职责拆开：

- engine 负责模型计算；
- service 负责任务编排；
- UI 负责展示；
- parser/exporter 负责读写结果；
- explanation 层负责把模型结果转成生物学可读说明。

### 5.4 缺少模型网络可视化

当前前端能展示结果表和日志，但还不能把分泌模型以“路径、瓶颈、基因扰动”的形式可视化。

研发用户需要看到：

- 目标蛋白从翻译到胞外分泌经过哪些步骤；
- 哪些模块可能限制产量；
- KO/OE 候选影响的是哪些反应、复合体或分泌步骤；
- 模型推荐和生物学机制之间的关系。

因此后续需要新增分泌模型可视化模块。

## 6. 推荐目录结构

新增或演进为：

```text
app/
  core/
    target_inputs.py
    pichia_model.py
    perturbation.py
    graph_models.py
  engines/
    __init__.py
    base.py
    matlab_pichia_engine.py
    python_pichia_engine.py
  services/
    pichia_target_simulation.py
    perturbation_screening.py
    bottleneck_explanation.py
    secretion_graph.py
  adapters/
    matlab.py
    lp_solver.py
    mat_loader.py
    graph_renderer.py
  ui/
    views/
      pichia_targets.py
      pichia_perturbations.py
      pichia_network.py
```

保留：

```text
Code/pcSecPichia/
Model/pcSecPichia.mat
Enzymedata/pcSecPichia/
Data/pcSecPichia/
```

这些作为 Python 迁移期间的参考实现和数据来源。

## 7. 核心接口设计

### 7.1 Engine 接口

新增 `app/engines/base.py`：

```python
class PcSecPichiaEngine(Protocol):
    def simulate_target(self, request: TargetSimulationRequest) -> TargetSimulationResult:
        ...

    def screen_gene_knockouts(self, request: GeneKnockoutScreenRequest) -> GeneKnockoutScreenResult:
        ...

    def screen_capacity_changes(self, request: CapacityScreenRequest) -> CapacityScreenResult:
        ...
```

短期只要求实现 `simulate_target()`。

中期再实现 KO/OE 相关接口。

### 7.2 目标蛋白输入

新增 `app/core/target_inputs.py`：

```text
TargetProteinInput
LeaderCandidateInput
TargetConstructInput
TargetProteinInputProvider
CsvTargetProteinInputProvider
BuiltinTargetProteinInputProvider
```

字段建议：

```text
target_id
display_name
mature_sequence
leader_sequence
signal_peptide_sequence
through_er
disulfide_sites
n_glycosylation_sites
o_glycosylation_sites
transmembrane
gpi_site
localization
cotranslation
evidence_note
risk_note
```

OPN 迁移后只是一个内置 target preset。

hLF 应作为第二个内置 target preset，但必须补全成熟序列和 PTM 参数后再进入模型。

新增蛋白的接入方式：

```text
BuiltinTargetProteinRegistry
CsvTargetProteinInputProvider
FastaTargetSequenceProvider
ExperimentalTargetMetadataProvider
```

首轮可以只实现 builtin + CSV，但接口要允许未来接入实验数据库或网页上传。

### 7.3 模型数据结构

新增 `app/core/pichia_model.py`：

```text
PcSecPichiaModel
ReactionIndex
GeneRule
EnzymeData
```

需要承载：

```text
rxns
mets
S
lb
ub
c
b
genes
rules
grRules
rxnGeneMat
```

以及酶数据：

```text
enzyme
subunit
subunit_stoichiometry
kcat
proteins
proteinMWs
proteinLength
proteinExtraMW
rxns
rxnscoef
label
```

### 7.4 分泌网络图接口

新增 `app/core/graph_models.py`：

```text
GraphNode
GraphEdge
SecretionNetworkGraph
PerturbationGraphOverlay
```

节点类型：

```text
target_protein
process_module
reaction
enzyme_complex
gene
metabolite
constraint
```

边类型：

```text
passes_through
catalyzed_by
encoded_by
constrains
perturbed_by
increases_flux
decreases_flux
```

可视化服务不直接画图，而是输出标准 graph model。具体渲染可以由 Plotly、PyVis、Cytoscape.js 或其它前端适配器完成。

## 8. 阶段化实施计划

## 阶段 0：基线冻结

目标：确保之后每一步都能和当前 MATLAB 结果对齐。

任务：

- 记录当前可运行的 MATLAB OPN smoke 参数；
- 保存一组代表性 LP 和 `.out`；
- 记录 objective value；
- 记录输入 target CSV；
- 记录模型和 enzymedata 文件路径；
- 给当前 Python app 跑测试。

验收：

- `python -m pytest -q` 通过；
- 有一份固定的 MATLAB baseline 结果；
- baseline 中至少包含 OPN 的一个候选，例如 `OPN_PPA_DDDK18`；
- 后续 Python engine 以此为对齐标准。

建议输出：

```text
local_runs/baselines/pichia_opn_glc_mu0p10/
  input.json
  matlab_generated.lp
  matlab_soplex.out
  summary.json
```

## 阶段 1：抽象 Engine 边界

目标：先重构接口，不改变现有行为。

任务：

- 新增 `app/engines/base.py`；
- 新增 `app/engines/matlab_pichia_engine.py`；
- 把当前 `OpnSimulationService` 中调用 MATLAB 的逻辑包装进 `MatlabPichiaEngine`；
- service 层只依赖 engine 接口；
- UI 不直接知道 MATLAB 或 Python。

验收：

- Streamlit 现有仿真入口行为不变；
- 仍可调用 MATLAB 生成 LP；
- 测试不减少；
- 原有 `local_opn_pichia_glc.m` 不改或只做兼容性小改。

风险：

- 过早抽象会显得繁琐；
- 但这是后续替换 Python engine 的必要边界。

## 阶段 2：通用目标蛋白输入

目标：从 OPN 专用输入升级为 OPN/hLF/其它目标蛋白通用输入。

任务：

- 新增 `app/core/target_inputs.py`；
- 将 `app/core/opn_inputs.py` 中通用部分迁移过去；
- 保留 `opn_inputs.py` 作为兼容 wrapper；
- 新增 hLF 输入 fixture；
- 新增 CSV 导入 schema；
- 新增 target construct 导出能力；
- 调整 CDS/PichiaCLM 相关服务，让其依赖通用 target input。
- 新增 target registry，允许未来添加其它目标蛋白。
- 为每个 target 标记参数状态：
  - `draft`
  - `ready_for_model`
  - `validated`

验收：

- OPN 当前测试全部通过；
- hLF 可以被读取为 target input；
- hLF 如果 PTM 参数未确认，UI/服务必须显示“参数待确认”，不能静默默认；
- 生成的 target protein row 与 pcSecPichia 格式兼容。
- 新增第三个测试蛋白 fixture 不需要改 engine 代码。

关键测试：

- OPN 从旧接口迁到新接口后输出 CSV 完全一致；
- hLF 输入缺失必要 PTM 参数时返回中文错误；
- CSV provider 能接入外部候选 leader。
- target registry 可以列出所有可用目标蛋白。

## 阶段 3：Python 读取 pcSecPichia 模型

目标：不用 MATLAB，读取已有 `.mat` 模型和酶数据。

任务：

- 新增 `app/adapters/mat_loader.py`；
- 读取 `Model/pcSecPichia.mat`；
- 读取 `Enzymedata/pcSecPichia/*.mat`；
- 把 MATLAB struct 转为 Python dataclass；
- 支持稀疏矩阵；
- 建立 reaction/metabolite/gene 索引；
- 实现 `change_rxn_bounds()`；
- 实现 `set_objective()`。

验收：

- 能读取模型字段：
  - `rxns`
  - `mets`
  - `S`
  - `lb`
  - `ub`
  - `c`
  - `b`
  - `genes`
  - `grRules`
  - `rxnGeneMat`
- 能按 reaction id 修改上下界；
- 能保存/导出基础模型摘要；
- 不依赖 MATLAB。

关键测试：

- `BIOMASS`、`Ex_glc_D`、`Ex_o2` 等关键 reaction 能被定位；
- 修改 bounds 后数值符合预期；
- 模型维度与 MATLAB struct 一致。

## 阶段 4：基础 LP 求解

目标：Python 能对读取的模型做最基础 LP 求解。

任务：

- 新增 LP solver adapter；
- 优先支持 SciPy HiGHS 或外部 SoPlex；
- 把 COBRA 风格模型转成 LP；
- 支持最大化/最小化；
- 输出 objective、status、关键通量。

验收：

- Python 能求解基础 pcSecPichia 模型；
- 对固定培养基条件，求解状态为 optimal 或能解释不可行原因；
- 与 MATLAB/SoPlex 结果在容许误差内一致。

建议优先级：

1. 直接用 `scipy.optimize.linprog(method="highs")` 做最小闭环；
2. 如遇到数值差异，再保留 SoPlex LP 路线；
3. 后续可支持 optlang/cobra 工具链，但不要一开始引入过重依赖。

## 阶段 5：复现 glucose 条件

目标：复现 `local_pichia_ref_glc.m` 的核心逻辑。

任务：

- Python 实现 `setMediaPP` 的 glucose 相关部分；
- 设置：
  - `Ex_glc_D`
  - `Ex_o2`
  - `Ex_glyc`
  - `Ex_meoh`
  - `BIOMASS`
  - `BIOMASS_glyc`
  - `BIOMASS_meoh`
- 固定 `mu`；
- 生成或直接求解 LP；
- 解析 objective。

验收：

- Python 和 MATLAB 在同一 `mu` 下结果一致；
- 至少完成 `mu=0.10`；
- 输出对非计算机用户可读。

## 阶段 6：移植 `writeLPGlc` 核心约束

目标：从普通 FBA 进入 pcSec/proteome-constrained 仿真。

任务：

- 阅读并拆解 `Code/pcSecPichia/CoreFunction/writeLPGlc.m`；
- 将约束拆成 Python 函数：
  - metabolic enzyme coupling；
  - secretion enzyme coupling；
  - total enzyme mass constraint；
  - unmodeled protein constraint；
  - unmodeled ER protein constraint；
  - ribosome / proteasome / degradation constraints；
  - optional misfolding constraints；
- 每个约束函数输出结构化 LP row，而不是直接拼字符串；
- 保留 LP 导出，便于和 MATLAB 文本 diff。

验收：

- Python 生成 LP 的变量数、约束数、关键约束名称与 MATLAB 对齐；
- 对同一输入，SoPlex objective 与 MATLAB LP 接近；
- 可以定位每条约束来自哪个模型模块。

风险：

- 这是最关键也是最容易出错的阶段；
- 不建议和 `addTargetProtein` 同时做，应该先让参考模型 LP 对齐。

## 阶段 7：目标蛋白添加

目标：移植 `addTargetProtein` 相关能力，使 OPN/hLF 可作为目标蛋白进入模型。

任务：

- 阅读并拆解 `Code/pcSecPichia/secPart/addTargetProtein.m`；
- 支持 fakeProteinInfo / target input schema；
- 添加 target protein translation / translocation / maturation / secretion reactions；
- 生成 target exchange reaction；
- 生成 target enzymedata；
- 调用或移植 `SimulateRxnKcatCoef`；
- 与 `CombineEnzymedata` 结果对齐。

验收：

- OPN 候选能生成与 MATLAB 等价的 target exchange；
- hLF 在参数完整后也能进入同一流程；
- target exchange 固定 production ratio 后，模型可求解；
- Python LP 与 MATLAB `local_opn_pichia_glc.m` 对齐。

关键边界：

- hLF 的二硫键、糖基化等参数必须确认；
- 参数未确认时不允许静默进入推荐结果。

## 阶段 8：单基因敲除筛选

目标：给出模型预测的 KO 候选。

任务：

- 解析 `grRules/rxnGeneMat`；
- 实现 gene deletion；
- 对每个 gene 批量求解；
- 固定最低生长要求；
- 记录目标蛋白通量变化；
- 过滤生长不可行或严重生长损伤候选。

输出字段：

```text
gene_id
gene_name
affected_reactions
baseline_target_flux
ko_target_flux
target_flux_change_percent
growth_feasible
growth_margin
recommendation_level
risk_note
```

验收：

- 能对 pcSecPichia 模型执行单基因 KO screen；
- 输出至少包含候选排序；
- 每个候选能解释影响了哪些 reaction；
- 对不可行结果有清晰说明。

## 阶段 9：过表达 / 容量放宽筛选

目标：给出模型预测的 OE 或分泌 capacity 候选。

任务：

- 定义 capacity perturbation，而不是简单 gene overexpression；
- 支持以下扰动：
  - enzyme kcat multiplier；
  - secretion complex capacity multiplier；
  - machinery complex capacity multiplier；
  - degradation / folding / translocation 模块放宽；
- 从 `enzymedata.subunit` 映射到 gene/subunit；
- 输出单基因或组合候选。

输出字段：

```text
perturbation_id
module
enzyme_or_complex
subunit_genes
capacity_multiplier
baseline_target_flux
perturbed_target_flux
target_flux_change_percent
growth_feasible
recommendation_level
interpretation
wet_lab_note
```

验收：

- 能识别目标蛋白分泌提升最敏感的 enzyme/complex/module；
- 能映射到 subunit genes；
- 能明确区分“单基因建议”和“复合体组合建议”；
- 不把模型预测写成实验定论。

## 阶段 10：瓶颈解释层

目标：把模型输出翻译成生物学专家能看懂的解释。

任务：

- 从 LP 结果提取关键约束余量或 tight constraints；
- 标记瓶颈属于：
  - 代谢底物；
  - 总蛋白分配；
  - ER folding；
  - Golgi processing；
  - translocation；
  - degradation / misfolding；
  - ribosome / translation；
- 给出中文解释和下一步实验建议。

验收：

- 每个推荐 KO/OE 候选都有解释；
- 页面不要求用户理解 LP、dual value 或 MATLAB；
- 下载结果中保留原始数值和解释文本。

## 阶段 11：分泌模型网络可视化

目标：让生物学专家能看到类似代谢网络的分泌模型图，但避免展示全量复杂网络。

原则：

- 不画全量 GEM 大图；
- 只画目标蛋白相关的分泌路径、瓶颈模块和扰动解释；
- 图是解释工具，不是装饰图；
- 图中所有节点都能追溯到模型 reaction、enzyme、gene 或 constraint。

三层图：

### 第一层：分泌路径总览

展示：

```text
目标蛋白
翻译
ER translocation
folding / disulfide
N-glycosylation
O-glycosylation
Golgi processing
vesicle transport
extracellular secretion
```

用途：

- 让用户理解目标蛋白在模型里经历哪些步骤；
- 对比 OPN、hLF 等不同蛋白的 PTM 负担。

### 第二层：瓶颈模块图

展示：

- 哪些模块约束接近饱和；
- 哪些模块被模型判定为产量限制因素；
- 不同 leader 或目标蛋白之间瓶颈是否不同。

颜色建议：

```text
绿色：余量充足
黄色：接近瓶颈
红色：强瓶颈
灰色：未参与或证据不足
```

### 第三层：基因/复合体扰动图

展示：

```text
候选基因 -> enzyme/complex -> reaction/module -> target flux impact
```

用途：

- 解释为什么推荐某个 KO/OE 候选；
- 区分单基因候选和复合体组合候选；
- 帮助实验同事设计组合过表达。

任务：

- 新增 `app/core/graph_models.py`；
- 新增 `app/services/secretion_graph.py`；
- 从 target simulation result 构建 pathway graph；
- 从 KO/OE screening result 构建 overlay；
- 新增 Streamlit 网络图页面；
- 提供 PNG/SVG/JSON 下载。

验收：

- OPN 和 hLF 都能生成分泌路径图；
- KO/OE 候选能在图上高亮；
- 图中节点点击后能看到 reaction/gene/enzyme 证据；
- 页面说明“这是模型解释图，不是完整细胞代谢网络”。

## 阶段 12：Streamlit UI 集成

目标：让生物学专家可以使用 Python engine。

页面建议：

- 目标蛋白管理；
- leader 候选管理；
- 基础分泌仿真；
- KO 筛选；
- OE/容量筛选；
- 瓶颈解释；
- 分泌网络可视化；
- 结果下载。

交互原则：

- 默认隐藏高级参数；
- 每个页面说明输入、输出、适用边界；
- 所有模型预测都标注“需要湿实验验证”；
- 下载 CSV/XLSX/FASTA/JSON；
- 保留原始 LP 和求解日志供研发人员追溯。

验收：

- 用户可选择 OPN 或 hLF；
- 用户可选择 leader 候选；
- 用户可运行基础仿真；
- 用户可运行 KO/OE screen；
- 用户可查看分泌路径和瓶颈网络图；
- 用户可下载结果；
- 页面不再要求用户启动 MATLAB。

## 9. 模块边界和依赖方向

为了保证高内聚、低耦合，依赖方向必须固定：

```text
ui -> services -> engines/adapters/core
engines -> core/adapters
services -> core
adapters -> external systems
core -> no app/service/ui imports
```

禁止：

- `core` import `streamlit`；
- `core` import MATLAB adapter；
- `engine` import UI；
- visualization service 直接调用 solver；
- target input provider 直接修改模型；
- KO/OE screening 直接读写 Streamlit session state。

允许：

- service 编排多个 engine/adapters；
- adapter 包装外部系统；
- UI 调用 service；
- graph renderer 只消费 graph model。

## 10. 测试策略

### 单元测试

- target input schema；
- `.mat` loader；
- bounds 修改；
- GPR 解析；
- gene deletion；
- enzyme capacity perturbation；
- LP matrix 构造；
- SoPlex/HiGHS 输出解析。
- secretion graph model 构造；
- target registry 扩展；
- graph overlay 高亮。

### 对齐测试

- Python 读取模型维度与 MATLAB 一致；
- Python 生成 LP 与 MATLAB baseline 关键约束一致；
- objective value 在容许误差内一致；
- OPN smoke 对齐；
- hLF smoke 在参数确认后加入。

### UI 测试

- 页面能加载目标蛋白；
- 缺失 hLF PTM 参数时给出中文提示；
- KO/OE 结果表可下载；
- 结果解释不乱码；
- Streamlit health check 通过。
- 分泌网络页面能渲染 OPN/hLF 示例图。

## 11. 依赖建议

建议新增 Python 依赖：

```text
numpy
scipy
pandas
h5py
openpyxl
pydantic
```

可选依赖：

```text
cobra
optlang
highspy
networkx
```

建议先使用 `scipy.optimize.linprog(method="highs")` 或外部 SoPlex，避免一开始把 Python COBRA 工具链也变成新的复杂依赖。

网络可视化可选依赖：

```text
plotly
pyvis
streamlit-cytoscapejs
```

建议首轮优先使用 `networkx` 生成 graph model，再用 Plotly 或 Cytoscape.js 渲染。

## 12. 风险和应对

### 风险 1：MATLAB 与 Python 数值不一致

应对：

- 保留 MATLAB baseline；
- 每阶段只对齐一个 smoke case；
- 保存 LP 文本用于 diff；
- 不直接跳到 KO/OE screen。

### 风险 2：hLF 参数不完整

应对：

- hLF 作为 target preset 时标记参数来源；
- 缺失二硫键/糖基化参数时不允许进入正式推荐；
- 支持“草案模式”和“正式模型模式”。

### 风险 3：过表达解释过度

应对：

- 输出“模型预测”；
- 将单基因和复合体组合分开；
- 把 wet-lab validation note 作为强制字段。

### 风险 4：范围膨胀

应对：

- 第一版只做毕赤酵母；
- 第一版只做 glucose；
- 第一版只做 OPN/hLF；
- 第一版不做论文复现。

### 风险 5：可视化图过于复杂

应对：

- 不画全量代谢网络；
- 默认只展示目标蛋白相关分泌路径；
- 使用三层图逐步展开；
- 高级 reaction/gene 细节放到点击详情或下载文件。

### 风险 6：架构重新耦合

应对：

- 每个新增模块先写接口和测试；
- 禁止 UI 直接调用 adapter；
- engine 返回稳定 result schema；
- graph model 与具体渲染库分离；
- code review 时把依赖方向作为检查项。

## 13. 推荐提交拆分

建议按小提交推进：

1. `docs(pichia): add python refactor plan`
2. `refactor(engine): add pichia engine interface`
3. `refactor(targets): generalize target protein inputs`
4. `feat(pichia): load pcSecPichia matlab model`
5. `feat(pichia): solve baseline glucose LP`
6. `feat(pichia): reproduce target protein simulation`
7. `feat(screening): add gene knockout screen`
8. `feat(screening): add capacity perturbation screen`
9. `feat(graph): add secretion network graph model`
10. `feat(ui): expose pichia python engine workflow`

每个提交都应有对应测试，不把大迁移压成一个无法审查的提交。

## 14. 第一批具体任务清单

建议下一轮从这里开始：

1. 新增 `app/engines/base.py`；
2. 新增 `app/engines/matlab_pichia_engine.py`；
3. 将当前 OPN MATLAB 调用包进 engine；
4. 新增 `app/core/target_inputs.py`；
5. 将 OPN 输入迁到通用 target schema；
6. 新增 hLF target draft；
7. 新增 `app/adapters/mat_loader.py`；
8. 读取 `Model/pcSecPichia.mat` 并输出模型摘要；
9. 编写模型摘要测试；
10. 保存 MATLAB OPN smoke baseline；
11. 新增 `app/core/graph_models.py` 的空模型和测试；
12. 明确 UI/service/engine/core/adapters 的依赖规则。

第一批验收：

- 不改变现有 Streamlit 启动命令；
- 不改变当前 OPN smoke 行为；
- `python -m pytest -q` 通过；
- 有通用 target input schema；
- 有 pcSecPichia 模型读取摘要；
- 有 graph model 的最小数据结构；
- 后续可以开始 Python LP 求解。

## 15. 决策摘要

这次重构的核心不是“去 MATLAB 化”本身，而是让 pcSecPichia 变成一个可以服务研发决策的 Python 工具。

优先顺序应为：

1. 可部署；
2. 可复现；
3. 可解释；
4. 可扩展到 OPN/hLF；
5. 可筛选 KO/OE 候选；
6. 可用网络图解释模型预测；
7. 最后才考虑完整替代 MATLAB 论文工作流。

只要坚持这个顺序，迁移工作量就是可控的，而且每个阶段都能产生实际价值。
