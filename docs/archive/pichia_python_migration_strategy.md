# pcSecPichia Python 移植与研发用途决策记录

日期：2026-06-22

## 1. 背景

当前 `pcSecYeastSpecies` 项目包含三种酵母的 pcSec 模型：

- *Saccharomyces cerevisiae* / 酿酒酵母
- *Komagataella phaffii* / 毕赤酵母
- *Kluyveromyces marxianus* / 马克斯克鲁维酵母

原始建模和仿真主体仍然是 MATLAB + COBRA/RAVEN + SoPlex。Python 部分目前主要负责：

- Streamlit 前端展示
- 结果浏览
- 健康检查
- 日志和 SoPlex 输出解析
- OPN 输入生成接口
- PichiaCLM / SigScout 等外部工具的衔接

近期讨论的核心问题是：是否有必要把这个项目移植到 Python，以及移植后能否服务真实研发需求。

## 2. 总体结论

不建议把整个 `pcSecYeastSpecies` 全量翻译成 Python。

更合理的目标是：

> 构建一个 **毕赤酵母 pcSecPichia 目标蛋白分泌优化 Python 引擎**。

这个引擎不以复现全部论文图为目标，而是面向研发问题：

- 给定目标蛋白，例如 OPN、hLF；
- 给定信号肽/leader 候选；
- 在毕赤酵母分泌模型中评估目标蛋白分泌能力；
- 进一步筛选可能提高产量的敲除/过表达候选基因或分泌通路模块。

## 3. 为什么只做毕赤酵母

如果全量迁移三种酵母，需要处理 `Code/` 下约 253 个 MATLAB 文件，以及大量论文复现脚本。

如果只保留毕赤酵母，核心范围收缩到 `Code/pcSecPichia/`，约 75 个 MATLAB 文件。更重要的是，我们还可以进一步不迁移所有文件，而是优先迁移真实研发需要的运行链路。

可先舍弃：

- 酿酒酵母和马克斯克鲁维酵母模型；
- 跨物种比较；
- 全部论文 figure 复现；
- 温度依赖、人源化糖基化等非首轮功能；
- 完整 `buildModel_pcSecPichia.m` 重建流程。

优先保留：

- 已构建好的 `Model/pcSecPichia.mat`；
- `Enzymedata/pcSecPichia/` 中的酶和分泌 machinery 数据；
- glucose 条件下的目标蛋白分泌仿真；
- OPN、hLF 以及后续目标蛋白输入；
- 基因扰动筛选能力。

## 4. 为什么不能只做 OPN

只做 OPN 会把系统做成一次性工具，不适合后续研发。

实际需求已经包括 hLF，因此应抽象成通用目标蛋白输入：

```text
target_id
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
```

这样 OPN、hLF、HSA、PHO 或未来其它蛋白都可以复用同一套 Python runtime。

当前 `app/core/opn_inputs.py` 已经为 OPN 留出了 provider-style 输入边界。后续建议把它推广为更通用的 target protein input schema，而不是继续维护 OPN 专用接口。

## 5. 移植后能否回答敲除/过表达问题

可以，但要注意模型结论的边界。

移植后可以输出：

- 模型预测的敲除候选；
- 模型预测的过表达或容量放宽候选；
- 候选改造对目标蛋白分泌通量的影响；
- 候选改造对生长率的影响；
- 相关反应、酶复合体、分泌 machinery 的解释。

但这些结果不能直接等同于湿实验产量结论。它们应作为实验设计前的计算筛选和排序依据。

## 6. 可行性依据

本地检查显示，`Model/pcSecPichia.mat` 中保留了基因和反应关联字段：

```text
genes
rules
grRules
rxnGeneMat
```

这意味着可以基于 GPR / gene-reaction mapping 做基因敲除模拟。

`enzymedataSEC_PP.mat` 和 `enzymedataMachine_PP.mat` 中保留了分泌与 machinery 信息：

```text
enzyme
subunit
subunit_stoichiometry
kcat
proteins
proteinMWs
proteinLength
proteinPST
proteinExtraMW
```

这意味着可以把模型瓶颈映射回分泌通路复合体、subunit genes 和潜在过表达对象。

## 7. 敲除筛选怎么做

敲除相对更直接，推荐先实现。

基本流程：

1. 读取 pcSecPichia 模型；
2. 选择目标蛋白，例如 OPN 或 hLF；
3. 固定培养基、氧气、碳源、生长速率下限；
4. 设定目标函数为最大化目标蛋白 exchange；
5. 对每个候选基因做单基因删除；
6. 根据 `grRules/rxnGeneMat` 找到受影响反应；
7. 将受影响反应上限/下限置零或按 GPR 规则更新；
8. 重新求解；
9. 比较目标蛋白分泌通量、生长影响和可行性。

输出示例：

```text
基因 X 删除后，hLF 分泌目标值提高 18%，生长约束仍可满足。
基因 Y 删除后，分泌提高 25%，但模型生长不可行，不建议作为首轮候选。
```

优先输出字段：

```text
gene_id
affected_reactions
target_flux_change_percent
growth_feasible
objective_value
risk_note
recommendation_level
```

## 8. 过表达筛选怎么做

过表达比敲除复杂，不能简单理解成把某个反应上限调大。

在 pcSec 模型中，很多候选对象是分泌 machinery 或酶复合体。单独过表达一个 subunit 可能不增加复合体容量，因为其它 subunit、总蛋白约束或 ER/Golgi 处理能力仍可能成为瓶颈。

因此建议分两层做。

第一层：容量放宽 / sensitivity screen

- 放宽某个 enzyme complex 的 kcat 或容量约束；
- 放宽某个分泌模块约束；
- 放宽 ER folding、Golgi processing、translocation、degradation 等模块；
- 观察目标蛋白分泌是否提高。

第二层：映射到基因

- 从 `enzymedata.subunit` 找到相关 subunit genes；
- 判断是单基因候选还是复合体组合候选；
- 输出“建议测试的基因或组合”，而不是武断给出单基因结论。

输出示例：

```text
模型提示 sec_xxx_complex 是 hLF 分泌瓶颈。
相关 subunit genes 为 A/B/C。
建议优先测试 A+B+C 组合过表达；单独 A 可能不足。
```

## 9. 推荐 Python 移植架构

建议引入一个可切换的 engine 层：

```text
app/engines/
  base.py
  matlab_engine.py
  python_pichia_engine.py
```

接口示例：

```text
simulate_target_protein(...)
screen_gene_knockouts(...)
screen_capacity_overexpression(...)
explain_bottlenecks(...)
```

这样前端和服务层不需要知道底层是 MATLAB 还是 Python。早期可以保留 MATLAB engine 作为参考标准，Python engine 逐步替换。

## 10. 推荐迁移阶段

### 阶段 0：定义验收基线

选择一个最小可对齐场景：

- 毕赤酵母；
- glucose 条件；
- 固定 `mu`；
- 目标蛋白为 OPN 或 hLF；
- 输出 LP 和 SoPlex 求解结果；
- 与 MATLAB 的 `local_opn_pichia_glc.m` 或等价 hLF 脚本对齐。

### 阶段 1：Python 读取模型并做基础求解

实现：

- 读取 `Model/pcSecPichia.mat`；
- 映射 `rxns/mets/S/lb/ub/c/b`；
- 实现 `change_rxn_bounds`；
- 支持基础 LP 求解；
- 与 MATLAB 生成的基础 LP/objective 对齐。

### 阶段 2：复现 glucose 条件目标蛋白仿真

优先复现：

- `setMediaPP`
- glucose / oxygen / glycerol / methanol bounds；
- biomass 固定；
- target protein exchange；
- `writeLPGlc` 的核心约束。

### 阶段 3：通用目标蛋白输入

把 OPN 专用输入升级为通用 target protein schema。

至少支持：

- OPN
- hLF

此阶段不要追求全部目标蛋白，只要保证 schema 能扩展。

### 阶段 4：单基因敲除筛选

实现：

- GPR 解析；
- 单基因删除；
- 反应失活；
- 批量求解；
- 候选排序；
- 对生长损伤做过滤。

### 阶段 5：过表达 / 容量放宽筛选

实现：

- enzyme / complex / secretion machinery 级别的容量放宽；
- subunit gene 映射；
- 输出单基因或组合过表达建议；
- 标注模型证据等级和湿实验风险。

### 阶段 6：Streamlit 展示

前端应面向实验和研发人员，而不是 MATLAB 用户。

建议页面：

- 目标蛋白管理：OPN、hLF、其它目标蛋白；
- 候选 leader 管理；
- 基础分泌仿真；
- 敲除候选筛选；
- 过表达/瓶颈模块筛选；
- 结果解释和下载。

## 11. 首轮不建议做的事

首轮不要做：

- 三物种全量迁移；
- 完整论文复现；
- 全部 figure 脚本 Python 化；
- 温度依赖 ETC 全量迁移；
- humanized glycosylation 全量迁移；
- 直接承诺真实发酵产量。

这些可以作为后续增强，而不是第一阶段验收标准。

## 12. 推荐的下一步

下一步最有价值的是写一个具体实施计划：

1. 新增 `PcSecPichiaModel` Python 数据结构；
2. 新增 `.mat` 模型加载器；
3. 新增 `PichiaTargetProteinInput` 通用 schema；
4. 把 OPN 输入迁到通用 schema；
5. 补 hLF 成熟序列和 PTM 参数；
6. 实现 `python_pichia_engine` 的基础 glucose 仿真；
7. 与 MATLAB OPN smoke 对齐；
8. 再进入敲除/过表达筛选。

最小验收标准：

- 不启动 MATLAB，也能读取 `pcSecPichia.mat`；
- Python 能生成或求解一个与 MATLAB smoke 等价的 glucose LP；
- OPN 和 hLF 都能作为目标蛋白输入；
- 结果能说明目标蛋白分泌通量、可行性和主要瓶颈；
- 后续能扩展到 gene KO / overexpression screen。

## 13. 关键边界

Python 移植的目标不是让模型“更科学”，而是让模型更可用、更易部署、更容易被研发团队和生物学专家使用。

模型给出的敲除/过表达候选应被表述为：

> 模型预测的优先实验候选。

而不是：

> 已证明能提高真实发酵产量的改造方案。

最终结论仍需要小试、发酵数据和产物检测验证。
