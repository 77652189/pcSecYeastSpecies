# pcSecPichia 当前架构与需求

状态：active  
最后更新：2026-07-09

## 目标

当前 Python 工作流服务于一个具体研发问题：帮助研发同事围绕 Pichia 中目标蛋白产量提升，提出、解释和复核 KO/OE 改造候选。

系统输出用于模型内相对比较、候选排序和研发讨论，不承诺真实发酵产量、mg/L 绝对值或实验成功率。

## 当前主线

- 目标蛋白：项目内置参考目标，以及用户输入的 custom target。
- 上游输入：密码子优化、信号肽筛选、目标蛋白序列和构建设计。
- 核心问题：在现有 pcSecPichia 模型中比较分泌负担、生长权衡、KO/OE 改造方向和证据等级。
- 当前重点：把模型可执行性、外部同源证据、表型证据和实验复核边界分开表达。
- 当前 BLAST/RBH 产品目标：把本地同源 cache 提升为 Streamlit 中可查看、筛选、导出和解释的“基因命名标准化 + 同源规则迁移评估”功能；CLI / scripts 只作为证据生成层，不是最终用户入口。
- 当前外部数据库目标：通过受控联网 fetcher 拉取 UniProt / NCBI / SGD external reference cache，并从 yeast-GEM / BiGG / BioModels / ModelSEED 等外部模型资源提取候选 reaction/GPR evidence，用于命名标准化、KO/OE gene function 补充和同源规则迁移校对；页面加载、核心仿真和 KO/OE screen 不依赖实时网络状态。该主线已完成 Round 1-9，小批量 smoke 产物仍保留在 ignored `local_runs/`，后续需要人工复核后再决定是否提升为稳定科学资产。

## 架构分层

```text
Streamlit UI / API facade
  -> app/services
    -> python_pichia engine
      -> loading / media / targets
      -> constraints / simulation / shadow_lp
      -> screens / analysis / reports
      -> gene evidence / gene catalog / homology evidence
```

职责边界：

- `python_pichia/` 是核心 engine，放模型加载、目标蛋白构建、约束、求解、KO/OE screen、证据分级和报告逻辑。
- `app/services/` 是 facade，负责请求响应、后台任务、缓存路径、错误和 warning 汇总，不实现核心科学判断。
- `app/ui/` 只做 Streamlit 展示和表单，不直接承载模型判断。
- `Code/`、`Model/`、`Enzymedata/`、`Results/` 是 MATLAB/reference 资产目录，默认只读。
- `local_runs/` 存放本地运行产物、缓存、报告和验证证据，默认不提交。

数据与结果目录边界：`Results/` 保留为 legacy MATLAB results，只读参考，不作为当前 Python / Streamlit 输出目录；当前运行产物统一进入 ignored `local_runs/`。

## 已有能力

- 加载 pcSecPichia 参考模型输入。
- 构建内置参考目标、custom target 的分泌路径。
- 支持培养基和碳源条件设置，包括 mixed-carbon probe。
- 支持固定生长率分泌能力求解、growth tradeoff、protein cost summary。
- 支持 KO/OE 候选预览、screen rows、yield recommendation 和 report 输出。
- 支持 gene evidence、gene catalog、phenotype evidence 和 recommendation tier。
- 支持 genome-wide / catalog-level KO/OE screen 工具。
- 支持 Shadow LP constrained solve 路径，用于承载 pcSec 语义约束的可替代求解后端。
- 已有 BLAST/RBH 本地 homology cache builder，可生成 homology cache、name audit、rule-transfer audit、summary，并可合并 external name crosscheck 字段。
- 已有 Streamlit 页面“基因命名与同源规则审计”，通过 `app/services` 只读 cache，支持查看、筛选、导出和解释同源命名/规则迁移结果。
- 在线外部数据库证据层已完成受控 fetcher、cache schema、命名校对、gene function evidence、external GPR candidate evidence、Streamlit 手动刷新/状态展示，以及 KO/OE screen/report/LLM fact pack 透传。

## COBRApy / Shadow LP 状态

当前不能说“COBRApy 完全替代原有逻辑”。更准确的状态是：

- `solve_shadow_secretion_capacity(...)` 已经作为并行入口实现。
- Shadow LP 使用结构化 constraint builders 和 `ScipyHighsBackend` 作为大模型默认求解后端。
- `CobraOptlangBackend` 保留为 tiny LP / 语义验证用途，不作为 full pcSec 默认后端。
- 内置参考目标的默认 fixed-growth capacity 已和 reference path 高精度对齐。
- 旧入口 `solve_secretion_capacity(...)` 仍然默认走 reference / 原 pcSec 路径。
- reference solver 目前仍可用于 validation / comparison boundary。

下一步若要推进替代，应先做 opt-in service/backend toggle，而不是立刻切换默认路径。

## 证据边界

- 数据库注释可信：UniProt / NCBI / SGD / annotation 字段只说明注释来源和命名一致性，不能单独证明表型效果。
- 模型 GPR 可执行：说明 KO 或 reaction proxy 能在模型中运行，不等于实验可行。
- OE reaction proxy：只是 reaction-level capacity proxy，不是完整 gene-level expression simulation。
- 表型证据：只由明确 curated phenotype evidence 决定，且 KO / OE 必须按 intervention 分开。
- 同源证据：BLAST/RBH 只能说明跨物种候选关系，不能直接变成模型可操作 gene。
- 命名标准化和同源规则迁移评估必须同时区分 sequence-level homology evidence、external/database name evidence、current model gene_index operability 和 manual review status。
- 联网证据必须保留 provenance、retrieved_at、source_version、query 和 warning；外部名称匹配不得覆盖内部 `gene_id` 或直接改变 `recommendation_tier`。
- 外部 GPR 规则和 reaction association 只能作为候选证据；只有映射到当前 Pichia GEM 的 reaction/gene 后，才可标记为 `model_gpr_executable`。
- `experiment_calibrated` 只用于 host、target/context、intervention 严格匹配的高置信表型证据。

## 非目标

- 不做三物种 MATLAB 项目全量重写。
- 不重建论文全部 figure pipeline。
- 不自动生成新 MATLAB baseline，除非单独立项。
- 不把本地缓存、LP 输出、solver output 或 UI 运行产物提交进源码。
- 不把模型输出解释成绝对产量、mg/L 或实验成功率。

## 当前文档入口

- [文档索引](README.md)
- [下一步计划](pichia_next_plan.md)
- [BLAST/RBH 同源映射架构](pichia_homology_crosswalk_architecture.md)
- [在线外部数据库证据层架构](pichia_online_external_reference_architecture.md)
- [数据与结果治理策略](data_and_results_policy.md)
