# pcSecPichia 当前需求与架构

状态：active

## 用户需求

当前用户是生物研发组同事。文档和系统设计都应服务于一个实际目标：提高毕赤酵母中目标蛋白的产量，并把模型结果转化为可讨论、可复核、可安排实验的候选方向。

当前目标蛋白分为两类：

- 人乳铁蛋白 hLF。
- 骨桥蛋白 OPN。

已经完成或已有外部输入：

- 目标蛋白序列和构建设计。
- 密码子优化。
- 信号肽筛选。
- 部分候选分泌工程基因、反应代理和外部注释证据。

接下来要支持的核心工作：

- KO/OE 模拟：帮助研发同事比较候选敲除和过表达方向。
- MATLAB 原有功能迁移：把既有 pcSecPichia 能力迁入 Python engine，形成更容易在本地 UI/API 中使用的工作流。
- 结果解释：区分模型可执行性、数据库注释、外部表型证据和实验建议等级。

## 科学边界

- 当前输出用于模型内相对排序和解释，不预测绝对 mg/L 产量。
- KO/OE 结果不是实验成功率承诺。
- OE gene 当前通常是 reaction-level capacity proxy，不等同于完整基因表达调控模拟。
- UniProt、KEGG、common name 等数据库注释只说明注释来源，不能单独证明表型效果。
- GPR 可执行只说明模型能运行 gene deletion 或 reaction proxy，不能自动升级为实验校准结论。
- 表型证据只来自明确 curated/fixture evidence，且 KO 与 OE 必须按 intervention 分开判断。
- `experiment_calibrated` 只用于同 host、同 target/context、同 intervention 的高置信表型证据。

## 当前系统分层

```text
Streamlit UI / FastAPI facade
  -> app/services
    -> python_pichia engine
      -> loading / media / targets / secretion_plan
      -> constraints / simulation / screens / analysis / reports
      -> gene evidence / gene catalog / rule overlay
```

职责边界：

- `python_pichia/` 是核心 engine，放模型加载、培养基条件、目标蛋白构建、约束、求解、KO/OE screen、分析、报告和证据分级逻辑。
- `app/services/` 是 facade，负责请求响应、后台任务、缓存路径、错误和 warning 汇总，不写核心模型算法。
- `app/ui/` 是 Streamlit 展示和表单，不直接承载科学判断。
- `app/api/` 是 experimental API facade，不绕过 service/engine 边界。
- `Code/`、`Model/`、`Enzymedata/`、`Results/` 是原始 MATLAB/reference 目录，不作为 Python 修改对象。
- `local_runs/` 只放本地运行产物、缓存和验证证据，不提交为源码。

## 数据与结果目录边界

- `Data/`、`Model/`、`Enzymedata/` 是稳定科学输入资产目录，只读使用，新增内容必须明确说明科学来源和提交理由。
- `Results/` 保留为 legacy MATLAB results，只读参考历史模拟/分析结果，不作为当前 Python 或 Streamlit 输出目录。
- `local_runs/` 是当前 Python corrected pipeline、Streamlit 工作台、MATLAB harness、LP diff、缓存和本地验证证据的默认输出位置。
- 新生成的 LP、solver output、CSV/JSON/Markdown 报告和缓存不得进入 `Data/`、`Model/`、`Enzymedata/` 或 `Results/`，除非另行声明为科学资产。

## 当前能力

- 支持 OPN、hLF 和 custom target 的 Python corrected 分泌仿真。
- 支持基础培养基和碳源条件，包括 mixed-carbon probe。
- 支持 small-grid growth tradeoff 和目标蛋白成本解释。
- 支持 KO/OE 候选预览、screen rows、yield recommendation 和报告透传。
- 支持 KO/OE phenotype evidence tier：
  - `model_executable`
  - `evidence_supported`
  - `experiment_calibrated`
  - `manual_review_required`
  - `not_recommended_growth_risk`
- 支持 Streamlit 本地工作台和最小 API facade。

## 当前非目标

- 不做三物种 MATLAB 项目全量迁移。
- 不重建论文全部 figure。
- 不自动生成新 MATLAB baseline，除非单独立项。
- 不做全模型 KO/OE 大规模筛选，除非后续把计算和解释边界单独定义清楚。
- 不把本地缓存、LP 输出、求解结果或 UI 运行产物提交进源码。
