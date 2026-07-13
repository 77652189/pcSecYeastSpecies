# pcSecPichia 数据与结果治理策略

状态：active  
最后更新：2026-07-13

## 目录职责

- `Data/`：人工确认、可追溯、允许长期维护的稳定输入和 curated scientific data。
- `Model/`：人工确认的 GEM / pcSec 模型资产。
- `Enzymedata/`：人工确认的酶容量和资源约束资产。
- `Results/`：legacy MATLAB results，只读参考，不作为当前 Python / Streamlit 输出目录。
- `local_runs/`：运行产物、临时缓存、外部下载、smoke、报告和待复核实验导入，默认 ignored。
- `docs/archive/`：已完成计划、阶段验证和被 active 文档吸收的历史设计，默认不进入公开版本控制。

## 实验反馈数据

实验数据分为三层：

1. 原始导入：仪器导出、研发记录和未清洗表格进入 `local_runs/experiment_feedback/inbox/`。
2. 标准化缓存：完成 schema、单位、重复和 prediction linkage 检查后，进入 `local_runs/experiment_feedback/validated/`。
3. 稳定 curated 数据：只有人工确认来源、权限、脱敏、单位和上下文后，才可单独 checkpoint 提升到 `Data/` 下的正式目录。

不得把自由文本“提高/降低”当作唯一实验结果。至少保留 target、host、intervention、培养条件、测量方法、单位、重复、误差和原始来源。

失败、无效、不可测和阴性结果必须保留，不能只收集成功实验。

## 外部数据与模型

- BLAST/RBH、在线数据库响应、外部 GEM、MEMOTE 报告和 GPR mapping 默认进入 `local_runs/`。
- 外部数据必须保留 source URL、version、retrieved_at、query、hash、license/provenance 和 warning。
- 外部名称、同源关系或 GPR 不得自动覆盖内部 `gene_id` 和当前模型规则。
- 将外部模型或 cache 提升到 `Data/Model/Enzymedata` 前，必须单独 review 和 checkpoint。

## OE capacity 参数与映射

- gene-enzyme-reaction mapping、kcat、分子量、复合体亚基、基线丰度和剂量映射都必须保留单位、来源、版本、模型指纹、置信度和 warning。
- 当前 pcSecPichia 模型 GPR 与本地 enzyme data 是“当前模型可执行性”的权威来源；外部 iPichia/ecPichia、UniProt、BRENDA、SABIO-RK 或同源转移只作为候选证据。
- promoter、copy number、induction mode 等实验标签没有经过审核的 dose mapping 时，只能保留为类别输入，不能自动换算成单一 expression multiplier。
- 缺失参数使用显式区间和 low/nominal/high 场景；不得用未标注的默认值伪装成测量真值。
- 外部下载、参数候选、mapping audit 和 Phase 2 screen 输出默认写入 `local_runs/oe_capacity/`。
- 只有人工复核 license、provenance、映射和参数后，才能通过独立 checkpoint 提升到稳定科学资产目录。

## LLM 数据边界

- LLM 只读取程序生成的 fact pack，不直接遍历运行目录。
- 只有用户明确触发时才调用外部 LLM API。
- fact pack 应只包含完成报告所需的模型结果和证据，不包含 API key、环境配置或无关实验原始文件。
- 最终报告必须同时通过程序 validator 和 Judge。
- LLM 不能修改模型、recommendation tier、curated evidence 或实验记录。

## 密钥和本地配置

- `.env`、`.env.*`、API key、认证文件和个人路径不得提交。
- 可提交不含真实值的 `.env.example`。
- 日志、manifest、异常信息和测试输出不得打印密钥或完整认证 header。

## 提交规则

- LP、solver output、BLAST cache、外部下载、Streamlit run、LLM report 和 smoke 默认进入 `local_runs/`。
- 大文件进入 Git 前必须说明来源、可再生成性、license、人工复核状态和为什么不是运行产物。
- `Data/Model/Enzymedata/Results` 的任何修改都必须明确声明为科学资产变更。
- 不提交包含真实实验人员身份、客户信息或未经批准的原始实验数据。

## 保护检查

每个涉及模型、筛查、证据、实验反馈或报告的任务结束前运行：

```powershell
git diff --name-only -- Code Model Enzymedata Results requirements.txt python_pichia\pyproject.toml
```

预期为空，除非本轮明确声明相应边界变更。

同时检查：

```powershell
git status --short
git diff --check
git ls-files -- .env '.env.*'
```

## 归档规则

适合归档：已完成阶段计划、一次性可行性验证、长排查记录和已被当前文档吸收的函数设计。

active 根目录只保留当前架构、下一步计划、数据治理和索引。归档文档不能成为执行当前任务的必要前置条件。
