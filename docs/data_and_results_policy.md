# pcSecPichia 数据与结果治理策略

状态：active  
最后更新：2026-07-07

## 目录职责

- `Data/`：稳定输入、参考数据、外部注释、目标蛋白输入和人工确认的科学数据。
- `Model/`：人工确认的 GEM / pcSec 模型资产。
- `Enzymedata/`：人工确认的酶约束资产。
- `Results/`：legacy MATLAB results，只读参考历史结果，不作为当前 Python / Streamlit 默认输出。
- `local_runs/`：当前 Python、Streamlit、MATLAB harness、LP diff、cache 和验证证据的默认运行产物目录，保持 ignored。
- `docs/archive/`：历史计划、阶段性验证记录和已被 active 文档吸收的长文档。

## 提交规则

- 新生成的 LP、solver output、BLAST cache、MATLAB harness output、Streamlit run result、Markdown report 默认进入 `local_runs/`。
- 新增大文件前必须说明来源、是否可再生成、是否应进入 Git LFS，以及为什么不是运行产物。
- `Data/`、`Model/`、`Enzymedata/`、`Results/` 的新增或修改必须被明确声明为科学资产变更。
- 本轮 BLAST/RBH 同源映射 cache 不写入 `Data/`，先作为 `local_runs/` 产物验证。
- Streamlit 同源审计页面默认只读 `local_runs/` 或未来人工提升的稳定 cache；页面打开、筛选和导出不得默认运行 BLAST、联网查询或生成大文件。
- 若未来要把同源 crosswalk 升级为稳定资产，应单独建立 checkpoint，并说明来源、版本、生成命令和人工复核状态。

## 保护检查

每个涉及模型、screen、cache 或报告的任务结束前都应运行：

```powershell
git diff --name-only -- Code Model Enzymedata Results requirements.txt python_pichia\pyproject.toml
```

预期为空，除非该任务明确要求修改这些边界文件。

## 归档规则

归档不是删除历史，而是把不再作为当前入口的长文档移到 `docs/archive/`。

适合归档：

- 已完成阶段验证记录。
- 已被当前架构文档吸收的旧计划。
- 长篇排查记录，但当前下一步不再依赖其细节。
- 与当前任务无关的历史候选或旧设计。

不应归档：

- 当前架构入口。
- 下一步计划。
- 数据边界和提交规则。
- 当前任务的正式设计文档。

## 后续迁移边界

如需迁移历史 `Results/`、启用 Git LFS、发布 GitHub Release artifact 或接入外部存储，应单独立项并先做 checkpoint。迁移前不得删除或移动 legacy MATLAB 结果。
