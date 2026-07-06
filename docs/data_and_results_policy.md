# pcSecPichia 数据与结果治理策略

状态：active

## 目录职责

- `Data/` 只放稳定输入、参考数据、外部注释、目标蛋白输入和人工确认的候选数据；不放运行输出。
- `Model/` 只放人工确认的 GEM 和 pcSec 模型资产；不放脚本生成的临时模型。
- `Enzymedata/` 只放人工确认的酶约束资产；不放本地求解缓存。
- `Results/` 是 legacy MATLAB results，只读参考历史结果，不是当前 Python 或 Streamlit 的默认输出目录。
- `local_runs/` 是当前 Python、Streamlit、MATLAB harness、LP diff、缓存和验证证据的默认运行产物目录，并保持 ignored。
- `python_pichia/local_runs/` 若被局部工具创建，也视为运行产物目录，并保持 ignored。

## 提交规则

- 新生成的 LP、solver output、MATLAB harness output、Streamlit run result、缓存和盘点文件默认进入 `local_runs/`。
- 新增大文件前必须说明来源、是否可再生成、是否应进入 Git LFS，以及为什么不是运行产物。
- `Data/`、`Model/`、`Enzymedata/`、`Results/` 的新增或修改必须被明确声明为科学资产变更。
- 第一轮治理不移动历史 MATLAB 数据、不清理 Git 历史、不改变 legacy browser 对 `Results/` 的只读访问。

## 当前盘点

本 checkpoint 生成的只读盘点文件位于：

- `local_runs/data_inventory_tracked_files.txt`
- `local_runs/git_status_ignored_snapshot.txt`
- `local_runs/data_inventory_by_extension.csv`
- `local_runs/data_path_audit_matches.txt`

这些文件是本地治理证据，不提交到 Git。

## 后续迁移边界

如需迁移历史 `Results/`、启用 Git LFS、发布 GitHub Release artifact 或外部存储，应单独立项并先做 checkpoint。迁移前不得删除或移动 legacy MATLAB 结果。
