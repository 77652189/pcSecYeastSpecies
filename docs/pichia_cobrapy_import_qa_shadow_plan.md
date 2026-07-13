# COBRApy 导入、GEM QA 与 Shadow LP 产品化计划

状态：implemented_with_follow_up
最后更新：2026-07-13

## 当前验收状态

- COBRApy import probe 已实现 optional dependency、artifact cache 输入、结构计数、objective/GPR 语义检查和 libSBML 对照状态；真实 toy SBML 已在隔离 COBRApy 环境跑通。
- GEM basic QA 已实现，MEMOTE 适配改用其公开 `test_model` / `snapshot_report` API；默认开发环境未安装 MEMOTE，因此真实 MEMOTE 评分仍需在专用 QA 环境验收。
- Shadow LP cross-check 已有 core、service facade 和 Streamlit 入口；hLF 完整约束层 smoke 与 reference 对齐，相对差约 `2.85e-7`。
- 默认 KO/OE、report 和核心仿真不依赖 COBRApy/MEMOTE，运行产物继续写入 ignored `local_runs/`。
- 尚未完成任意 custom target 的 Shadow LP preparation；当前正式入口只支持已有 built-in target。该能力需要单独设计 prepared-target 输入契约，不能通过 UI 文本框伪装完成。

## 决策摘要

COBRApy 不作为 full pcSec 默认求解后端。当前推荐方向是精准引入：

1. 用 COBRApy 辅助外部 GEM/SBML 导入和语义核对。
2. 用 MEMOTE/COBRApy 做外部 GEM 离线 QA。
3. 将 Shadow LP cross-check 产品化为 Streamlit 可触发的验证能力。

这三件事服务于同一个目标：增强外部模型接入、模型质量审计和求解一致性验证，但不替换当前 `ScipyHighsBackend` / reference path 的默认生产语义。

## 共同边界

- `python_pichia` 承载核心模型、导入、QA、Shadow LP 和验证逻辑。
- `app/services` 只做 facade：读取状态、启动任务、汇总结果和透传警告。
- `app/ui` 只做 Streamlit 展示和用户触发，不实现科学判断。
- 不修改 `Code/`、`Model/`、`Enzymedata/`、`Results/`。
- 新的运行产物、外部模型 QA 报告、cross-check 报告和 smoke 输出写入 ignored `local_runs/`。
- COBRApy / MEMOTE 先作为 optional tooling，不进入默认求解路径；缺依赖时必须稳定返回 unavailable/skip，而不是让页面或核心仿真失败。
- 外部 GEM/GPR、MEMOTE 分数、COBRApy 导入结果和 Shadow LP 对照结果都不能自动提升 `recommendation_tier`，也不能声明 mg/L、绝对产量或实验成功率。

## 主线 A：COBRApy 辅助外部 GEM/SBML 导入

目标：在已有 `external_refs` 和 SBML GPR 解析链路旁增加可选 COBRApy import probe，用于核对外部 GEM 的 reaction、metabolite、gene、GPR 和 objective 语义。

适用范围：

- 外部模型 artifact cache 中已有 SBML 或可读模型文件。
- 需要判断 `python-libsbml` 解析和 COBRApy 解析是否一致。
- 需要给人工复核提供 import diagnostics，而不是直接改写当前 Pichia GEM。

输出建议：

- `local_runs/external_model_gpr_inventory/<run>/cobrapy_import_probe/manifest.json`
- `local_runs/external_model_gpr_inventory/<run>/cobrapy_import_probe/model_import_summary.tsv`
- `local_runs/external_model_gpr_inventory/<run>/cobrapy_import_probe/report.md`

核心字段：

- `model_id`
- `artifact_path`
- `cobrapy_available`
- `import_status`
- `reaction_count`
- `metabolite_count`
- `gene_count`
- `gpr_count`
- `objective_reaction`
- `id_sanitization_warnings`
- `libsbml_comparison_status`
- `manual_review_required`

非目标：

- 不把 COBRApy model 直接作为当前 pcSecPichia 生产模型。
- 不把外部 GPR 原样写入当前模型。
- 不在页面加载时自动联网或自动 import 大模型。

## 主线 B：MEMOTE / COBRApy 外部 GEM QA

目标：为外部 GEM 建立离线质量门，帮助判断某个外部模型是否值得进入人工复核、GPR 映射或后续 production cache。

适用范围：

- `model_inventory` 中已确认 artifact 可下载或已人工放入本地 cache。
- 需要比较 iPichia、ecPichia、Kp.1.0、iAUKM、Yeast8/Yeast9 等来源质量。
- 需要把 QA 结论作为 evidence source metadata，而不是作为模型 truth。

输出建议：

- `local_runs/external_model_gpr_inventory/<run>/gem_qa/manifest.json`
- `local_runs/external_model_gpr_inventory/<run>/gem_qa/model_qa_summary.tsv`
- `local_runs/external_model_gpr_inventory/<run>/gem_qa/report.md`
- 若 MEMOTE 可用，可保存 MEMOTE 原始 HTML/JSON 到同一 run 目录。

核心字段：

- `model_id`
- `artifact_path`
- `qa_backend`
- `backend_available`
- `qa_status`
- `memote_score`
- `stoichiometric_consistency_status`
- `annotation_score`
- `gpr_coverage`
- `blocked_reaction_count`
- `dead_end_metabolite_count`
- `manual_review_reasons`

非目标：

- 不把 MEMOTE 分数直接变成 KO/OE 推荐依据。
- 不要求所有开发环境默认安装 MEMOTE。
- 不把 QA 产物提交进源码仓库。

## 主线 C：Shadow LP Cross-check Service

目标：把已有 Shadow LP validation 从内部探索能力提升为可被研发人员使用的验证入口，用于对已保存的 hLF / OPN / custom target 运行结果做一致性对照和质量门判断。

推荐产品形态：

- Streamlit 中提供“对当前/已保存结果运行求解一致性验证”的按钮。
- 默认不在每次 screen 后自动运行，避免拖慢全基因组 KO/OE 工作流。
- 结果写入 `local_runs/shadow_lp_cross_check/<run>/`，页面只读取 manifest/report。

输出建议：

- `cross_check_manifest.json`
- `cross_check_summary.tsv`
- `cross_check_report.md`
- 可选：`reference_vs_shadow_diff.json`

核心字段：

- `target_id`
- `screen_run_id`
- `reference_capacity`
- `shadow_capacity`
- `absolute_diff`
- `relative_diff`
- `within_tolerance`
- `constraint_layer`
- `backend`
- `solver_status`
- `warnings`

非目标：

- 不把 Shadow LP 立即切成默认 solver。
- 不在 service/UI 中重新实现约束构建。
- 不把 cross-check 失败自动解释成实验不可行；它只说明模型/求解一致性需要复核。

## 执行轮次

### Round 0：边界审计和文档校准

目标：

- 确认当前 COBRApy、Shadow LP、external GEM/GPR、Streamlit 页面和测试入口状态。
- 确认 dirty worktree 中与本主线无关的文件不混入。
- 确认 optional dependency 策略：缺 COBRApy/MEMOTE 时稳定 unavailable/skip。

允许改动：

- 只允许修改 active docs 和测试白名单。

验收：

```powershell
python -m pytest -q tests\test_docs_active_boundary.py
git diff --name-only -- Code Model Enzymedata Results requirements.txt python_pichia\pyproject.toml
```

### Round 1：COBRApy import probe 契约

目标：

- 在 `python_pichia` 中定义外部 GEM import probe 的输入、输出和 unavailable 行为。
- 支持从 artifact cache 中读取本地模型文件并生成 import diagnostics。
- 若 COBRApy 不可用，稳定返回 `backend_available=false`，不影响默认流程。

允许改动：

- `python_pichia/src/pcsec_pichia/external_refs/`
- `scripts/`
- focused tests

验收：

```powershell
python -m pytest -q python_pichia\tests\test_external_model_import_probe.py
python -m compileall -q python_pichia\src scripts
git diff --name-only -- Code Model Enzymedata Results requirements.txt python_pichia\pyproject.toml
```

### Round 2：GEM QA / MEMOTE 离线质量门

目标：

- 增加外部 GEM QA 报告生成能力。
- 在 MEMOTE 不可用时保留 COBRApy/basic QA 或 unavailable 结果。
- 把 QA 结果作为 external model metadata，不接入 recommendation tier。

允许改动：

- `python_pichia/src/pcsec_pichia/external_refs/`
- `scripts/`
- focused tests

验收：

```powershell
python -m pytest -q python_pichia\tests\test_external_model_gem_qa.py
python -m compileall -q python_pichia\src scripts
git diff --name-only -- Code Model Enzymedata Results requirements.txt python_pichia\pyproject.toml
```

### Round 3：Shadow LP cross-check 核心服务

目标：

- 在 `python_pichia` 中提供读取已保存结果并运行 Shadow LP 对照的核心函数。
- 在 `app/services` 中提供 facade，不在 service 中实现科学判断。
- 结果写入 `local_runs/`，manifest/report 可被页面读取。

允许改动：

- `python_pichia/src/pcsec_pichia/analysis/shadow_lp/`
- `app/services/`
- focused tests

验收：

```powershell
python -m pytest -q python_pichia\tests\test_shadow_lp_cross_check.py tests\test_pichia_shadow_cross_check_service_contract.py
python -m compileall -q app python_pichia\src scripts
git diff --name-only -- Code Model Enzymedata Results requirements.txt python_pichia\pyproject.toml
```

### Round 4：Streamlit 入口和用户工作流

目标：

- 在合适的 Streamlit 页面增加“运行/查看 Shadow LP cross-check”的入口。
- 页面只触发 service、读取 manifest/report、展示状态和 warning。
- 检查页面入口、导航注册、旧模块引用、`st.session_state` 和缓存键。

允许改动：

- `app/ui/`
- `app/services/`
- UI/service contract tests

验收：

```powershell
python -m pytest -q tests\test_pichia_shadow_cross_check_service_contract.py tests\test_pichia_homology_audit_ui_contract.py tests\test_pichia_secretion_service_contract.py
python -m compileall -q app python_pichia\src scripts
git diff --name-only -- Code Model Enzymedata Results requirements.txt python_pichia\pyproject.toml
```

### Round 5：端到端 smoke 和收束

目标：

- 运行一个小型外部 GEM import/QA smoke。
- 运行一个 hLF 或 OPN Shadow LP cross-check smoke。
- 更新 active docs 的当前状态，避免继续追加大段设计。

验收：

```powershell
python scripts\build_external_model_gpr_inventory.py --output-dir local_runs\external_model_gpr_inventory\cobrapy_qa_smoke
python scripts\cache_external_model_artifacts.py --inventory-dir local_runs\external_model_gpr_inventory\cobrapy_qa_smoke --output-dir local_runs\external_model_gpr_inventory\cobrapy_qa_smoke_artifacts
python -m pytest -q python_pichia\tests\test_external_model_import_probe.py python_pichia\tests\test_external_model_gem_qa.py python_pichia\tests\test_shadow_lp_cross_check.py
python -m pytest -q tests\test_pichia_shadow_cross_check_service_contract.py tests\test_pichia_secretion_service_contract.py
python -m compileall -q app python_pichia\src scripts
git diff --name-only -- Code Model Enzymedata Results requirements.txt python_pichia\pyproject.toml
```

## 结束条件

本主线完成时，应满足：

- COBRApy import probe 可用于外部 GEM artifact 的离线诊断。
- GEM QA 可在 MEMOTE 可用或不可用两种环境下稳定输出状态。
- Shadow LP cross-check 可从 Streamlit/service 触发或读取，并生成研发可读报告。
- 默认 KO/OE screen、recommendation、report 和核心仿真不依赖 COBRApy/MEMOTE。
- 保护目录和依赖声明没有非预期 diff。

当前结论：前三项主能力已落地；真实 MEMOTE 环境验收和 custom target preparation 作为后续条件保留，不据此宣称 COBRApy 已替代 full pcSec 默认求解逻辑。
