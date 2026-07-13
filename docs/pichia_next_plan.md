# pcSecPichia 下一步计划

状态：active  
最后更新：2026-07-13

## 工作原则

- 一次只推进一个可验证轮次；每轮完成实现、focused tests、review、commit/push 后再进入下一轮。
- Python 核心逻辑放在 `python_pichia`，`app/services` 只做 facade，`app/ui` 只做展示。
- 不修改 `Code/`、`Model/`、`Enzymedata/`、`Results/`，除非明确声明为科学资产变更。
- 新运行产物、cache、报告和验证证据默认写入 ignored `local_runs/`。
- 输出只作为模型内相对解释、候选排序和研发复核依据，不声明 mg/L 绝对产量或实验成功率。

## 当前完成状态

- BLAST/RBH 本地同源审计已经形成可用链路：homology cache、name audit、rule-transfer audit、Streamlit 只读展示和 KO/OE 解释字段透传。
- Shadow LP 已作为 reference constrained solve 的替代候选完成 compare/validation 路径；纯 COBRApy/optlang 不作为 full pcSec 默认后端。
- 在线外部数据库证据层 Round 1-9 已完成：schema/cache IO、query builder、受控 clients、name resolution、KO/OE gene function evidence、external GPR candidate evidence、service/Streamlit、screen/report/LLM fact pack 透传和端到端 smoke 均已落地。
- 外部 GEM/GPR 资源优先级审计与导入主线 Round A-G 已完成：external model inventory、artifact cache/manual-required manifest、SBML GPR parse、GPR source priority、当前 Pichia GEM 映射报告、service/Streamlit/report/fact pack 透传和本地 Round G smoke 均已落地。最近 smoke 产物在 ignored `local_runs/external_model_gpr_inventory/round_g_smoke_20260709/`，不作为稳定科学资产提交。
- COBRApy import probe、GEM basic QA、Shadow LP cross-check core/service/Streamlit 已落地；hLF cross-check 与 reference 高精度对齐。MEMOTE 仍需专用环境真实验收，Shadow LP custom target preparation 尚未实现。
- 当前剩余边界是外部真实 GEM/MEMOTE 验收、custom target preparation、production cache 人工复核，以及是否提升稳定科学资产。

## 已完成主线：COBRApy 导入、GEM QA 与 Shadow LP 产品化

设计规格：[COBRApy 导入、GEM QA 与 Shadow LP 产品化计划](pichia_cobrapy_import_qa_shadow_plan.md)

目标：不把 COBRApy / optlang 切成 full pcSec 默认后端，而是把 COBRApy 用于外部 GEM/SBML import probe，把 MEMOTE/COBRApy 用于外部 GEM 离线 QA，并把 Shadow LP cross-check 做成 Streamlit/service 可触发的验证能力。

状态：核心 Round 0-5 已实现并完成 focused tests 与 hLF smoke。后续不重复实现现有模块，只补真实 MEMOTE 环境验收和 custom target preparation。

## 已完成主线：在线外部数据库证据层

设计规格：[在线外部数据库证据层架构](pichia_online_external_reference_architecture.md)

目标：受控联网获取 UniProt / NCBI / SGD 以及外部 GEM/SBML/API 证据，用于命名标准化、KO/OE gene function 补充和 external GPR candidate 评估；页面加载、核心仿真和 KO/OE screen 不依赖实时网络。KEGG 保留为后续可选扩展来源。

状态：已按 Round 1-9 完成并通过 focused tests、compileall、保护目录检查和小批量联网 smoke。当前 smoke 产物保留在 ignored `local_runs/`，不作为稳定科学资产提交。

### Dependency bundle

这条主线允许在 Round 1 一次性新增完整依赖，避免后续轮次反复修改依赖文件：

```text
httpx
tenacity
biopython
python-libsbml
pydantic
```

用途：

- `httpx`：统一官方数据库 HTTP client、timeout、headers、sync/async 扩展。
- `tenacity`：retry/backoff，处理数据库限流、瞬时网络失败和 5xx。
- `biopython`：NCBI Entrez、GenBank/FASTA/sequence record 解析。
- `python-libsbml`：外部 SBML/GEM 的 reaction、geneProduct、GPR 解析。
- `pydantic`：cache manifest、service payload、外部记录 schema validation；root `requirements.txt` 已有，但 `python_pichia/pyproject.toml` 需要同步声明。

不纳入本主线：

- `highspy`：属于 Shadow LP backend，不属于 external evidence 主线。
- `memote`：属于 GEM QA，可后续离线质量审计。
- `pyhmmer`：属于远缘同源/结构域增强，等 BLAST/RBH 和外部命名稳定后再评估。
- `bioservices`：多数据库封装依赖和许可边界复杂，先不用。
- `requests-cache`：本项目使用显式 JSONL/manifest cache，不用隐式 HTTP cache。

### Round 0：文档冻结和现状审计

改动范围：

- 只允许整理 active docs 和测试白名单。
- 确认 dirty worktree 中 unrelated AI report/env 改动不混入本主线。
- 确认当前 homology cache、gene catalog、KO/OE 输出字段和 Streamlit 页面入口。

验收：

```powershell
python -m pytest -q tests\test_docs_active_boundary.py
git diff --name-only -- Code Model Enzymedata Results requirements.txt python_pichia\pyproject.toml
```

### Round 1：schema 和 cache IO

改动范围：

- 允许修改 `requirements.txt` 和 `python_pichia/pyproject.toml`，一次性加入 dependency bundle 中缺失的依赖；本主线后续轮次原则上不再追加新库。
- 新增 `python_pichia/src/pcsec_pichia/external_refs/schema.py`。
- 新增 `python_pichia/src/pcsec_pichia/external_refs/cache_io.py`。
- 定义 external reference、gene function、reaction association、GPR candidate、manifest 等 dataclass。
- 实现 JSONL/manifest 读写、schema validation、重复键和 provenance 检查。
- 不联网，不接 UI，不改 KO/OE 输出。

验收：

```powershell
python -m pytest -q python_pichia\tests\test_external_refs_schema.py
python -m compileall -q python_pichia\src
git diff --name-only -- Code Model Enzymedata Results
```

### Round 2：query builder

改动范围：

- 新增 `external_refs/queries.py`。
- 从 homology cache、name audit、gene catalog、KO/OE candidate rows 生成 `ExternalReferenceQuery`。
- 区分 Pichia gene、SCE homolog、model gene、external accession 四类查询来源。
- 实现稳定去重、query fingerprint、source_context 和 warning。
- 不联网。

验收：

```powershell
python -m pytest -q python_pichia\tests\test_external_refs_schema.py python_pichia\tests\test_external_refs_queries.py
python -m compileall -q python_pichia\src
git diff --name-only -- Code Model Enzymedata Results requirements.txt python_pichia\pyproject.toml
```

### Round 3：受控在线 clients

改动范围：

- 新增 `external_refs/clients.py`、`uniprot.py`、`sgd.py`，可选 `ncbi.py`。
- 实现 timeout、retry、rate limit、user agent、API key env、source URL、retrieved_at、raw hash、失败记录。
- 默认单元测试使用 mock/fake HTTP，不依赖实时网络。
- 新增脚本 `scripts/build_pichia_external_reference_cache.py` 或改造已有脚本到新 schema。

验收：

```powershell
python -m pytest -q python_pichia\tests\test_external_refs_clients.py python_pichia\tests\test_external_refs_schema.py
python -m compileall -q python_pichia\src scripts
git diff --name-only -- Code Model Enzymedata Results requirements.txt python_pichia\pyproject.toml
```

联网 smoke 单独运行：

```powershell
python scripts\build_pichia_external_reference_cache.py --sources uniprot,sgd --limit 10 --output-dir local_runs\pichia_external_reference_cache\smoke
```

### Round 4：name resolution 和外部命名合并

改动范围：

- 新增 `external_refs/name_resolution.py` 和 `merge.py`。
- 把 external reference 合并到 name audit，输出 `external_match_confirmed`、`external_alias_confirmed`、`external_locus_confirmed`、`external_conflict`、`external_reference_missing`。
- 不覆盖内部 `gene_id`、RBH 事实、model gene_index 事实。
- 不改变 KO/OE `recommendation_tier`。

验收：

```powershell
python -m pytest -q python_pichia\tests\test_external_refs_merge.py python_pichia\tests\test_external_refs_name_resolution.py
python -m pytest -q python_pichia\tests\test_homology_crosswalk.py
python -m compileall -q python_pichia\src scripts
git diff --name-only -- Code Model Enzymedata Results requirements.txt python_pichia\pyproject.toml
```

### Round 5：KO/OE gene function evidence

改动范围：

- 新增 `external_refs/gene_function.py`。
- 从 external records 中抽取 protein name、function、EC、GO、pathway、orthology、reviewed status。
- 生成 `ExternalGeneFunctionEvidence` 和 `KoOeExternalGeneEvidence`。
- 接入 KO/OE candidate explanation fields，但不改变 recommendation tier。

验收：

```powershell
python -m pytest -q python_pichia\tests\test_external_refs_gene_function.py
python -m pytest -q python_pichia\tests\test_screens_entrypoints.py python_pichia\tests\test_yield_improvement_recommendations.py
python -m compileall -q python_pichia\src
git diff --name-only -- Code Model Enzymedata Results requirements.txt python_pichia\pyproject.toml
```

### Round 6：external GPR candidate evidence

改动范围：

- 新增 `external_refs/gpr_sources.py` 和 `external_refs/gpr_candidates.py`。
- 从外部 GEM/SBML/API 或缓存中读取 reaction association 和 gene rule。
- 将外部 reaction/gene rule 尝试映射到当前 Pichia GEM reaction/gene。
- 区分 `external_gpr_candidate`、`model_gpr_confirmed`、`reaction_mapping_required`、`gene_mapping_required`、`conflicting_gpr_sources`。
- 不把外部 GPR 原样写入当前模型。

验收：

```powershell
python -m pytest -q python_pichia\tests\test_external_refs_gpr_candidates.py
python -m pytest -q python_pichia\tests\test_screens_entrypoints.py
python -m compileall -q python_pichia\src scripts
git diff --name-only -- Code Model Enzymedata Results requirements.txt python_pichia\pyproject.toml
```

### Round 7：service 和 Streamlit

改动范围：

- 新增 `app/services/pichia_external_reference_service.py`。
- Streamlit 同源审计页面显示 external cache 状态、source counts、retrieved_at、name conflicts、gene function、GPR candidate、manual review reasons。
- 页面加载、筛选、导出不自动联网。
- 可提供手动刷新按钮或生成明确脚本命令。
- 更新 service/UI contract tests。

验收：

```powershell
python -m pytest -q tests\test_pichia_external_reference_service_contract.py tests\test_pichia_homology_audit_service_contract.py tests\test_pichia_homology_audit_ui_contract.py
python -m compileall -q app python_pichia\src scripts
git diff --name-only -- Code Model Enzymedata Results requirements.txt python_pichia\pyproject.toml
```

### Round 8：KO/OE screen、report 和 LLM fact pack 集成

改动范围：

- KO/OE preview rows、screen rows、yield recommendation rows、report fact pack 透传 external evidence fields。
- 历史结果缺字段时保持兼容。
- LLM summary 只能读取 fact pack，不直接读取任意运行目录；Judge 继续校验原始输入和总结一致性。
- 不把 external evidence 升级为 `experiment_calibrated`。

验收：

```powershell
python -m pytest -q python_pichia\tests\test_screens_entrypoints.py python_pichia\tests\test_yield_improvement_recommendations.py
python -m pytest -q tests\test_pichia_secretion_service_contract.py tests\test_screen_report_llm_contract.py
python -m compileall -q app python_pichia\src scripts
git diff --name-only -- Code Model Enzymedata Results requirements.txt python_pichia\pyproject.toml
```

### Round 9：端到端验收和收束

改动范围：

- 运行小批量联网 smoke，生成 `local_runs/pichia_external_reference_cache/<run_name>/`。
- 将 external cache 合并到 homology/name/rule-transfer audit。
- 在 Streamlit 中确认可查看、筛选、导出。
- 更新当前文档的“已完成/剩余边界”，不再新增大段设计。

验收：

```powershell
python scripts\build_pichia_external_reference_cache.py --sources uniprot,sgd --limit 10 --output-dir local_runs\pichia_external_reference_cache\smoke
python -m pytest -q python_pichia\tests\test_external_refs_schema.py python_pichia\tests\test_external_refs_clients.py python_pichia\tests\test_external_refs_merge.py python_pichia\tests\test_external_refs_name_resolution.py python_pichia\tests\test_external_refs_gene_function.py python_pichia\tests\test_external_refs_gpr_candidates.py
python -m pytest -q tests\test_pichia_external_reference_service_contract.py tests\test_pichia_homology_audit_service_contract.py tests\test_pichia_homology_audit_ui_contract.py tests\test_screen_report_llm_contract.py
python -m compileall -q app python_pichia\src scripts
git diff --name-only -- Code Model Enzymedata Results requirements.txt python_pichia\pyproject.toml
```

## 并行但非当前主线

### External GEM/GPR cache curation

当前外部 GEM/GPR 主线仍只把 iPichia、ecPichia、Kp.1.0、iAUKM、Yeast8/Yeast9、BioModels 和 GPRuler 作为可审计 evidence source。下一步如果要进入生产数据，应先人工复核 `local_runs/` 中的 inventory、manual download 记录、SBML parse、priority 和 mapping report，再决定是否提升为稳定 cache；不得把 external GPR candidate 当作当前 Pichia GEM 的 `model_gpr_executable`。

### Shadow LP service/backend toggle

把 `solve_shadow_secretion_capacity(...)` 接入 service / pipeline 的 opt-in 后端选择，让 reference 和 shadow 可双轨对比。默认仍为 reference，不立刻切换生产默认路径。

### KO/OE evidence summary consolidation

把 homology evidence、external evidence、model GPR executability、OE reaction proxy、phenotype evidence 合成同一张解释表，但所有核心判断仍放在 `python_pichia`，service/UI 只透传。

## 暂不做

- 不把 BLAST/RBH 命中自动写入 curated catalog。
- 不把同源命中直接当作 Pichia 模型 gene_id。
- 不把外部 GPR 原样写入当前 Pichia GEM。
- 不让页面加载、核心仿真或 KO/OE screen 自动联网；联网只允许在明确 builder/manual refresh 中发生。
- 不把外部 annotation、EC、GO、pathway 或 GPR candidate 当作 phenotype evidence。
- 不把 COBRApy / optlang 设为 full pcSec 默认后端。
- 不迁移历史 `Results/`、不改 Git LFS 或仓库历史。
