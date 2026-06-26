# pcSecPichia Python 验证命令

日期：2026-06-26
状态：active validation checklist

## 每轮固定检查

```powershell
git diff --name-only -- Code Model Enzymedata Results
```

期望输出为空。默认不启动 MATLAB，不生成新 MATLAB baseline，不回退 Python corrected 默认条件。

## 日常快测

用于 app/service/UI/docs 等小改动。

```powershell
python -m py_compile app\__init__.py app\services\pichia_secretion_service.py app\ui\views\simulation.py
python -m pytest -q tests\test_pichia_secretion_service_contract.py
git diff --name-only -- Code Model Enzymedata Results
```

## Engine 聚焦测试

只跑受影响模块，不跑全量。

```powershell
python -m pytest -q python_pichia\tests\test_target_entrypoints.py
python -m pytest -q python_pichia\tests\test_secretion_plan_entrypoints.py
python -m pytest -q python_pichia\tests\test_constraints_entrypoints.py
python -m pytest -q python_pichia\tests\test_simulation_entrypoints.py
python -m pytest -q python_pichia\tests\test_reports_entrypoints.py
python -m pytest -q python_pichia\tests\test_alignment_entrypoints.py
python -m pytest -q python_pichia\tests\test_pipeline_entrypoints.py
```

`test_pipeline_entrypoints.py` 默认只跑快速 contract；真实求解慢测需要显式打开：

```powershell
$env:PCSEC_RUN_SLOW_PIPELINE_TESTS="1"
python -m pytest -q python_pichia\tests\test_pipeline_entrypoints.py
Remove-Item Env:\PCSEC_RUN_SLOW_PIPELINE_TESTS
```

`test_screens_entrypoints.py` 里的真实 KO/OE 求解 smoke 也默认跳过，需要显式打开：

```powershell
$env:PCSEC_RUN_SLOW_SCREEN_TESTS="1"
python -m pytest -q python_pichia\tests\test_screens_entrypoints.py
Remove-Item Env:\PCSEC_RUN_SLOW_SCREEN_TESTS
```

## UI/API 聚焦测试

```powershell
python -m py_compile app\services\pichia_secretion_service.py app\ui\views\simulation.py app\api\pichia_secretion_api.py
python -m pytest -q tests\test_pichia_secretion_service_contract.py tests\test_pichia_fastapi_entrypoints.py
```

Streamlit 手动检查只在需要看交互时进行：

```powershell
powershell -ExecutionPolicy Bypass -File scripts\start_pcSecYeastSpecies_lan.ps1
```

当前项目使用 8502；8501 可能属于其它项目。

## Release 前验证

```powershell
python -m compileall app python_pichia tests
python -m pytest -q tests\test_pichia_secretion_service_contract.py tests\test_pichia_fastapi_entrypoints.py
python -m pytest -q python_pichia\tests\test_target_entrypoints.py python_pichia\tests\test_secretion_plan_entrypoints.py python_pichia\tests\test_constraints_entrypoints.py
python -m pytest -q python_pichia\tests\test_simulation_entrypoints.py python_pichia\tests\test_reports_entrypoints.py python_pichia\tests\test_alignment_entrypoints.py python_pichia\tests\test_pipeline_entrypoints.py
git diff --name-only -- Code Model Enzymedata Results
```

如果要包含真实 pipeline 或 screen 求解，在 release 验证中额外运行 slow pipeline / slow screen 命令。

## 最近一次非 slow rehearsal

日期：2026-06-26
分支：`codex/pichia-python-draft-engine`
状态：通过 review-ready focused gate。

已运行：

```powershell
python -m compileall app python_pichia tests
python -m pytest -q tests\test_pichia_secretion_service_contract.py tests\test_review_package_boundaries.py tests\test_docs_active_boundary.py tests\test_streamlit_startup_scripts.py
python -m pytest -q python_pichia\tests\test_target_entrypoints.py python_pichia\tests\test_secretion_plan_entrypoints.py python_pichia\tests\test_alignment_entrypoints.py python_pichia\tests\test_pipeline_entrypoints.py
git diff --name-only -- Code Model Enzymedata Results
```

结果：

- `compileall` 通过。
- app/service/docs/startup 聚焦测试：35 passed。
- engine target/secretion/alignment/pipeline 聚焦测试：32 passed, 5 skipped。
- `Code/`、`Model/`、`Enzymedata/`、`Results/` diff 为空。
- 未启动 MATLAB，未打开 slow pipeline/screen/probe 环境变量。

旧 probe 迁移回归也属于慢任务，默认跳过；需要显式打开：

```powershell
$env:PCSEC_RUN_SLOW_PROBE_TESTS="1"
python -m pytest -q python_pichia\tests\test_probe_migration.py
Remove-Item Env:\PCSEC_RUN_SLOW_PROBE_TESTS
```

## 慢任务规则

以下任务不属于日常门禁：

- MATLAB harness / baseline 生成。
- 新 LP artifact 生成。
- 行级 LP diff。
- 全模型 KO/OE 批量筛选。
- 长网格 growth tradeoff。

需要运行这些任务时，先说明命令、输出目录、预计影响，并把产物写入 `local_runs/` 下的新子目录。
