# pcSecPichia Python 独立包迁移进度

日期：2026-06-23

本文档只记录 `python_pichia/` 独立包的迁移状态。原始 MATLAB 项目仍作为 baseline 和论文复现参考完整保留。

## 当前边界

已建立独立 Python 项目目录：

```text
python_pichia/
  src/pcsec_pichia/
  tests/
  docs/
  pyproject.toml
```

第一批迁移代码来自当前已验证的 Python pcSecPichia 迁移模块，已放入：

- `python_pichia/src/pcsec_pichia/core/`
- `python_pichia/src/pcsec_pichia/adapters/`
- `python_pichia/src/pcsec_pichia/engines/`

当前仍保留旧 `app/` 迁移代码，以避免一次性破坏 Streamlit 和既有测试。后续应逐步让 `app/` 只作为 UI/service 调用层，核心 pcSecPichia 计算能力迁入 `python_pichia/`。

## 已迁入能力

第一批迁入：

- `ProjectPaths` 仓库路径定位；
- pcSecPichia `.mat` 模型读取；
- enzymedata 读取和合并数据结构；
- target / leader 通用输入 schema；
- target protein build plan；
- amino-acid stoichiometry；
- target reaction extension helpers；
- LP writer / parser / alignment；
- SoPlex output parser 和 Docker SoPlex solver adapter；
- `PythonPichiaEngine` 当前已迁移 secretion route 的 LP package / smoke 入口；
- `PichiaTargetSimulationService` 作为独立包 service 边界，负责把稳定的 route 名称分发到 engine 方法，并支持 OPN smoke 所需的可选 `production_ratio` 参数；
- `TargetRouteDefinition` registry 已开始承接 route 元数据、validation role 和 pipeline step 顺序，后续可替换为真正的 `SecretionPipeline` 执行器。

## 本轮对齐证据

本轮在 `local_runs/` 下补齐了 DSB+GPI synthetic extracellular route 的 MATLAB baseline 和 Python smoke 证据。

MATLAB baseline：

```text
local_runs/pichia_python/dsb_gpi_secretory_baseline/generate_dsb_gpi_matlab_baseline.m
local_runs/pichia_python/dsb_gpi_secretory_baseline/Simulation_dilutionDSB_GPI80_TEST_LEADER_mu0p10_media4_dsb_gpi_secretory_noMisfoldEq_noRiboEq_PP.lp
local_runs/pichia_python/dsb_gpi_secretory_baseline/Simulation_dilutionDSB_GPI80_TEST_LEADER_mu0p10_media4_dsb_gpi_secretory_noMisfoldEq_noRiboEq_PP.lp.float.out
local_runs/pichia_python/dsb_gpi_secretory_baseline/Simulation_dilutionDSB_GPI80_TEST_LEADER_mu0p10_media4_dsb_gpi_secretory_noMisfoldEq_noRiboEq_PP.lp.float_summary.json
```

Python output：

```text
local_runs/pichia_python/dsb_gpi_secretory_smoke/python_DSB_GPI80_TEST_LEADER_mu0p10_dsb_gpi_secretory_target_reference_constraints.lp
local_runs/pichia_python/dsb_gpi_secretory_smoke/python_DSB_GPI80_TEST_LEADER_mu0p10_dsb_gpi_secretory_target_reference_constraints_summary.json
local_runs/pichia_python/dsb_gpi_secretory_smoke/python_DSB_GPI80_TEST_LEADER_mu0p10_dsb_gpi_secretory_target_reference_constraints.lp.float.out
local_runs/pichia_python/dsb_gpi_secretory_smoke/python_DSB_GPI80_TEST_LEADER_mu0p10_dsb_gpi_secretory_target_reference_constraints.lp.float_summary.json
local_runs/pichia_python/dsb_gpi_secretory_smoke/dsb_gpi_matlab_python_objective_comparison.json
```

对齐结果：

```text
focused DSB+GPI tests: 8 passed
MATLAB SoPlex status: problem is solved [optimal]
Python SoPlex status: problem is solved [optimal]
MATLAB objective: 3.00021479e-03
Python objective: 3.00021479e-03
objective difference: 0.0
```

这只证明 synthetic `DSB>0, GPI>0, NG=0, OG=0, Transmembrane=0, localization=e` 小路线在当前 smoke 条件下已对齐，不代表完整 `addTargetProtein`、hLF、KO/OE 或真实产量预测已经完成。

## 未触碰原始目录

本轮修改应继续满足：

```powershell
git diff --name-only -- Code Model Enzymedata Results
```

输出为空。

## 下一步

1. 将阶段 7 验收口径收敛为 OPN 真实 target smoke 对齐 + hLF ready/draft 参数闸门。
2. 继续把 route 分发从逐方法硬编码迁到 `TargetRouteDefinition` / 后续 `SecretionPipeline`。
3. 只给已有 synthetic route 补缺失 alignment/objective 证据，不把理论 PTM 笛卡尔积设为阶段 8 前置条件。
4. 逐步把旧 `tests/` 中与 pcSecPichia engine 相关的测试迁移到 `python_pichia/tests/`。
5. 继续让 Streamlit / `app/services` 通过 `pcsec_pichia` service 调用独立包能力，而不是继续扩展 `app/` 内部 engine。
6. 后续新增 OPN、hLF 或其它目标蛋白时，必须通过通用 target schema/provider 接入。

## 阶段 7：DSB+GPI+OG route

本轮在独立包中新增 synthetic `DSB+GPI+OG extracellular` route：

```text
ThroughER=1
Cotranslation=0
DSB=2
GPI=1
OG=2
NG=0
Transmembrane=0
localization=e
```

新增包内能力：
- `TargetDsbGpiOgSecretoryBranchExtensionResult`；
- `add_dsb_gpi_og_secretory_branch_reactions()`；
- `add_target_dsb_gpi_og_and_misfolding_reactions()`；
- `add_target_dsb_gpi_og_er_golgi_maturation_transport_reactions()`；
- `PythonPichiaEngine.build_current_supported_dsb_gpi_og_secretory_target_model()`；
- `PythonPichiaEngine.write_current_supported_dsb_gpi_og_secretory_target_lp_alignment_package()`；
- `PythonPichiaEngine.run_current_supported_dsb_gpi_og_secretory_target_glucose_smoke()`；
- `PichiaTargetSimulationService.write_alignment_package("dsb_gpi_og_secretory", ...)`。

证据文件：

```text
local_runs/pichia_python/dsb_gpi_og_secretory_baseline/generate_dsb_gpi_og_matlab_baseline.m
local_runs/pichia_python/dsb_gpi_og_secretory_baseline/Simulation_dilutionDSB_GPI_OG80_TEST_LEADER_mu0p10_media4_dsb_gpi_og_secretory_noMisfoldEq_noRiboEq_PP.lp
local_runs/pichia_python/dsb_gpi_og_secretory_baseline/Simulation_dilutionDSB_GPI_OG80_TEST_LEADER_mu0p10_media4_dsb_gpi_og_secretory_noMisfoldEq_noRiboEq_PP.lp.float.out
local_runs/pichia_python/dsb_gpi_og_secretory_baseline/Simulation_dilutionDSB_GPI_OG80_TEST_LEADER_mu0p10_media4_dsb_gpi_og_secretory_noMisfoldEq_noRiboEq_PP.lp.float_summary.json
local_runs/pichia_python/dsb_gpi_og_secretory_smoke/python_DSB_GPI_OG80_TEST_LEADER_mu0p10_dsb_gpi_og_secretory_target_reference_constraints.lp
local_runs/pichia_python/dsb_gpi_og_secretory_smoke/python_DSB_GPI_OG80_TEST_LEADER_mu0p10_dsb_gpi_og_secretory_target_reference_constraints_summary.json
local_runs/pichia_python/dsb_gpi_og_secretory_smoke/python_DSB_GPI_OG80_TEST_LEADER_mu0p10_dsb_gpi_og_secretory_target_reference_constraints.lp.float.out
local_runs/pichia_python/dsb_gpi_og_secretory_smoke/python_DSB_GPI_OG80_TEST_LEADER_mu0p10_dsb_gpi_og_secretory_target_reference_constraints.lp.float_summary.json
local_runs/pichia_python/dsb_gpi_og_secretory_smoke/dsb_gpi_og_matlab_python_lp_math_alignment.json
local_runs/pichia_python/dsb_gpi_og_secretory_smoke/dsb_gpi_og_matlab_python_soplex_objective_comparison.json
```

对齐结果：

```text
focused DSB+GPI+OG tests: 8 passed
MATLAB LP math alignment: match, constraint differences 0, bound differences 0
MATLAB SoPlex status: problem is solved [optimal]
Python SoPlex status: problem is solved [optimal]
MATLAB objective: 2.96485685e-03
Python objective: 2.96485685e-03
objective difference: 0.0
```

这只证明 synthetic `DSB>0, GPI>0, OG>0, NG=0, Transmembrane=0, localization=e` 小路线在当前 smoke 条件下与 MATLAB baseline 对齐；不代表真实产量预测、hLF、KO/OE、瓶颈解释或其它 Stage 7 route 已完成。

## 阶段 7：DSB+GPI+NG route 状态修正

当前代码已经包含 synthetic `DSB+GPI+NG extracellular` route：

```text
ThroughER=1
Cotranslation=0
DSB=2
GPI=1
NG=2
OG=0
Transmembrane=0
localization=e
```

已有包内能力：
- `TargetDsbGpiNgSecretoryBranchExtensionResult`；
- `add_dsb_gpi_ng_secretory_branch_reactions()`；
- `PythonPichiaEngine.write_current_supported_dsb_gpi_ng_secretory_target_lp_alignment_package()`；
- `PythonPichiaEngine.run_current_supported_dsb_gpi_ng_secretory_target_glucose_smoke()`；
- `PichiaTargetSimulationService.write_alignment_package("dsb_gpi_ng_secretory", ...)`；
- `TargetRouteDefinition(key="dsb_gpi_ng_secretory", validation_role="synthetic_regression", ...)`。

已有本地产物：

```text
local_runs/pichia_python/dsb_gpi_ng_secretory_baseline/
local_runs/pichia_python/dsb_gpi_ng_secretory_smoke/
```

尚缺证据：

```text
dsb_gpi_ng_matlab_python_soplex_objective_comparison.json
```

因此当前结论是“代码和 LP 产物已有，但完整 objective 对齐证据待补”，不能再写成“无独立 route”，也不能写成“完整对齐完成”。
