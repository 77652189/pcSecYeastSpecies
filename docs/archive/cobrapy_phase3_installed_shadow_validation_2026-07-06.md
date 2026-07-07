# COBRApy Phase 3 Installed Shadow Validation

Date: 2026-07-06

Scope: `pcSecYeastSpecies / python_pichia`

## Summary

Phase 3 created an isolated local-only COBRApy environment and used it to run the existing real-model shadow FBA parity harness against `Model/pcSecPichia.mat`.

Result: `completed_within_tolerance`.

All configured objective cases solved with both the current SciPy/HiGHS path and the optional COBRApy shadow path. Objective values were within the configured tolerance for every case.

This result only validates base GEM FBA conversion parity for stoichiometry, bounds, and objective reactions. It is not a pcSec constraint migration, not a KO/OE evidence source, not a recommendation tier source, and not an mg/L or experimental success prediction.

## Isolated Environment

Environment path:

```text
local_runs/cobrapy_shadow_env/
```

This environment is under `local_runs/` and is ignored/local-only.

Installed in the isolated environment:

| package | status/version |
| --- | --- |
| cobra | 0.31.1 |
| optlang | 1.9.1 |
| libsbml | 5.21.1 |
| swiglpk | available |
| scipy | 1.18.0 |
| h5py | 3.16.0 |
| openpyxl | 3.1.5 |

Default environment check remained clean:

```text
cobra_default_installed=False
optlang_default_installed=False
swiglpk_default_installed=False
libsbml_default_installed=False
```

No COBRApy dependency was added to `requirements.txt` or `python_pichia/pyproject.toml`.

## Harness Command

The harness was run from the isolated venv with `PYTHONPATH` pointed at `python_pichia/src`:

```powershell
$env:PYTHONPATH = "python_pichia\src"
local_runs\cobrapy_shadow_env\Scripts\python.exe -m pcsec_pichia.analysis.cobrapy_shadow_baseline --output-dir local_runs\cobrapy_shadow_phase3
```

Output artifacts:

```text
local_runs/cobrapy_shadow_phase3/cobrapy_shadow_baseline_summary.json
local_runs/cobrapy_shadow_phase3/cobrapy_shadow_baseline_report.md
local_runs/cobrapy_shadow_phase3/console.log
```

These are local validation artifacts and should not be committed.

## Model Summary

The harness loaded the current real model from:

```text
Model/pcSecPichia.mat
```

Model summary:

| field | value |
| --- | ---: |
| model_id | splited_model |
| reactions | 29026 |
| metabolites | 20195 |
| genes | 1025 |
| stoichiometric shape | 20195 x 29026 |
| rxn-gene shape | 29026 x 1025 |

## Objective Parity Results

| case | current objective | COBRApy shadow objective | abs diff | rel diff | within tolerance |
| --- | ---: | ---: | ---: | ---: | --- |
| BIOMASS_maximize | 1.6152188242557852 | 1.6152188262259553 | 1.9701700271212985e-09 | 1.219754249475104e-09 | true |
| Ex_glc_D_maximize | -0.08781250000008711 | -0.08781249999999992 | 8.719414079649823e-14 | 9.929581870054005e-13 | true |
| Ex_glyc_maximize | 0.0 | 0.0 | 0.0 | 0.0 | true |
| Ex_meoh_maximize | 0.0 | 0.0 | 0.0 | 0.0 | true |
| Ex_o2_maximize | 0.0 | 0.0 | 0.0 | 0.0 | true |

Run status:

```text
completed_within_tolerance
```

## Key Flux Diff Summary

Objective values matched within tolerance for all cases. Some key fluxes differed even when objective values were equal:

| case | max key flux diff | nonzero key flux diffs |
| --- | ---: | --- |
| BIOMASS_maximize | 1.2355531708863055e-08 | `BIOMASS=1.9701700271212985e-09`, `Ex_o2=1.2355531708863055e-08` |
| Ex_glc_D_maximize | 2.9187763317395365e-13 | `Ex_glc_D=8.719414079649823e-14`, `Ex_o2=2.9187763317395365e-13` |
| Ex_glyc_maximize | 59.99999999999964 | `Ex_glc_D=8.594999999999848`, `Ex_o2=59.99999999999964` |
| Ex_meoh_maximize | 59.99999999999964 | `Ex_glc_D=8.594999999999848`, `Ex_o2=59.99999999999964` |
| Ex_o2_maximize | 8.59499999999993 | `Ex_glc_D=8.59499999999993` |

The large key flux differences in zero-objective exchange cases should be treated as likely non-unique optimum / flux degeneracy behavior. They are not evidence of different production capacity and must not be interpreted as target protein yield, mg/L output, or experiment success probability.

## Review Fixes During Phase 3

The installed COBRApy run exposed two conversion-boundary issues that were fixed in the optional shadow adapter:

- SciPy sparse arrays such as `csc_array` are normalized to `csc_matrix` before column access.
- COBRApy object ids are sanitized for whitespace while original metabolite/reaction names remain available for interpretation; reaction flux outputs are mapped back to original reaction ids.

These fixes remain inside `pcsec_pichia.adapters.cobrapy_shadow` and do not change default pipeline, screen, report, service, or UI behavior.

## Explicit Non-Claims

This Phase 3 validation does not:

- migrate pcSec protein/secretion constraints to COBRApy;
- replace current SciPy/HiGHS or MATLAB-aligned paths;
- validate KO/OE GPR planning semantics;
- validate OE reaction proxy semantics;
- provide phenotype evidence or recommendation tiers;
- affect Streamlit UI, FastAPI/service, pipeline, screen, or report defaults;
- predict mg/L, absolute production, fermentation outcome, or experiment success rate.

## Phase 4 Recommendation

Phase 4 is reasonable only as another opt-in developer validation slice.

Recommended Phase 4 scope:

- keep COBRApy optional and local/dev-only;
- compare additional base GEM objective/bound scenarios;
- add focused handling for solver-status conventions and flux degeneracy notes;
- do not connect COBRApy shadow results to KO/OE recommendation, phenotype tiers, report ranking, or user-facing production claims.

Do not start Phase 4 by replacing pcSec constraints or default solving behavior.
