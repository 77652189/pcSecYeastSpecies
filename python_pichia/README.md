# pcSecPichia Python Engine

This directory contains the Python engine for the current Pichia productivity
workflow. The user-facing goal is to help biology R&D colleagues reason about
hLF and OPN production improvement candidates, especially KO/OE interventions
and relevant migrated MATLAB functionality.

The original MATLAB project remains reference material. Do not move, delete, or
rewrite these original directories:

- `../Code/`
- `../Model/`
- `../Enzymedata/`
- `../Results/`

Generated LP files, cache files, solver outputs, and local evidence artifacts
belong under `../local_runs/`.

## Scope

Current scope:

- load pcSecPichia model inputs from the reference MATLAB data;
- build OPN, hLF, and custom target protein plans;
- apply corrected medium/carbon-source conditions;
- run draft secretion simulations, growth tradeoff probes, and cost summaries;
- run small KO/OE candidate screens and evidence-aware recommendations;
- write structured summaries, candidate tables, and reports for app/service use.

Out of scope for this package right now:

- full three-species MATLAB workflow replacement;
- paper figure reproduction;
- automatic new MATLAB baseline generation;
- absolute mg/L production prediction;
- full-model KO/OE batch screening unless separately scoped.

## Layout

```text
python_pichia/
  src/pcsec_pichia/
    core/
    adapters/
    engines/
  tests/
  pyproject.toml
```

The package boundary is intentional:

- `loading`, `media`: reference model inputs and medium/carbon-source setup;
- `targets`, `secretion_plan`: target protein and secretory route planning;
- `constraints`, `simulation`: model constraints and solver-facing simulation;
- `screens`, `analysis`, `reports`: KO/OE screens, interpretation, and outputs;
- `services`: gene evidence, gene catalog, and rule overlay support;
- `tests`: package-local regression and focused validation tests.

## Validation

Use focused tests from the repository root. A common gate is:

```powershell
python -m compileall -q app python_pichia\src\pcsec_pichia
python -m pytest -q python_pichia\tests\test_screens_entrypoints.py python_pichia\tests\test_yield_improvement_recommendations.py
python -m pytest -q tests\test_pichia_secretion_service_contract.py
git diff --name-only -- Code Model Enzymedata Results
```

The last command must output nothing.

The active project-level docs are:

- `../docs/pichia_current_architecture_and_requirements.md`
- `../docs/pichia_next_plan.md`
