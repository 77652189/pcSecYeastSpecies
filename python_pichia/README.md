# pcSecPichia Python Engine

This directory is the independent Python migration area for the Pichia pastoris `pcSecPichia` workflow.

The original MATLAB project remains the source of truth for baseline alignment and paper reproduction. Do not move, delete, or rewrite these original directories:

- `../Code/`
- `../Model/`
- `../Enzymedata/`
- `../Results/`

Python code in this package reads the original model/data as input, writes generated LP files and solver outputs to `../local_runs/`, and records alignment evidence against MATLAB baselines.

## Scope

Current migration scope:

- read `Model/pcSecPichia.mat`;
- read pcSecPichia enzymedata;
- build target-protein plans through a generic target/leader schema;
- generate LP files for migrated secretion routes;
- solve LP files through a replaceable solver adapter;
- compare Python outputs with MATLAB baseline LP and SoPlex outputs.

Not yet complete:

- full MATLAB workflow replacement;
- hLF production-ready target simulation;
- KO/OE screening;
- bottleneck explanation;
- expert-facing standalone UI.

## Layout

```text
python_pichia/
  src/pcsec_pichia/
    core/
    adapters/
    engines/
  tests/
  docs/
  pyproject.toml
```

The package boundary is intentional:

- `core`: domain data structures and model objects;
- `adapters`: file formats, MATLAB `.mat` loading, LP parsing/writing, SoPlex process wrappers;
- `engines`: pcSecPichia calculation orchestration;
- `tests`: package-local regression and alignment tests;
- `docs`: package-local migration notes.

## Validation

Run from the repository root:

```powershell
python -m compileall python_pichia app scripts tests
python -m pytest -q
git diff --name-only -- Code Model Enzymedata Results
```

The last command must output nothing.

Only routes with MATLAB baseline alignment should be marked as restored. Routes without baseline alignment are draft or partial migration only.
