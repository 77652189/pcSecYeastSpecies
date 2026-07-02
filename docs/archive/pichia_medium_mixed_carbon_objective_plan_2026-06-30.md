# pcSecPichia mixed-carbon medium objective plan

## Current conclusion

`glucose_glycerol` and `glycerol_methanol` are active Python draft medium
boundaries built from model-native exchange reactions. They are not old MATLAB
baseline conditions.

The MATLAB harness can reproduce the same exchange bounds, but legacy
`writeLPGlc` still optimizes `Ex_glc_D`. That is useful for artifact comparison,
but it is not a complete mixed-carbon cost objective.

## Product behavior for now

- Keep glucose as the corrected default reference.
- Keep non-glucose and mixed-carbon conditions available as explicit optional
  draft boundaries.
- Surface `scientific_status` and medium warnings in summary/report/UI.
- Do not interpret glucose+glycerol as AOX1 methanol induction biology.
- Do not claim MATLAB fully aligned mixed-carbon results until a matching
  objective or uptake-budget formulation is defined.

## Recommended V1 objective design

Add a separate opt-in formulation, without changing the corrected default:

- `carbon_objective_mode = "single_exchange"`
  - current behavior; objective is one exchange such as `Ex_glc_D`.
- `carbon_objective_mode = "weighted_carbon_uptake"`
  - minimize a weighted sum of active carbon uptake reactions.
  - example variables: `Ex_glc_D`, `Ex_glyc`, `Ex_meoh`.
  - weights should be explicit and auditable, not implicit.
- `carbon_objective_mode = "fixed_uptake_budget"`
  - fix or cap active carbon uptake reactions and maximize target secretion.
  - suitable for comparing feed recipes or experimental designs.

## Implemented probe

`run_mixed_carbon_objective_probe(...)` now implements an opt-in
`weighted_carbon_uptake` formulation. It introduces non-negative uptake-cost
auxiliary variables so the objective minimizes uptake cost instead of rewarding
positive exchange flux.

Current validation artifact:

- `local_runs/mixed_carbon_objective_probe_2026-06-30/`

Observed behavior:

- Equal glucose/glycerol weights keep the optimum on glucose uptake.
- Mild glucose penalty also keeps the optimum on glucose uptake.
- A strong glucose penalty, tested as `Ex_glc_D=10.0` and `Ex_glyc=1.0`,
  shifts the optimum to glycerol uptake for both OPN and hLF.

Interpretation: the mixed-carbon objective is functional and can express carbon
preference assumptions, but the weights are scientific assumptions and need
experimental or literature support before product-level claims.

## Validation gates

1. Bound parity: Python and MATLAB harness use identical carbon exchange and
   biomass-variant bounds.
2. Objective parity: Python and MATLAB use the same objective semantics.
3. Target parity: OPN and hLF project target parameters match the current
   Python target definitions.
4. LP diff: objective, bounds, row labels, RHS, sparsity, and coefficients are
   compared before numerical claims.
5. Scientific annotation: reports distinguish growth boundary, induction
   boundary, mixed-feed boundary, and legacy MATLAB artifact compatibility.

## Open scientific choices

- Whether glucose+glycerol should have fixed uptake ratios.
- Whether glycerol+methanol should model an induction phase, a transition phase,
  or a steady mixed-feed phase.
- Whether AOX1 promoter repression/induction should be represented as a
  separate regulatory constraint instead of a medium condition.
