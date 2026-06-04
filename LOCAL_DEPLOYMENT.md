# Local deployment notes

Repository: https://github.com/77652189/pcSecYeastSpecies

Remote `main` HEAD recorded during deployment:

```text
ef8d719dc5c9bb21f3fb576a94519aee2ac946f0
```

This local copy was provided as a downloaded worktree, so there is no `.git`
directory and `git rev-parse HEAD` is not available locally.

## MATLAB checks

Run from the repository root:

```matlab
startup_pcsec_local('yeast')
local_verify_phase1
local_smoke_sce_glc
```

Or run the MATLAB checks from PowerShell after MATLAB is installed:

```powershell
.\run_matlab_checks.ps1
```

If MATLAB is installed but not in `PATH`, set `MATLAB_EXE` first:

```powershell
$env:MATLAB_EXE = "C:\Program Files\MATLAB\R2026a\bin\matlab.exe"
.\run_matlab_checks.ps1
```

MATLAB itself is not bundled with this repository. Install MATLAB R2020b or
later from MathWorks Downloads:

```text
https://www.mathworks.com/downloads/
```

MathWorks documents that downloading MATLAB requires signing in to a
MathWorks Account with an eligible license.

`local_verify_phase1` loads `Model/pcSecYeast.mat` and saves:

```text
local_runs/phase1/Fig1b_ModelComparisonPpa.png
```

`local_smoke_sce_glc` writes one glucose LP and runner script to:

```text
local_runs/SCE_GLC_smoke/
```

Current validation status:

```text
local_verify_phase1: passed
local_smoke_sce_glc: passed for LP generation; MATLAB unique(...,'rows') cell
compatibility warning fixed in LP writers
local_soplex_solver_smoke: passed
real SCE glucose LP solve: passed with Docker SoPlex 6.0.4 using generated
sub_1.sh parameters for mu=0.10; output contains `problem is solved [optimal]`
and an Objective value.
```

## SoPlex

The generated `sub_1.sh` uses `SOPLEX_BIN` when it is set, otherwise it runs
`soplex` from `PATH`.
It defaults to `SOPLEX_READMODE=0`, which matches the generated floating-point
LP files; set `SOPLEX_READMODE=1` only if you intentionally want rational input
reading.

Example:

```bash
export SOPLEX_BIN=soplex
bash sub_1.sh
```

## Preflight and setup helpers

From PowerShell:

```powershell
.\local_preflight.ps1
.\setup_wsl_soplex.ps1
```

If WSL Ubuntu is not installed, Docker Desktop can run SoPlex instead:

```powershell
.\run_soplex_docker.ps1 -TimeoutSeconds 300
```

`run_soplex_docker.ps1` verifies every `.lp.out` target listed in `sub_1.sh`
and fails if SoPlex does not report `problem is solved [optimal]` or if the
output lacks an Objective value.

To verify the Docker SoPlex solver independently of the biological model LP:

```powershell
.\local_soplex_solver_smoke.ps1
```

## Streamlit local web app

Install Python dependencies if needed:

```powershell
pip install -r requirements.txt
```

Start the local/LAN Streamlit app:

```powershell
.\run_streamlit.ps1
```

The default URL is:

```text
http://localhost:8501
```

For another computer on the same LAN, open:

```text
http://<this-workstation-ip>:8501
```

The local Docker image has been built as:

```text
pcsec-soplex:24.04
```

It contains Ubuntu 24.04 and SoPlex 6.0.4 from the Ubuntu package repository.

MATLAB toolboxes are expected under:

```text
C:\Users\63097\MATLAB\toolboxes\
```

The local toolbox folders are:

```text
C:\Users\63097\MATLAB\toolboxes\cobratoolbox
C:\Users\63097\MATLAB\toolboxes\RAVEN
```
