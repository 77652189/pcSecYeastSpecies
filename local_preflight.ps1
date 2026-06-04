param(
    [string]$RepoRoot = $PSScriptRoot
)

$ErrorActionPreference = "Stop"

function Test-Command {
    param([string]$Name)
    $cmd = Get-Command $Name -ErrorAction SilentlyContinue
    if ($cmd) {
        [PSCustomObject]@{
            Name = $Name
            Status = "OK"
            Detail = $cmd.Source
        }
    } else {
        [PSCustomObject]@{
            Name = $Name
            Status = "Missing"
            Detail = ""
        }
    }
}

function Test-File {
    param([string]$Path)
    $fullPath = Join-Path $RepoRoot $Path
    [PSCustomObject]@{
        Name = $Path
        Status = if (Test-Path -LiteralPath $fullPath) { "OK" } else { "Missing" }
        Detail = $fullPath
    }
}

function Test-Directory {
    param([string]$Path)
    [PSCustomObject]@{
        Name = $Path
        Status = if (Test-Path -LiteralPath $Path -PathType Container) { "OK" } else { "Missing" }
        Detail = $Path
    }
}

function Find-MatlabExecutable {
    $cmd = Get-Command matlab -ErrorAction SilentlyContinue
    if ($cmd) {
        return $cmd.Source
    }

    foreach ($root in @("C:\Program Files\MATLAB", "C:\Program Files (x86)\MATLAB")) {
        if (-not (Test-Path -LiteralPath $root)) {
            continue
        }

        $candidate = Get-ChildItem -LiteralPath $root -Directory -ErrorAction SilentlyContinue |
            Sort-Object Name -Descending |
            ForEach-Object { Join-Path $_.FullName "bin\matlab.exe" } |
            Where-Object { Test-Path -LiteralPath $_ } |
            Select-Object -First 1

        if ($candidate) {
            return $candidate
        }
    }

    return $null
}

Write-Host "== Commands =="
Test-Command git
$matlabExe = Find-MatlabExecutable
if ($matlabExe) {
    [PSCustomObject]@{
        Name = "matlab"
        Status = "OK"
        Detail = $matlabExe
    }
} else {
    [PSCustomObject]@{
        Name = "matlab"
        Status = "Missing"
        Detail = ""
    }
}
Test-Command wsl
$windowsSoplex = Test-Command soplex
if ($windowsSoplex.Status -eq "Missing") {
    [PSCustomObject]@{
        Name = "soplex (Windows PATH)"
        Status = "Optional"
        Detail = "Not installed; Docker image is used for local SoPlex runs"
    }
} else {
    [PSCustomObject]@{
        Name = "soplex (Windows PATH)"
        Status = "OK"
        Detail = $windowsSoplex.Detail
    }
}
Test-Command docker

Write-Host "`n== Required repository files =="
Test-File "README.md"
Test-File "startup_pcsec_local.m"
Test-File "local_verify_phase1.m"
Test-File "local_smoke_sce_glc.m"
Test-File "run_matlab_checks.ps1"
Test-File "run_soplex_docker.ps1"
Test-File "local_soplex_solver_smoke.ps1"
Test-File "docker/soplex/Dockerfile"
Test-File "Model/pcSecYeast.mat"
Test-File "Enzymedata/pcSecYeast/enzymedata_SCE.mat"
Test-File "Code/Figures/Fig1b_ModelComparisonPpa.m"
Test-File "Results/CSource/CSource_res.xlsx"

Write-Host "`n== WSL distributions =="
if (Get-Command wsl -ErrorAction SilentlyContinue) {
    & wsl.exe -l -v
} else {
    Write-Host "wsl.exe not found."
}

Write-Host "`n== MATLAB toolbox directory =="
$toolboxRoot = Join-Path $env:USERPROFILE "MATLAB/toolboxes"
if (Test-Path -LiteralPath $toolboxRoot) {
    Get-ChildItem -LiteralPath $toolboxRoot | Select-Object Name,FullName
} else {
    Write-Host "Missing: $toolboxRoot"
}

Write-Host "`n== Required MATLAB toolboxes =="
Test-Directory (Join-Path $toolboxRoot "cobratoolbox")
Test-Directory (Join-Path $toolboxRoot "RAVEN")
Test-File "startup_pcsec_local.m"

Write-Host "`n== Docker SoPlex image =="
if (Get-Command docker -ErrorAction SilentlyContinue) {
    $image = docker images pcsec-soplex --format "{{.Repository}}:{{.Tag}} {{.Size}}"
    if ($image) {
        [PSCustomObject]@{
            Name = "pcsec-soplex"
            Status = "OK"
            Detail = $image
        }
    } else {
        [PSCustomObject]@{
            Name = "pcsec-soplex"
            Status = "Missing"
            Detail = "Run: docker build -t pcsec-soplex:24.04 docker/soplex"
        }
    }
} else {
    [PSCustomObject]@{
        Name = "pcsec-soplex"
        Status = "Missing"
        Detail = "docker.exe not found"
    }
}
