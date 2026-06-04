param(
    [string]$MatlabExe = $env:MATLAB_EXE,
    [switch]$SmokeOnly
)

$ErrorActionPreference = "Stop"

function Find-Matlab {
    param([string]$ExplicitPath)

    if ($ExplicitPath -and (Test-Path -LiteralPath $ExplicitPath)) {
        return (Resolve-Path -LiteralPath $ExplicitPath).Path
    }

    $cmd = Get-Command matlab -ErrorAction SilentlyContinue
    if ($cmd) {
        return $cmd.Source
    }

    $roots = @(
        "C:\Program Files\MATLAB",
        "C:\Program Files (x86)\MATLAB"
    )

    foreach ($root in $roots) {
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

    throw "MATLAB executable not found. Install MATLAB R2020b+ or set MATLAB_EXE to matlab.exe."
}

$repoRoot = $PSScriptRoot
$matlab = Find-Matlab $MatlabExe

if ($SmokeOnly) {
    $matlabCommand = "cd('$repoRoot'); startup_pcsec_local('yeast'); local_smoke_sce_glc;"
} else {
    $matlabCommand = "cd('$repoRoot'); startup_pcsec_local('yeast'); local_verify_phase1; local_smoke_sce_glc;"
}

Write-Host "Using MATLAB: $matlab"
& $matlab -batch $matlabCommand

