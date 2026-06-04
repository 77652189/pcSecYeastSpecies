param(
    [string]$RunDir = "local_runs/SCE_GLC_smoke",
    [string]$Image = "pcsec-soplex:24.04",
    [int]$TimeoutSeconds = 0
)

$ErrorActionPreference = "Stop"

$repoRoot = $PSScriptRoot
$dockerfile = Join-Path $repoRoot "docker/soplex/Dockerfile"
$runPath = Join-Path $repoRoot $RunDir

if (-not (Get-Command docker -ErrorAction SilentlyContinue)) {
    throw "docker.exe was not found. Install/start Docker Desktop first."
}

if (-not (Test-Path -LiteralPath $dockerfile)) {
    throw "Missing Dockerfile: $dockerfile"
}

if (-not (Test-Path -LiteralPath $runPath -PathType Container)) {
    throw "Missing run directory: $runPath. Run local_smoke_sce_glc in MATLAB first."
}

if (-not (Test-Path -LiteralPath (Join-Path $runPath "sub_1.sh"))) {
    throw "Missing sub_1.sh in $runPath. Run local_smoke_sce_glc in MATLAB first."
}

$runner = Join-Path $runPath "sub_1.sh"
$runnerText = Get-Content -LiteralPath $runner -Raw
$outputFiles = [regex]::Matches($runnerText, '>\s*([^\s]+\.lp\.out)') |
    ForEach-Object { $_.Groups[1].Value } |
    Select-Object -Unique

if (-not $outputFiles) {
    throw "No .lp.out target found in $runner"
}

docker build -t $Image (Join-Path $repoRoot "docker/soplex")

$mountPath = (Resolve-Path -LiteralPath $runPath).Path
if ($TimeoutSeconds -gt 0) {
    docker run --rm -v "${mountPath}:/work" -w /work $Image sh -lc "timeout $TimeoutSeconds bash sub_1.sh"
} else {
    docker run --rm -v "${mountPath}:/work" -w /work $Image bash sub_1.sh
}

foreach ($outputFile in $outputFiles) {
    $outputPath = Join-Path $runPath $outputFile
    if (-not (Test-Path -LiteralPath $outputPath)) {
        throw "Missing SoPlex output: $outputPath"
    }

    $outputText = Get-Content -LiteralPath $outputPath -Raw
    if ($outputText -notmatch 'problem is solved \[optimal\]') {
        throw "SoPlex did not report optimal in $outputPath"
    }

    if ($outputText -notmatch 'Objective value') {
        throw "SoPlex output did not contain an objective value in $outputPath"
    }

    Write-Host "Verified SoPlex output: $outputPath"
}
