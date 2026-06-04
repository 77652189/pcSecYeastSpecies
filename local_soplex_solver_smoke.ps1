param(
    [string]$Image = "pcsec-soplex:24.04"
)

$ErrorActionPreference = "Stop"

$repoRoot = $PSScriptRoot
$runPath = Join-Path $repoRoot "local_runs/soplex_solver_smoke"
New-Item -ItemType Directory -Force -Path $runPath | Out-Null

$lpPath = Join-Path $runPath "solver_smoke.lp"
$outPath = Join-Path $runPath "solver_smoke.lp.out"

@"
Maximize
obj: X1
Subject To
C1: X1 <= 1
Bounds
0 <= X1 <= +infinity
End
"@ | Set-Content -LiteralPath $lpPath -Encoding ASCII

docker build -t $Image (Join-Path $repoRoot "docker/soplex") | Out-Host

$mountPath = (Resolve-Path -LiteralPath $runPath).Path
docker run --rm -v "${mountPath}:/work" -w /work $Image `
    sh -lc "soplex -s0 -g5 -t60 -f1e-9 -o1e-9 -x -q -c --int:readmode=0 solver_smoke.lp > solver_smoke.lp.out"

$out = Get-Content -LiteralPath $outPath -Raw
if ($out -notmatch "problem is solved \[optimal\]" -or $out -notmatch "Objective value") {
    throw "SoPlex smoke solve did not report optimal objective. See $outPath"
}

Write-Host "SoPlex solver smoke passed: $outPath"

