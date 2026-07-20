param(
    [string]$Address = "0.0.0.0",
    [int]$Port = 8502
)

$ErrorActionPreference = "Stop"

$repoRoot = $PSScriptRoot
$appPath = Join-Path $repoRoot "app/ui/streamlit_app.py"

if (-not (Test-Path -LiteralPath $appPath)) {
    throw "Missing Streamlit app: $appPath"
}

if (-not (Get-Command streamlit -ErrorAction SilentlyContinue)) {
    throw "streamlit was not found. Install dependencies with: pip install -r requirements.txt"
}

Set-Location $repoRoot

# `streamlit run` (the console-script shim) puts the shim's own directory on
# sys.path[0], not the repo root, so `import app.ui.common` fails with
# ModuleNotFoundError. `python -m streamlit` puts the current directory (the
# repo root, since we just Set-Location'd there) on sys.path[0] instead, which
# is what `app.*` absolute imports need.
$pathEntries = @($repoRoot, (Join-Path $repoRoot "python_pichia\src"))
$existing = @()
if (-not [string]::IsNullOrWhiteSpace($env:PYTHONPATH)) {
    $existing = $env:PYTHONPATH -split [IO.Path]::PathSeparator
}
foreach ($entry in $pathEntries) {
    if ((Test-Path -LiteralPath $entry) -and ($existing -notcontains $entry)) {
        $existing = @($entry) + $existing
    }
}
$env:PYTHONPATH = ($existing -join [IO.Path]::PathSeparator)
python -m streamlit run $appPath --server.address $Address --server.port $Port
