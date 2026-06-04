param(
    [string]$Address = "0.0.0.0",
    [int]$Port = 8501
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
streamlit run $appPath --server.address $Address --server.port $Port
