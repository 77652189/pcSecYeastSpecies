param(
    [string]$Distro = "Ubuntu-24.04"
)

$ErrorActionPreference = "Stop"

Write-Host "Target WSL distro: $Distro"

$installed = (& wsl.exe -l -q) -replace "`0", "" | Where-Object { $_ -eq $Distro }
if (-not $installed) {
    Write-Host "Installing $Distro. If this stalls, install it manually with: wsl --install -d $Distro"
    & wsl.exe --install -d $Distro --no-launch
}

Write-Host "Initializing/checking $Distro..."
& wsl.exe -d $Distro -u root -- sh -lc "cat /etc/os-release"

Write-Host "Installing SoPlex package if available..."
& wsl.exe -d $Distro -u root -- sh -lc "apt-get update && (apt-cache search '^soplex$' || true) && apt-get install -y soplex"

Write-Host "Verifying SoPlex..."
& wsl.exe -d $Distro -u root -- sh -lc "command -v soplex && soplex --version"

