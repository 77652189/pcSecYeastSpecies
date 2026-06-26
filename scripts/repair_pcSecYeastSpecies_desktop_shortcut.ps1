$ErrorActionPreference = "Stop"

$projectRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$launcherPath = Join-Path $projectRoot "start_pcSecYeastSpecies_lan.bat"
$desktopPath = [Environment]::GetFolderPath("Desktop")

if (-not (Test-Path -LiteralPath $launcherPath)) {
    throw "Missing LAN launcher: $launcherPath"
}

$existingShortcut = Get-ChildItem -LiteralPath $desktopPath -Filter "*pcSecYeastSpecies*8502*.lnk" -ErrorAction SilentlyContinue | Select-Object -First 1
if ($existingShortcut) {
    $shortcutPath = $existingShortcut.FullName
} else {
    $shortcutPath = Join-Path $desktopPath "pcSecYeastSpecies LAN 8502.lnk"
}

$shell = New-Object -ComObject WScript.Shell
$shortcut = $shell.CreateShortcut($shortcutPath)
$shortcut.TargetPath = $launcherPath
$shortcut.WorkingDirectory = $projectRoot
$shortcut.Arguments = ""
$shortcut.Description = "Start pcSecYeastSpecies Streamlit LAN service on port 8502"

try {
    $shortcut.Save()
} catch {
    Write-Host "Failed to save the desktop shortcut." -ForegroundColor Red
    Write-Host "This usually means the current sandbox cannot write to Desktop." -ForegroundColor Yellow
    Write-Host "Run this command in a normal PowerShell window:" -ForegroundColor Yellow
    Write-Host "  powershell -NoProfile -ExecutionPolicy Bypass -File `"$PSCommandPath`"" -ForegroundColor Yellow
    throw
}

Write-Host "Desktop shortcut repaired:" -ForegroundColor Green
Write-Host "  $shortcutPath" -ForegroundColor Green
Write-Host "Target:" -ForegroundColor Green
Write-Host "  $launcherPath" -ForegroundColor Green