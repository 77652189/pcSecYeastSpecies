param([switch]$ConfigureFirewallOnly)

$ErrorActionPreference = "Stop"

$projectName = "pcSecYeastSpecies"
$port = 8502
$ruleName = "pcSecYeastSpecies Streamlit 8502 LAN"
$projectRoot = Resolve-Path (Join-Path $PSScriptRoot "..")
$appPath = "app/ui/streamlit_app.py"

function Test-IsAdministrator {
    $identity = [Security.Principal.WindowsIdentity]::GetCurrent()
    $principal = New-Object Security.Principal.WindowsPrincipal($identity)
    return $principal.IsInRole([Security.Principal.WindowsBuiltInRole]::Administrator)
}

function Test-FirewallRuleReady {
    $rule = Get-NetFirewallRule -DisplayName $ruleName -ErrorAction SilentlyContinue | Select-Object -First 1
    if (-not $rule -or $rule.Enabled -ne "True" -or $rule.Direction -ne "Inbound" -or $rule.Action -ne "Allow") {
        return $false
    }

    $portFilter = $rule | Get-NetFirewallPortFilter
    $addressFilter = $rule | Get-NetFirewallAddressFilter
    return ($portFilter.Protocol -eq "TCP" -and $portFilter.LocalPort -eq "$port" -and $addressFilter.RemoteAddress -eq "LocalSubnet")
}

function Ensure-FirewallRule {
    if (-not (Test-IsAdministrator)) {
        Write-Host "Administrator permission is required to configure the Windows firewall rule. Requesting UAC..." -ForegroundColor Yellow
        $args = @("-NoProfile", "-ExecutionPolicy", "Bypass", "-File", "`"$PSCommandPath`"", "-ConfigureFirewallOnly")
        Start-Process -FilePath "powershell.exe" -ArgumentList $args -Verb RunAs -Wait
        return
    }

    $rules = Get-NetFirewallRule -DisplayName $ruleName -ErrorAction SilentlyContinue
    if (-not $rules) {
        New-NetFirewallRule -DisplayName $ruleName -Direction Inbound -Action Allow -Protocol TCP -LocalPort $port -RemoteAddress LocalSubnet -Profile Any -Enabled True | Out-Null
    } else {
        $rules | Set-NetFirewallRule -Enabled True -Direction Inbound -Action Allow -Profile Any
        $rules | Get-NetFirewallPortFilter | Set-NetFirewallPortFilter -Protocol TCP -LocalPort $port
        $rules | Get-NetFirewallAddressFilter | Set-NetFirewallAddressFilter -RemoteAddress LocalSubnet
    }
}

function Get-PortOwner {
    $connection = Get-NetTCPConnection -LocalPort $port -State Listen -ErrorAction SilentlyContinue | Select-Object -First 1
    if (-not $connection) {
        return $null
    }
    return Get-CimInstance Win32_Process -Filter "ProcessId=$($connection.OwningProcess)" -ErrorAction SilentlyContinue
}

function Test-HealthOk {
    try {
        $response = Invoke-WebRequest -Uri "http://127.0.0.1:$port/_stcore/health" -UseBasicParsing -TimeoutSec 5
        return ($response.Content.Trim() -eq "ok")
    } catch {
        return $false
    }
}

function Write-LanUrls {
    $lanAddresses = Get-NetIPAddress -AddressFamily IPv4 |
        Where-Object { $_.IPAddress -notlike "127.*" -and $_.IPAddress -notlike "169.254.*" -and $_.InterfaceAlias -notlike "*WSL*" } |
        Select-Object -ExpandProperty IPAddress

    foreach ($address in $lanAddresses) {
        Write-Host "  http://$address`:$port" -ForegroundColor Green
    }
}

if ($ConfigureFirewallOnly) {
    Ensure-FirewallRule
    exit
}

if (-not (Test-FirewallRuleReady)) {
    Ensure-FirewallRule
}

Set-Location $projectRoot

$owner = Get-PortOwner
if ($owner) {
    $commandLine = [string]$owner.CommandLine
    $isThisApp = ($commandLine -like "*app/ui/streamlit_app.py*" -or $commandLine -like "*app\ui\streamlit_app.py*")
    if ($isThisApp -and (Test-HealthOk)) {
        Write-Host ""
        Write-Host "$projectName is already running. No need to start it again." -ForegroundColor Green
        Write-Host "LAN URLs:" -ForegroundColor Green
        Write-LanUrls
        Write-Host ""
        Write-Host "Local URL: http://127.0.0.1:$port" -ForegroundColor Green
        Write-Host ""
        exit 0
    }

    Write-Host ""
    Write-Host "Port $port is already in use. Cannot start $projectName." -ForegroundColor Red
    Write-Host "Owner PID: $($owner.ProcessId)" -ForegroundColor Yellow
    Write-Host "Owner command: $commandLine" -ForegroundColor Yellow
    Write-Host ""
    Write-Host "If this is an old project service, close that window or stop that process first." -ForegroundColor Yellow
    exit 1
}

$python = Join-Path $projectRoot ".venv\Scripts\python.exe"
if (-not (Test-Path $python)) {
    $python = "python"
}

Write-Host ""
Write-Host "$projectName will start on the LAN:" -ForegroundColor Green
Write-LanUrls
Write-Host ""
Write-Host "To stop the service, close this window or press Ctrl+C." -ForegroundColor Yellow
Write-Host ""

& $python -m streamlit run $appPath --server.address=0.0.0.0 --server.port=$port
