param(
    [switch]$Stop,
    [int]$MinAgeMinutes = 30,
    [string]$RepoRoot = ""
)

$ErrorActionPreference = "Stop"

if ([string]::IsNullOrWhiteSpace($RepoRoot)) {
    $RepoRoot = Resolve-Path (Join-Path $PSScriptRoot "..")
} else {
    $RepoRoot = Resolve-Path $RepoRoot
}
$repoRootText = ([string]$RepoRoot).TrimEnd('\') + '\'

$cwdReaderSource = @'
using System;
using System.Text;
using System.Runtime.InteropServices;

public static class PcSecProcessCwdReader {
  [StructLayout(LayoutKind.Sequential)]
  public struct PROCESS_BASIC_INFORMATION {
    public IntPtr Reserved1;
    public IntPtr PebBaseAddress;
    public IntPtr Reserved2_0;
    public IntPtr Reserved2_1;
    public IntPtr UniqueProcessId;
    public IntPtr Reserved3;
  }
  [StructLayout(LayoutKind.Sequential)]
  public struct UNICODE_STRING64 {
    public ushort Length;
    public ushort MaximumLength;
    public ulong Buffer;
  }
  [DllImport("ntdll.dll")]
  static extern int NtQueryInformationProcess(IntPtr ProcessHandle, int ProcessInformationClass, ref PROCESS_BASIC_INFORMATION ProcessInformation, int ProcessInformationLength, out int ReturnLength);
  [DllImport("kernel32.dll", SetLastError=true)]
  static extern IntPtr OpenProcess(uint processAccess, bool bInheritHandle, int processId);
  [DllImport("kernel32.dll", SetLastError=true)]
  static extern bool ReadProcessMemory(IntPtr hProcess, IntPtr lpBaseAddress, byte[] lpBuffer, int dwSize, out IntPtr lpNumberOfBytesRead);
  [DllImport("kernel32.dll", SetLastError=true)]
  static extern bool CloseHandle(IntPtr hObject);
  const uint PROCESS_QUERY_INFORMATION = 0x0400;
  const uint PROCESS_VM_READ = 0x0010;
  static ulong ReadUInt64(IntPtr h, ulong addr) {
    byte[] b = new byte[8]; IntPtr n;
    if (!ReadProcessMemory(h, (IntPtr)addr, b, b.Length, out n)) return 0;
    return BitConverter.ToUInt64(b, 0);
  }
  static UNICODE_STRING64 ReadUnicodeString(IntPtr h, ulong addr) {
    byte[] b = new byte[16]; IntPtr n;
    if (!ReadProcessMemory(h, (IntPtr)addr, b, b.Length, out n)) return new UNICODE_STRING64();
    return new UNICODE_STRING64 { Length=BitConverter.ToUInt16(b,0), MaximumLength=BitConverter.ToUInt16(b,2), Buffer=BitConverter.ToUInt64(b,8) };
  }
  static string ReadString(IntPtr h, UNICODE_STRING64 us) {
    if (us.Buffer == 0 || us.Length == 0 || us.Length > 32766) return "";
    byte[] b = new byte[us.Length]; IntPtr n;
    if (!ReadProcessMemory(h, (IntPtr)us.Buffer, b, b.Length, out n)) return "";
    return Encoding.Unicode.GetString(b);
  }
  public static string GetCwd(int pid) {
    IntPtr h = OpenProcess(PROCESS_QUERY_INFORMATION | PROCESS_VM_READ, false, pid);
    if (h == IntPtr.Zero) return "";
    try {
      PROCESS_BASIC_INFORMATION pbi = new PROCESS_BASIC_INFORMATION(); int ret;
      int status = NtQueryInformationProcess(h, 0, ref pbi, Marshal.SizeOf(typeof(PROCESS_BASIC_INFORMATION)), out ret);
      if (status != 0 || pbi.PebBaseAddress == IntPtr.Zero) return "";
      ulong peb = (ulong)pbi.PebBaseAddress.ToInt64();
      ulong processParameters = ReadUInt64(h, peb + 0x20);
      if (processParameters == 0) return "";
      return ReadString(h, ReadUnicodeString(h, processParameters + 0x38));
    } finally { CloseHandle(h); }
  }
}
'@

if (-not ([System.Management.Automation.PSTypeName]'PcSecProcessCwdReader').Type) {
    Add-Type -TypeDefinition $cwdReaderSource
}

$now = Get-Date
$processes = @(Get-CimInstance Win32_Process -Filter "name = 'python.exe'" -ErrorAction SilentlyContinue)
$processIds = @($processes | ForEach-Object { [int]$_.ProcessId })
$tcpOwners = @(Get-NetTCPConnection -ErrorAction SilentlyContinue | Select-Object -ExpandProperty OwningProcess -Unique)
$suspects = @()

foreach ($process in $processes) {
    $pidValue = [int]$process.ProcessId
    $commandLine = [string]$process.CommandLine
    $parentPid = [int]$process.ParentProcessId
    $parentMissing = $parentPid -gt 0 -and ($processIds -notcontains $parentPid) -and -not (Get-Process -Id $parentPid -ErrorAction SilentlyContinue)
    if ($process.CreationDate -is [datetime]) {
        $creationTime = $process.CreationDate
    } else {
        $creationTime = [Management.ManagementDateTimeConverter]::ToDateTime([string]$process.CreationDate)
    }
    $ageMinutes = ($now - $creationTime).TotalMinutes
    $cwd = [PcSecProcessCwdReader]::GetCwd($pidValue)
    $inRepo = $cwd.StartsWith($repoRootText, [StringComparison]::OrdinalIgnoreCase)
    $isStdinPython = $commandLine -match 'python(?:\.exe)?"?\s+-\s*$'
    $ownsTcp = $tcpOwners -contains $pidValue

    if ($inRepo -and $isStdinPython -and $parentMissing -and -not $ownsTcp -and $ageMinutes -ge $MinAgeMinutes) {
        $suspects += [pscustomobject]@{
            ProcessId = $pidValue
            ParentProcessId = $parentPid
            AgeMinutes = [math]::Round($ageMinutes, 1)
            WorkingDirectory = $cwd
            CommandLine = $commandLine
        }
    }
}

if ($suspects.Count -eq 0) {
    Write-Host "No stale pcSecYeastSpecies orphan python stdin processes found." -ForegroundColor Green
    exit 0
}

$suspects | Format-Table ProcessId, ParentProcessId, AgeMinutes, WorkingDirectory -AutoSize

if (-not $Stop) {
    Write-Host ""
    Write-Host "Dry run only. Re-run with -Stop to terminate these stale orphan compute processes." -ForegroundColor Yellow
    exit 0
}

foreach ($suspect in $suspects) {
    Stop-Process -Id $suspect.ProcessId -Force -ErrorAction Stop
    Write-Host "Stopped stale orphan python PID $($suspect.ProcessId)." -ForegroundColor Green
}
