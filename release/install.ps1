[CmdletBinding()]
param(
    [string]$Prefix = $(if ($env:RUSQSIEVE_PREFIX) {
        $env:RUSQSIEVE_PREFIX
    } else {
        Join-Path $env:ProgramFiles "rusqsieve"
    }),
    [switch]$Elevated
)

$ErrorActionPreference = "Stop"
$sourceRoot = Split-Path -LiteralPath $PSCommandPath -Parent
$principal = [Security.Principal.WindowsPrincipal]::new(
    [Security.Principal.WindowsIdentity]::GetCurrent()
)
$isAdministrator = $principal.IsInRole(
    [Security.Principal.WindowsBuiltInRole]::Administrator
)

if (-not $isAdministrator) {
    if ($Elevated) {
        throw "Administrator access was requested but was not granted."
    }

    Write-Host "Administrator access is required to install under $Prefix."
    $arguments = @(
        "-NoLogo"
        "-NoProfile"
        "-ExecutionPolicy"
        "Bypass"
        "-File"
        "`"$PSCommandPath`""
        "-Elevated"
        "-Prefix"
        "`"$Prefix`""
    )
    $process = Start-Process `
        -FilePath "powershell.exe" `
        -ArgumentList $arguments `
        -Verb RunAs `
        -Wait `
        -PassThru
    exit $process.ExitCode
}

$binDir = Join-Path $Prefix "bin"
$libDir = Join-Path $Prefix "lib"
$includeDir = Join-Path $Prefix "include"
$pkgConfigDir = Join-Path $libDir "pkgconfig"

New-Item -ItemType Directory -Force -Path @(
    $binDir
    $libDir
    $includeDir
    $pkgConfigDir
) | Out-Null

Copy-Item -Force -Path (Join-Path $sourceRoot "bin\*") -Destination $binDir
Copy-Item -Force -Path (Join-Path $sourceRoot "lib\*.lib") -Destination $libDir
Copy-Item -Force `
    -LiteralPath (Join-Path $sourceRoot "include\rusqsieve.h") `
    -Destination $includeDir

$pcSource = Join-Path $sourceRoot "lib\pkgconfig\rusqsieve.pc"
$pcDestination = Join-Path $pkgConfigDir "rusqsieve.pc"
$pcPrefix = $Prefix.Replace("\", "/")
$pc = [IO.File]::ReadAllText($pcSource)
$pc = [Text.RegularExpressions.Regex]::Replace(
    $pc,
    "(?m)^prefix=.*$",
    "prefix=$pcPrefix"
)
[IO.File]::WriteAllText(
    $pcDestination,
    $pc,
    [Text.Encoding]::ASCII
)

Write-Host "Installed rusqsieve under $Prefix"
Write-Host "Add $binDir to PATH to run qs-factor from any command prompt."
