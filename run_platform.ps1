$ErrorActionPreference = "Stop"
$projectRoot = $PSScriptRoot
$env:PYTHONPATH = Join-Path $projectRoot "src"
Set-Location -LiteralPath $projectRoot
python -m litrapid.platform_server --host 127.0.0.1 --port 8765
