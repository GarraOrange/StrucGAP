$TOOLS_DIR = "D:\Research\tools"
$PHANTOM_ZIP = "$TOOLS_DIR\phantomjs-2.1.1-windows.zip"
$PHANTOM_DIR = "$TOOLS_DIR\phantomjs-2.1.1-windows"
$PHANTOM_BIN = "$PHANTOM_DIR\bin"

New-Item -ItemType Directory -Force $TOOLS_DIR | Out-Null

$urls = @(
    "https://bitbucket.org/ariya/phantomjs/downloads/phantomjs-2.1.1-windows.zip",
    "https://mirrors.huaweicloud.com/phantomjs/phantomjs-2.1.1-windows.zip"
)

foreach ($url in $urls) {
    try {
        Write-Host "Trying download:" $url
        Invoke-WebRequest -Uri $url -OutFile $PHANTOM_ZIP
        break
    } catch {
        Write-Host "Failed:" $url
    }
}

if (!(Test-Path $PHANTOM_ZIP)) {
    throw "PhantomJS download failed."
}

Expand-Archive -Force $PHANTOM_ZIP -DestinationPath $TOOLS_DIR

$env:Path = "$PHANTOM_BIN;$env:Path"

$userPath = [Environment]::GetEnvironmentVariable("Path", "User")
if ($userPath -notlike "*$PHANTOM_BIN*") {
    [Environment]::SetEnvironmentVariable("Path", "$PHANTOM_BIN;$userPath", "User")
}

phantomjs --version

$ASSET_DIR = "D:\Research\tools\pyecharts-assets"
New-Item -ItemType Directory -Force $ASSET_DIR | Out-Null

Invoke-WebRequest `
  -Uri "https://assets.pyecharts.org/assets/echarts.min.js" `
  -OutFile "$ASSET_DIR\echarts.min.js"

Get-Item "$ASSET_DIR\echarts.min.js"
