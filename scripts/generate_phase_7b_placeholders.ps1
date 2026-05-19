# Step 12 Phase 7b.7 — Generate placeholder PNGs for the workflow-panel
# UI screenshots directory. Run from the repository root:
#
#     powershell -ExecutionPolicy Bypass -File scripts/generate_phase_7b_placeholders.ps1
#
# Each placeholder is a 1280x720 mid-grey PNG with two centered text
# rows: a bold red "PENDING MANUAL CAPTURE" banner and the target
# filename in monospace. The placeholders make Markdown references in
# the Phase 8 reports resolve to a real image rather than a broken
# link; real captures (per the README's procedure) overwrite the PNG
# files in place.
#
# Requires Windows PowerShell 5.1 or PowerShell 7+ on Windows
# (System.Drawing is .NET Framework on PS 5.1, .NET on PS 7 with
# `System.Drawing.Common` shipped via the SDK).

Add-Type -AssemblyName System.Drawing

$outDir = Join-Path $PSScriptRoot "..\docs\reports\step12_phase_7b_screenshots"
$outDir = [System.IO.Path]::GetFullPath($outDir)
if (-not (Test-Path $outDir)) {
    New-Item -ItemType Directory -Path $outDir -Force | Out-Null
}

$names = @(
    "01_panel_idle_pre_run.png",
    "02_phase_a_running_cycle_3_of_5.png",
    "03_phase_a_completed.png",
    "04_phase_b_hd_output_ready.png"
)

foreach ($name in $names) {
    $bmp   = New-Object System.Drawing.Bitmap 1280, 720
    $g     = [System.Drawing.Graphics]::FromImage($bmp)
    $g.SmoothingMode     = [System.Drawing.Drawing2D.SmoothingMode]::AntiAlias
    $g.TextRenderingHint = [System.Drawing.Text.TextRenderingHint]::AntiAlias
    $g.Clear([System.Drawing.Color]::FromArgb(192, 192, 192))

    $bannerFont = New-Object System.Drawing.Font "Arial", 56, ([System.Drawing.FontStyle]::Bold)
    $nameFont   = New-Object System.Drawing.Font "Consolas", 28
    $hintFont   = New-Object System.Drawing.Font "Arial", 18, ([System.Drawing.FontStyle]::Italic)

    $bannerBrush = New-Object System.Drawing.SolidBrush([System.Drawing.Color]::FromArgb(160, 32, 32))
    $nameBrush   = New-Object System.Drawing.SolidBrush([System.Drawing.Color]::FromArgb(48, 48, 48))
    $hintBrush   = New-Object System.Drawing.SolidBrush([System.Drawing.Color]::FromArgb(96, 96, 96))

    $banner = "PENDING MANUAL CAPTURE"
    $hint   = "Capture procedure: docs/reports/step12_phase_7b_screenshots/README.md"

    $bannerSize = $g.MeasureString($banner, $bannerFont)
    $nameSize   = $g.MeasureString($name,    $nameFont)
    $hintSize   = $g.MeasureString($hint,    $hintFont)

    $bannerX = (1280 - $bannerSize.Width) / 2
    $nameX   = (1280 - $nameSize.Width)   / 2
    $hintX   = (1280 - $hintSize.Width)   / 2

    $g.DrawString($banner, $bannerFont, $bannerBrush, $bannerX, 280)
    $g.DrawString($name,   $nameFont,   $nameBrush,   $nameX,   380)
    $g.DrawString($hint,   $hintFont,   $hintBrush,   $hintX,   460)

    $border = New-Object System.Drawing.Pen([System.Drawing.Color]::FromArgb(120, 120, 120)), 4
    $g.DrawRectangle($border, 2, 2, 1276, 716)

    $outPath = Join-Path $outDir $name
    $bmp.Save($outPath, [System.Drawing.Imaging.ImageFormat]::Png)

    $g.Dispose()
    $bmp.Dispose()
    $bannerFont.Dispose()
    $nameFont.Dispose()
    $hintFont.Dispose()
    $bannerBrush.Dispose()
    $nameBrush.Dispose()
    $hintBrush.Dispose()
    $border.Dispose()

    Write-Host "wrote $outPath"
}
