@echo off
setlocal

rem ===== Configuration =====
set "EXE=C:\Users\Nikita-PC\Desktop\G4mingle_neutrongen\rez2\RYZENBENCH\mingle.exe"
set "ARGS=beamOn_100000.mac"
set "OUTFILE=C:\Users\Nikita-PC\Desktop\G4mingle_neutrongen\rez2\RYZENBENCH\elapsed_ms.txt"
rem =========================

powershell -NoProfile -Command ^
  "$sw = [System.Diagnostics.Stopwatch]::StartNew(); " ^
  "& '%EXE%' %ARGS%; " ^
  "$exitCode = $LASTEXITCODE; " ^
  "$sw.Stop(); " ^
  "Set-Content -Path '%OUTFILE%' -Value $sw.ElapsedMilliseconds; " ^
  "Write-Host ('Elapsed time: ' + $sw.ElapsedMilliseconds + ' ms'); " ^
  "exit $exitCode"

exit /b %ERRORLEVEL%