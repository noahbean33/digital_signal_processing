@echo off
REM Build script for Signal Processing Library (Windows)
REM Requires MinGW GCC or Visual Studio compiler

echo ========================================
echo Signal Processing Library - Build Script
echo ========================================
echo.

REM Check if bin directory exists, create if not
if not exist bin mkdir bin

REM Compiler settings
set CC=gcc
set CFLAGS=-Wall -O2 -Iinclude
set SRC=src\signal_io.c src\dft.c src\signal_stats.c src\signal_conv.c src\signal_utils.c

echo Building examples...
echo.

REM Build DFT example
echo [1/3] Building example_dft...
%CC% %CFLAGS% -o bin\example_dft.exe examples\example_dft.c %SRC%
if %errorlevel% neq 0 (
    echo ERROR: Failed to build example_dft
    goto error
)
echo       SUCCESS: bin\example_dft.exe

REM Build convolution example
echo [2/3] Building example_convolution...
%CC% %CFLAGS% -o bin\example_convolution.exe examples\example_convolution.c %SRC%
if %errorlevel% neq 0 (
    echo ERROR: Failed to build example_convolution
    goto error
)
echo       SUCCESS: bin\example_convolution.exe

REM Build stats/utils example
echo [3/3] Building example_stats_utils...
%CC% %CFLAGS% -o bin\example_stats_utils.exe examples\example_stats_utils.c %SRC%
if %errorlevel% neq 0 (
    echo ERROR: Failed to build example_stats_utils
    goto error
)
echo       SUCCESS: bin\example_stats_utils.exe

echo.
echo ========================================
echo Build completed successfully!
echo ========================================
echo.
echo Executables are in the bin\ directory:
echo   - bin\example_dft.exe
echo   - bin\example_convolution.exe
echo   - bin\example_stats_utils.exe
echo.
echo To run examples:
echo   cd bin
echo   example_dft.exe
echo   example_convolution.exe
echo   example_stats_utils.exe
echo.
goto end

:error
echo.
echo ========================================
echo Build failed!
echo ========================================
echo.
echo Make sure GCC (MinGW) is installed and in your PATH.
echo To install MinGW: https://www.mingw-w64.org/
echo.
exit /b 1

:end
pause
