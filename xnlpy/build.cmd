@echo off

rem set build_dir=%CD%

rem finding python
for /f "tokens=*" %%a in ('python -c "import sys, os; print(os.path.dirname(sys.executable))"') do (
	set python_path="%%a\Lib\site-packages"
)

rem running build

echo.

echo ---------------------------------------------------------------------
echo Compiling the source
echo ---------------------------------------------------------------------
echo.

if exist build (
	rmdir /q /s build
)
python setup.py build --compiler=mingw32

echo.
echo ---------------------------------------------------------------------
echo Creating XNLPY
echo ---------------------------------------------------------------------
echo.

rem creating the folder xnlpy
echo Checking if xnlpy already exists...
rem cd %python_path%
if exist %python_path%\xnlpy (
	echo Deleting the file...
	rmdir /q /s %python_path%\xnlpy
)
echo Creating xnlpy...
mkdir %python_path%\xnlpy

echo.
echo ---------------------------------------------------------------------
echo Moving files
echo ---------------------------------------------------------------------
echo.

rem copying data
echo Moving files to the Python directory...
rem cd %build_dir%\build\lib.win32-3.8
echo Moving the compiled module...
copy build\lib.win32-3.8\xnlpy* "%python_path%\xnlpy"
rem cd ../..
echo Moving __init__.py...
copy __init__.py "%python_path%\xnlpy"

echo.
echo ---------------------------------------------------------------------
echo Done.
