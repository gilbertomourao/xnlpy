@echo off

set python_path=C:\Python\Python38-32\Lib\site-packages
set build_dir=%CD%

rem running build
if exist build (
	rmdir /q /s build
)
python setup.py build --compiler=mingw32

rem creating the folder xnlpy
echo Checking if xnlpy already existis...
cd %python_path%
if exist xnlpy (
	echo Deleting the file...
	rmdir /q /s xnlpy
)
echo Creating xnlpy...
mkdir xnlpy

rem copying data
echo Moving files to the Python directory...
cd %build_dir%\build\lib.win32-3.8
echo Moving the compiled module...
copy xnlpy* "%python_path%\xnlpy"
cd ../..
echo Moving __init__.py...
copy __init__.py "%python_path%\xnlpy"