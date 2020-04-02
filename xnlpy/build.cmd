set python_path=C:\Python\Python38-32\Lib\site-packages

rem running build
if exist build (
	rmdir /q /s build
)
python setup.py build --compiler=mingw32

rem creating the folder xnlpy
cd %python_path%
if exist xnlpy (
	rmdir /q /s xnlpy
)
mkdir xnlpy

rem copying data
cd C:\Users\chico\AppData\Roaming\Sublime Text 3\Packages\User\Python 3\git_xnlpy\xnlpy\build\lib.win32-3.8
copy xnlpy* "%python_path%\xnlpy"
cd ../..
copy __init__.py "%python_path%\xnlpy"