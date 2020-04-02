rem running build
python setup.py build --compiler=mingw32

rem creating the folder xnlpy
cd C:\Python\Python38-32\Lib\site-packages
mkdir xnlpy

rem copying data
cd C:\Users\chico\AppData\Roaming\Sublime Text 3\Packages\User\Python 3\git_xnlpy\xnlpy\build\lib.win32-3.8
copy calculus* "C:\Python\Python38-32\Lib\site-packages\xnlpy"
cd ../..
copy __init__.py "C:\Python\Python38-32\Lib\site-packages\xnlpy"