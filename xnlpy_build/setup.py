from distutils.core import setup, Extension
import glob
import pathlib

c_files = [file for file in glob.glob("core/*.c")]
here = pathlib.Path(__file__).parent.resolve()
long_description = (here / "../README.md").read_text(encoding="utf-8")

module1 = Extension(
	'xnlpy', # name
	# define_macros = [('PY_SSIZE_T_CLEAN', None)],
	# it's not necessary to put the Python lib here, because it's already on the PATH
	sources = c_files,
	extra_compile_args = ['-O3','-funroll-loops', '-march=native', '-mfpmath=sse'],
	)

setup(
	name = 'xnlpy',
	version = '1.0',
	license='MIT',
	description = 'XNL library for python!',
	author = 'Gilberto Jose Guimaraes de Sousa Mourao',
	author_email = 'gilbertojos.mourao@gmail.com',
	url = 'https://github.com/gilbertomourao/xnlpy',
	platforms = 'Windows 10, Ubuntu 22.04',
	long_description = long_description,
	ext_modules = [module1],
	include_package_data=True
	);