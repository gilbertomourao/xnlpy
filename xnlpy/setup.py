from distutils.core import setup, Extension
import glob

c_files = [file for file in glob.glob("src/*.c")]

module1 = Extension(
	'xnlpy', # name
	# define_macros = [('PY_SSIZE_T_CLEAN', None)],
	# it's not necessary to put the Python lib here, because it's already on the PATH
	sources = c_files,
	extra_compile_args = ["-O3"],
	)

setup(
	name = 'xnlpy',
	version = '1.0',
	description = 'XNL library for python!',
	author = 'Gilberto Jose Guimaraes de Sousa Mourao',
	author_email = 'gilbertojos.mourao@gmail.com',
	url = 'https://github.com/gilbertomourao/xnlpy',
	platforms = 'Windows 10',
	long_description = 'XNL is a C library that contains some numerical algorithms.',
	ext_modules = [module1]
	);