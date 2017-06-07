from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [
        Extension("find_ECE", ["find_ECE.pyx"]),
        ]

setup(
		name = 'c_ICE',
		cmdclass = {'build_ext': build_ext},
		ext_modules = ext_modules,
        author_email='etseng@pacificbiosciences.com',
        author='Elizabeth Tseng'
)


