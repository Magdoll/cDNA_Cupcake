from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [#Extension("BioReaders", ["BioReaders.pyx"]), \
        Extension("c_branch", ["c_branch.pyx"]),\
        Extension("intersection_unique", ["intersection_unique.pyx"])]

setup(
		name = 'tofu_branch',
        cmdclass = {'build_ext': build_ext},
		ext_modules = ext_modules,
        author_email='etseng@pacificbiosciences.com',
        author='Elizabeth Tseng'
)

