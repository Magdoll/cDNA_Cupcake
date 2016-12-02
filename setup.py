from setuptools import setup, Extension, find_packages
import sys
import numpy as np


__author__ = "etseng@pacb.com"
version = "2.1"

ext_modules = [
                Extension("cupcake.tofu.branch.intersection_unique",
                    ["cupcake/tofu/branch/intersection_unique.c"]),
                Extension("cupcake.tofu.branch.c_branch",
                         ["cupcake/tofu/branch/c_branch.c"]),
              ]


setup(
    name = 'cupcake',
    version=version,
    author='Elizabeth Tseng',
    author_email='etseng@pacb.com',
    ext_modules = ext_modules,
    include_dirs = [np.get_include()],
    zip_safe=False,
    packages = ['cupcake', 'cupcake.io', 'cupcake.tofu', 'cupcake.tofu.branch', 'cupcake.tofu.counting'],
    install_requires=[
        'biopython',
        'bx-python'
        ],
    scripts = ['cupcake/tofu/collapse_isoforms_by_sam.py',
               'cupcake/tofu/get_abundance_post_collapse.py',
               'cupcake/tofu/make_sam_by_isoforms.py',
               'cupcake/tofu/filter_by_count.py',
               'cupcake/tofu/filter_away_subset.py',
               'cupcake/tofu/counting/chain_samples.py'
               ],
    )
