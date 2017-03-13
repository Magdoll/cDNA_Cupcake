from setuptools import setup, Extension, find_packages
import sys
import numpy as np


__author__ = "etseng@pacb.com"
version = "2.2"

ext_modules = [
                Extension("cupcake.tofu.branch.intersection_unique",
                    ["cupcake/tofu/branch/intersection_unique.c"]),
                Extension("cupcake.tofu.branch.c_branch",
                         ["cupcake/tofu/branch/c_branch.c"]),
                Extension("cupcake.ice.find_ECE",
                    ["cupcake/ice/find_ECE.c"]),
              ]


setup(
    name = 'cupcake',
    version=version,
    author='Elizabeth Tseng',
    author_email='etseng@pacb.com',
    ext_modules = ext_modules,
    include_dirs = [np.get_include()],
    zip_safe=False,
    packages = ['cupcake', 'cupcake.io', 'cupcake.ice', 'cupcake.tofu', 'cupcake.tofu.branch'],
    install_requires=[
        'biopython',
        'bx-python'
        ],
    scripts = ['cupcake/tofu/collapse_isoforms_by_sam.py',
               'cupcake/tofu/get_abundance_post_collapse.py',
               'cupcake/tofu/make_sam_by_isoforms.py',
               'cupcake/tofu/filter_by_count.py',
               'cupcake/tofu/filter_away_subset.py',
               'cupcake/tofu/process_blasr_to_read_stat.py',
               'cupcake/tofu/process_read_stats_to_count.py'
               ],
    )
