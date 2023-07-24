from setuptools import setup, Extension, find_packages
import sys
import numpy as np
from Cython.Build import cythonize


__author__ = "etseng@pacb.com"
version = "29.0.0"

ext_modules = [
                Extension("cupcake.tofu.branch.intersection_unique",
                    ["cupcake/tofu/branch/intersection_unique.pyx"]),
                Extension("cupcake.tofu.branch.c_branch",
                         ["cupcake/tofu/branch/c_branch.pyx"]),
                Extension("cupcake.ice.find_ECE",
                    ["cupcake/ice/find_ECE.pyx"]),
              ]


setup(
    name = 'cupcake',
    version=version,
    author='Elizabeth Tseng',
    author_email='etseng@pacb.com',
    ext_modules = cythonize(ext_modules, language_level = "2"),
    include_dirs = [np.get_include()],
    zip_safe=False,
    packages = ['cupcake', 'cupcake.io', 'cupcake.ice',
                'cupcake.tofu', 'cupcake.tofu.branch', 'cupcake.tofu.counting',
                'phasing', 'phasing.io'],
    setup_requires=[
        'numpy',
        'cython'
        ],
    install_requires=[
        'biopython',
        'bx-python>=0.7.3',
        'numpy',
        'bcbio-gff',
        'sklearn',
        'pysam'
        ],
    scripts = ['cupcake/tofu/collapse_isoforms_by_sam.py',
               'cupcake/tofu/get_abundance_post_collapse.py',
               'cupcake/tofu/filter_by_count.py',
               'cupcake/tofu/filter_away_subset.py',
               'cupcake/tofu/fusion_finder.py',
               'cupcake/tofu/fusion_collate_info.py',
			   'cupcake/tofu/color_bed12_post_sqanti.py',
               'cupcake/tofu/counting/chain_samples.py',
               'cupcake/tofu/counting/chain_fusion_samples.py',
               'cupcake/tofu/counting/summarize_sample_GFF_junctions.py',
               'cupcake/tofu/counting/scrub_sample_GFF_junctions.py',
               'cupcake/tofu/simple_stats_post_collapse.py',
               'phasing/create_fake_genome.py',
               'phasing/run_phaser.py',
			   'phasing/mag_phaser.py',
               'phasing/utils/paint_bam_post_phaser.py'
               ],
    )
