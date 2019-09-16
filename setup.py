from setuptools import setup, Extension, find_packages
import sys
import numpy as np


__author__ = "etseng@pacb.com"
version = "8.5"

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
    packages = ['cupcake', 'cupcake.io', 'cupcake.ice',
                'cupcake.tofu', 'cupcake.tofu.branch', 'cupcake.tofu.counting',
                'cupcake2', 'cupcake2.io', 'cupcake2.ice2', 'cupcake2.tofu2',
                'phasing', 'phasing.io'],
    install_requires=[
        'biopython',
        'bx-python==0.7.3'
        ],
    scripts = ['cupcake/tofu/collapse_isoforms_by_sam.py',
               'cupcake/tofu/get_abundance_post_collapse.py',
               'cupcake/tofu/make_sam_by_isoforms.py',
               'cupcake/tofu/filter_by_count.py',
               'cupcake/tofu/filter_away_subset.py',
               'cupcake/tofu/process_blasr_to_read_stat.py', 
	       'cupcake/tofu/process_read_stats_to_count.py',
               'cupcake/tofu/get_counts_by_barcode.py',
               'cupcake/tofu/fusion_finder.py',
               'cupcake/tofu/counting/chain_samples.py',
               'cupcake/tofu/counting/chain_fusion_samples.py',
               'cupcake/tofu/counting/summarize_sample_GFF_junctions.py',
               'cupcake/tofu/counting/scrub_sample_GFF_junctions.py',
               'cupcake2/tofu2/ice_pbdagcon2.py',
               'cupcake2/tofu2/run_preCluster.py',
               'cupcake2/tofu2/run_IceInit2.py',
               'cupcake2/tofu2/run_IceIterative2.py',
               'cupcake2/tofu2/run_IcePartial2.py',
               'cupcake2/tofu2/run_IceArrow2.py',
               'cupcake2/io/SeqSplitter.py',
               'cupcake2/tofu2/generate_batch_cmd_for_preCluster_out.py',
               'cupcake2/tofu2/generate_batch_cmd_for_polishing.py',
               'cupcake2/tofu2/collect_IceIterative2_result.py',
               'cupcake2/tofu2/picking_up_ice2.py',
               'phasing/create_fake_genome.py',
               'phasing/run_phaser.py'
               ],
    )
