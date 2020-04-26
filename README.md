# cDNA_Cupcake

![logo](https://github.com/Magdoll/images_public/blob/master/logos/Cupcake_logo.png)

Last Updated: 04.26.2020

**cDNA_Cupcake** is a miscellaneous collection of Python and R scripts used for analyzing sequencing data. Most of the scripts only require [Biopython](http://biopython.org/wiki/Download). For scripts that require additional libraries, it will be specified in documentation.

Current version: 12.0.0

**If you use Python 2.7, please use the branch [Py2_v8.7.x](https://github.com/Magdoll/cDNA_Cupcake/tree/Py2_v8.7.x). Note that however, the Py2 version will NOT have any new features, only bug fixes! Users are encouraged to migrate to Python 3.7 now.**


## Python Requirements
* Python >= 3.7
* Biopython 


If you are running post Iso-Seq analysis you will need more install requirements. See [wiki](https://github.com/Magdoll/cDNA_Cupcake/wiki/Cupcake-ToFU%3A-supporting-scripts-for-Iso-Seq-after-clustering-step) for more details.

## How to use this repository

Since most of the scripts are independent (do not depend on each other), you can either clone the whole directory, or, if you are only interested in a specific script, just download that specific script to your local drive.

You can clone the GitHub repository, then add the GitHub repo path to your `$PATH` variable. The scripts are organized into different sub-directories (ex: `sequence/`, `rarefaction/` etc) so you will have to add them individually.

```
git clone https://github.com/Magdoll/cDNA_Cupcake.git
export PATH=$PATH:<path_to_Cupcake>/sequence/
export PATH=$PATH:<path_to_Cupcake>/rarefaction/
```


For any issues or bugs, please report to [Issues](https://github.com/Magdoll/cDNA_Cupcake/issues).

## Documentation

Please see [wiki](https://github.com/Magdoll/cDNA_Cupcake/wiki) for the latest maintained list of scripts.

A brief list of currently listed scripts are:

### Annotation and Rarefaction
* `make_file_for_subsampling_from_collapsed.py`: Prepare file for running subsampling (rarefaction curve).
* `subsample.py` and `subsample_with_category.py`: Running subsamping. Results can be plotted with Excel graphics and R, etc.

### Sequence Manipulation
* `get_seq_stats.py`: Summarize length distribution of a FASTA/FASTQ file.
* `rev_comp.py`: Reverse complement a sequence from command line.
* `fa2fq.py` and `fq2fa.py`: Convert between FASTA and FASTQ format.
* `get_seqs_from_list.py`: extract list of sequences given a fasta file and a list of IDs.
* `err_correct_w_genome.py`: generate fasta sequences given genom
* `sam_to_bam.py`: quick script to run SAM to BAM conversion. Assumes `samtools` is installed.
* `sam_to_gff3.py`: use BCBio and BioPython to convert SAM file into GFF3 format. 
* `group_ORF_sequences.py`: group identical ORF sequences from different isoforms.

### Cupcake ToFU: supporting scripts for Iso Seq after clustering step
* `collapse_isoforms_by_sam.py`: Collapse HQ isoform results to unique isoforms (based on genome alignment).
* `get_abundance_post_collapse.py`: Obtain count information post collapse to unique isoforms.
* `filter_by_count.py`: Filter collapse result by FL count information.
* `filter_away_subset.py`: Filter away 5' degraded isoforms.
* `chain_samples.py`: Chaining together multiple samples.
* `fusion_finder.py`: Finding fusion genes.


## Version Changes

2020.04.26 updated to v12.0.0. `clip_out_UMI_cellBC.py` now supports new mode `G5-clip`.

2020.03.12 updated to v11.0.0. major chaining bug fixed for `chain_samples.py`.

<details>
   <summary>Click to see older version logs</summary>

	2020.01.31 updated to v10.0.1. fixed typo in `collapse_isoforms_by_sam.py`

	2020.01.27 updated to v10.0.0. `chain_samples.py` now supported multithreading via `--cpus` option, also fixed chain bugs related to 3' differences.

	2019.12.30 updated to v9.2.0. fixed support for running `run_preCluster.py` that is used by Cogent.

	2019.12.18 updated to v9.1.1. changed `bcbiogff` dependency to `bcbio-gff` in `setup.py`.

	2019.12.18 updated to v9.1.0. added fasta support for `filter_away_subset.py` and `filter_by_count.py`.

	2019.12.09 updated to v9.0.3. fixed geneid display issue with PB.X.Y in GFF output for collapse.

	2019.10.02 updated to v9.0.2. updated collapse script series parameters to fit with isoseq3 output.

	2019.10.02 updated to v9.0.1. bug fix on fasta output in `dedup_FLNC_per_cluster.py`. removed pbcore dependency so Py3 fully usable for single cell scripts!

	2019.09.25 updated to v8.7.0. `clip_out_UMI_cellBC.py` supports unusual 10X 5' end single cell schema.
    
    2019.09.24 updated to v8.6. `cupcake.io.GFF.py` now supports `gene_id` write out.
    
    2019.09.16 updated to v8.5. fixed `collapse_isoforms_by_sam.py` incorrect behavior in fuzzy chain
    
    2019.08.20 updated to v8.4. `run_phaser.py` dependncy is pyvcf, not bio-vcf.
    
    2019.08.19 updated to v8.3. removed bug testing code in cupcake preclustering.
    
    2019.07.26 updated to v8.2. `subsample.py` now uses `min_fl_count` cutoff in count setup. also added `subsample_with_category.py`
    
    2019.07.25 updated to v8.1. `sam_to_gff3.py` modified to work with SQANTI2 (v3.3+) changes
    
    2019.07.02 updated to v8.0. `cupcake.io.GFF.GTF` now can handle missing transcript_name field
    
    2019.06.25 updated to v7.9. fixed minor dict issue with `demux_by_barcode_groups.py`
    
    2019.06.21 updated to v7.8. fixed `demux_by_barcode_groups.py` tab/space mixing issue.
    
    2019.06.20 updated to v7.7. fixed chromosome output error in `chain_fusion_samples.py`.
    
    2019.06.07 updated to v7.6. changed preClustering to include "tucked" sequences
    
    2019.06.03 updated to v7.5. added `summarize_byloci_results.py` and `collect_all_vcf.py` in phasing/ for IsoPhase.
    
    2019.05.28 updated to v7.4. fixed `select_loci_to_phase.py` to work for short genomes.
    
    2019.05.22 updated to v7.3. made `calc_expected_accuracy_from_fastq.py` and `filter_lq_isoforms.py` work for FLNC. Also added 0-bp exon filter for collapse script.
    
    2019.04.30 updated to v7.2. fixed warning/bug in `coordinate_mapper.py` by use `str()` instead of `.tostring()` for Bio.Seq objects.
    
    2019.04.30 updated to v7.1. added `group_ORF_sequences.py` for grouping ORF predictions.
    
    2019.04.08 updated to v7.0. fixed `summarize_sample_GFF_junctions.py` for newline error.
    
    2019.03.27 updated to v6.9. fixed `clip_out_UMI_cellBC.py` to properly handle 0-length UMIs or BCs (but not both).
    
    2019.03.19 updated to v6.8. fixed `phasing.io.SAMMPileUpReader.py` for cov 0 returns
    
    2019.03.14 updated to v6.7. added `sam_to_collapsed_gff.py`
    
    2019.03.11 updated to v6.6. temp support of lazy  BED reader in BED.py
    
    2019.02.25 updated to v6.5. fixed `filter_away_subset.py` to handle edge case where the shorter one is monoexonic.
    
    2019.01.31 updated to v6.4. fixed junction 6-field support in `scrub_sample_GFF_junctions.py`. 
    
    2019.01.30 updated to v6.3. fixed typo in `summarize_sample_GFF_junctions.py`.
    
    2019.01.12 updated to v6.2. added first version of IsoPhase scripts.
    
    2018.10.29 updated to v6.1. changed confusing param name in `chain_samples.py` to `--dun-merge-5-shorter`
    
    2018.10.29 updated to v6.0. added `demux_by_barcode_group.txt` for creating demultiplexed GFF (and FASTX) from demux count files.
    
    2018.10.15 updated to v5.11. `sam_to_gff3.py` updated to allow `source` param.
    
    2018.10.12 updated to v5.10. collapse scripts further handles isoseq3 with mapping formats.
    
    2018.08.30 updated to v5.9. have collapse script handle isoseq3 formats correctly in get_fl_from_id().
    
    2018.08.01 updated to v5.8. (also tagged as `cupcake_v5.8`) fixed `sam_to_gff3.py` to output GFF3 correctly, also refreshed BioReaders.py in sequence/ to be up-to-date with cupcake/io version.
    
    2018.07.16 updated to v5.7. added `sam_to_bam.py` and `sam_to_gff3.py` (requires BCBio)
    
    2018.07.13 updated to v5.6. fixed polyA length bug in make classify report for isoseq3.
    
    2018.06.29 updated to v5.4. collapse,fusion,abundance,demux now works with isoseq3 output. 
    
    2018.03.29 updated to v5.3. Update to work with pitchfork SA5.1
    
    2018.03.12 updated to v5.2. Fixed over-collapsing genes in collapse script. Now processing strands separately in correct manner.
    
    2017.11.06 updated to v4.1. pCS merge incorrect in `chain_samples.py`. Fixed.
    
    2017.10.31 updated to v4.0. pCS merge incorrect in `run_preCluster.py`. Fixed.
    
    2017.10.10 updated to v3.9. Merged pCS branch (`--dun_use_partial`) and cdunn's random seed.
    
    2017.09.25 updated to v3.7. Fixed minor printing error in `scrubbed.group.txt` for `scrub_sample_GFF_junctions.py`.
</details>
