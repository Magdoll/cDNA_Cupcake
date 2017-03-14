# cDNA_Cupcake

![cupcake](http://i87.photobucket.com/albums/k138/Elizabeth_Tseng/cupcake-778703_640_zpsonkeyb9w.png)

Last Updated: 03/14/2017

**cDNA_Cupcake** is a miscellaneous collection of Python and R scripts used for analyzing sequencing data. Most of the scripts only require [Biopython](http://biopython.org/wiki/Download). For scripts that require additional libraries, it will be specified in documentation.

Current version: 2.4

## Python Requirements
* Python >= 2.7
* Biopython 

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
* `parse_matchAnnot.py`: Parse matchAnnot results into summary format.
* `make_file_for_subsampling_from_collapsed.py`: Prepare file for running subsampling (rarefaction curve).
* `subsample.py`: Running subsamping. Results can be plotted with Excel graphics and R, etc.

### Sequence Manipulation
* `get_seq_stats.py`: Summarize length distribution of a FASTA/FASTQ file.
* `rev_comp.py`: Reverse complement a sequence from command line.
* `fa2fq.py` and `fq2fa.py`: Convert between FASTA and FASTQ format.
* `sort_fasta_by_len.py`: sort fasta file by length (increasing or decreasing).
* `get_seqs_from_list.py`: extract list of sequences given a fasta file and a list of IDs.
* `err_correct_w_genome.py`: generate fasta sequences given genom

### Sequence Simulation
* `simulate.py`: Simulate error in sequences.

### Cupcake ToFU: supporting scripts for Iso Seq after clustering step
* `collapse_isoforms_by_sam.py`: Collapse HQ isoform results to unique isoforms (based on genome alignment).
* `get_abundance_post_collapse.py`: Obtain count information post collapse to unique isoforms.
* `filter_by_count.py`: Filter collapse result by FL count information.
* `filter_away_subset.py`: Filter away 5' degraded isoforms.
* `chain_samples.py`: Chaining together multiple samples
