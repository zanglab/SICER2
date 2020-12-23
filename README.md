# SICER2
Redesigned and improved ChIP-seq broad peak calling tool SICER

[![Build Status](https://travis-ci.com/zanglab/SICER2.svg?branch=master)](https://travis-ci.com/zanglab/SICER2)

## Introduction
Chromatin immunoprecipitation combined with high-throughput sequencing (ChIP-seq) can be used to map binding sites of a protein of interest in the genome. Histone modifications usually occupy broad chromatin domains and result in diffuse patterns in ChIP-seq data that make it difficult to identify signal enrichment. SICER, a spatial clustering approach for the identification of ChIP-enriched regions, was developed for calling broad peaks from ChIP-seq data. 

Usability of the original SICER software has been affected by increased throughputs of ChIP-seq experiments over the years. We now present SICER2, a more user-friendly version of SICER that has been redisgned and streamlined to handle large ChIP-seq data sets. This new Python package supports multiple job submissions on cluster systems and parallel processing on multicore architectures.

For more information about the original SICER algorithm, please see,

“A clustering approach for identification of enriched domains from histone modification
 ChIP-Seq data” Chongzhi Zang, Dustin E. Schones, Chen Zeng, Kairong Cui, Keji Zhao, and
 Weiqun Peng, *Bioinformatics* 25, 1952 - 1958 (2009)

In addition, we present an alternative algorithm for identification of broad domains from ChIP-seq data called RECOGNICER. It uses a coarse-graining approach to identify broad domains on both fine and coarse scale. 

## Installation
### Quick Installation
Easiest way to install SICER2 is through `pip`. Simply open the terminal and type `pip install SICER2`.

### Requirements
#### Python Version
Unlike the original version of SICER, SICER2 runs in Python 3.
Please use Python 3 to install and run SICER 2.0.

#### Libraries
Numpy and Scipy are required to run SICER2. Please have these installed before installing SICER2.
This can be done by simply typing `pip install numpy scipy` under command line (if python2.7 is your default python version, use `pip3`).

#### C Compiler
C compiler is required to compile C codes that are part of the SICER2 package. This also means that python header files (e.g. Python.h) are needed. For Linux users, make sure to have python-dev installed. For Mac OS X users, it is recommended that you install Xcode.

#### BedTools
Lastly, if you would like to directly pass BAM files as input files for SICER2, you need to have *bedtools* installed. Please refer to this [link](http://bedtools.readthedocs.io/en/latest/) for more details on installing bedtools. This is not required if you will intend to only pass BED files as input files.

### Other Installations
For local installation, the source distribution file is available at Zang Lab website ([link](http://faculty.virginia.edu/zanglab/))


## Using SICER2
The terminal command to run SICER is `sicer`. The command to run RECOGNICER is `recognicer`.

Both `sicer` command and `recognicer` command take several flag-based arguments as inputs. Two arguments are required: (1) name/path of treatment file, (2) species of reference genome. Arguments for each command are explained below.

### SICER Arguments
##### -t/--treatment_file (Required)
The file must either be in BED or BAM format (note that you need *bedtools* installed to directly enter BAM files).
The file name can either the relative path or the absolute path of the file.

##### -c/--control_file (Optional)
Like the treatment file, control file must be in BED or BAM format and can be the relative path or the absolute path of the file. However, control library input is optional.

##### -s/--species (Required)
ex) `-s hg38`

##### -rt/--redundancy_threshold (Optional)
The number of copies of indentical reads allowed in a library. Default value is 1

##### -w/--window_size (Optional)
Resolution of SICER. Default value is 200 (bp)

##### -f/--fragment_size (Optional)
The amount of shift from the beginning of a read to the center of the DNA fragment represented by the read.
Default value is 150 (bp).

##### -egf/--effective_genome_faction (Optional)
Effective genome as fraction of the genome size. Default value is 0.74.

##### -fdr/--false_discovery_rate (Optional)
Remove all islands with an false-discovery-rate below cutoff. Default value is 0.01.

##### -g/--gap_size (Optional)
The minimum length of a "gap" such that neighboring window is an "island."
Please note that this value must be a multiple of the window size.
Default value is 600 (bp).

##### -e/--e_value (Optional)
E-value. Requires user input when no control library is provided. Default value is 1000

##### -o/--output_directory (Optional)
Path of the directory in which results will be stored. Default output directory is the current working directory.

##### -cpu/--cpu (Optional)
The number of CPU cores SICER program will use when executing multi-processing tasks. Optimal number of cores is the species' number of chromosomes. Default value is the maximum number of cores avaiable in the system.

##### --significant_reads (Optional)
Significant Reads: Type "--significant_reads" flag to have SICER produce a BED file of treatment reads filtered by significant islands and WIG file of filtered reads binned into windows.

### RECOGNICER Arguments
All of the arguments for RECOGNICER are identical to those of SICER except for `gap_size` and `e_value`.
Instead of these two arguments, RECOGNICER has two arguments called `step_size` and `step_score`.

##### -s_size/--step_size (Optional)
The number of windows in one graining unit. Default value is 3.

##### -s_score/--step_score (Optional)
The minimum number of positive elements in the graining unit to call the unit positive. Default value is 2.


## Using SICER2 for differential peak calling
The commands for differential peak calling are `sicer_df` and `recognicer_df`.  

#### Arguments
The arguments for both SICER and RECOGNICER differential peak calling are identical to those of the regular peak callings except for the following arguments specified below.

Also, differential peak calling has one additional argument called `----false_discovery_rate_df`

##### -t/--treatment_file (Required)
Two files must be given as input. The first file must be the knockout (KO) file and the second file must be the wild-type (WT) file.
Both files must either be in BED or BAM format.

##### -c/--control_file (Optional)
While optional, two files must be given as input if you decide to provide the input. The first file must be the control library corresponding to the knockout (KO) treatment file and the second file must be the control library corresponding to the wild-type (WT) treatment file. Both files must either be in BED or BAM format.

##### -fdr_df/--false_discovery_rate_df (Optional)
Cutoff for identification of significant changes been wild-type library and knockout library. Default value is 0.01.


## Example Use
1. Calling SICER with a control library.
*Default parameters are explicitly entered for the sake of demonstration.*

`sicer -t treatment.bed -c control.bed -s hg38 -w 200 -rt 1 -f 150 -egf 0.74 -fdr 0.01 -g 600 -e 1000`

2. Calling SICER without a control library

`sicer -t treatment.bed -s hg38`

3. Calling SICER with control libraries for differential peak calling.

`sicer_df -t treatment1.bed treatment2.bed -c control1.bed control2.bed -s hg38`

4. Calling SICER without control libraries for differential peak calling.

`sicer_df -t treatment1.bed treatment2.bed -s hg38`

## Adding your own species
To add a new species, the user has to edit the `SICER2/sicer/lib/GenomeData.py` file directly. To do so,
1. Clone SICER2 repository.
2. Edit `sicer/lib/GenomeData.py` as following:

    Create a list of all the chromosome numbers like "hg38_chroms"
    Create a dictionary that maps chromosome numbers to their length like "hg38_chrom_lengths"
    Update both "species_chroms" and "species_chrom_lengths" dictionaries so that species name is mapped to the list of chromosomes and the length dictionary.
    Use the species name used as key value in "species_chroms" and "species_chrom_lengths" as argument for "--species"


3. Once finished with editing GenomeData.py, run  `pip install -e .` in the top directory of the repo. This should install the user's local version of SICER2.

## Contact
Please contact Zang Lab at zang@virginia.edu.
