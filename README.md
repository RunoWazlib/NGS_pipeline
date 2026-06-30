# Introduction
The scripts contained in this directory serve as the source code for our Next-Generation Sequencing (ngs) data analysis pipeline. This analysis will follow these broad steps:
1. Optionally, trim raw reads to desired length or to bases with desired features (Not yet implemented)
2. Align data to a known reference using Bowtie2 local alignment
3. Determine nucleotide/mutational frequency, and depth per base data

These analysis files can then be piped into other downstream analysis, such as the Mathematica "Alignment Correlation.nb"

# Usage
First, clone this project. Next, in the root directory of this project, execute the following:
`bash setup.sh`
**OR**
`conda env create -f environment.yml`
`conda activate {name of environment : Default=ngs-pipeline}`
`pip install .`

This should install all scripts required for NGS data processing, alignment, and analysis.

To call on this pipeline, first generate a configuration file for your sequencing data by executing `config-builder` your data directory. This script will automatically generate configuration files (.json) for each file (or set of files in paired-end mode).

This `.json` config file controls...
1. What mode of alignment to run (i.e. Paired-end alignment or single read/merged read alignment) 
**Note:** you only need to include ONE mode, 'docs/config_template.json' includes both for the sake of example formatting

2. The path to your reference fasta file i.e. `ref.fasta`

3. The name of the directory where you want output files to be written to (any existing pipeline outputs will be overwritten). Note that the script-generated config files include the sample names as well as the custom output name, so ensure that your sample names are unique in some way if using that tool

4. The path(s) to the `.fastq` read files, which is dependent on the alignment mode (see `docs/config_template.json`)

5. Lastly, what of the complete analysis you want to generate (Note that some analysis methods require certain files to be present or generated in the output directory i.e. mpileup.txt)

Next, run `ngs-pipeline {config.json}` to analyze your data as listed in the config.

# Contributions
To contribute to this package, a user should write their additions in python 3.13 and deposit their scripts in an appropriate directory in `/src/pipeline` or `/tests`. From there, a user should add any needed flags and novel inputs to the `config_template.json` (+ a method to add them automatically in `config_builder.py`) so other users can update their config files. After this, a user should add their function to the `src/pipeline/main.py:main()` function so their addition(s) can be evaluated on subsequent runs.