Nm_Detector: version 1
# Nm-seq read processing pipeline

Scripts and pipelines provided in this repository aid to process Nm-seq read libraries. It contains all scripts to allow a self-assembled processing and additionally provides pipeline scripts that run the entire processing automatically.

# Requirements

To run this pipeline, your computer requires **40 GB of available memory (RAM)** to process larger genomes (e.g. human or mouse). Moreover, snakemake was used to facilitate the automated execution of all analysis steps. The easiest way to make use of the pipeline is to set up a python3 virtual environment and run the pipeline is this environment. 
Download/Provide all necessary files:

hg38.genome: http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
hg38.gtf: http://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.gtf.gz
chromsome size: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
bedGraphToBigwig: https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64.v385/bedGraphToBigWig
R packages: DESeq2, dplyr, tidyr

# snakemake
Snakemake is a workflow management system that helps to create and execute data processing pipelines. It requires python3 and can be most easily installed via the bioconda package of the python anaconda distribution.

conda create -n Nmseq-snakemake -c bioconda -c conda-forge --file requirements.txt python=3

# Activate the environment
  ```bash
  source activate Nmseq-snakemake
  ```
To exit the environment (after finishing the usage of the pipeline), just execute
  ```bash
  source deactivate
  ```
# Run the pipeline

# Configure input parameters

The working directory contains a file named `Nmseq_config.yaml`. It's the central file in which all user settings, paramter values and path specifications are stored. During a run, all steps of the pipeline will retrieve their paramter values from this file. It follows the yaml syntax (find more information about yaml and it's syntax [here](http://www.yaml.org/)) what makes it easy to read and edit. The main principles are:
  - everything that comes after a `#` symbol is considered as comment and will not be interpreted
  - paramters are given as key-value pair, with `key` being the name and `value` the value of any paramter

Before starting the pipeline, open the `Nmseq_config.yaml` configuration file and set all options according as required. This should at least include:
  - **name of the input directory** - where are your input fastq files stored
  - **name of the output directory** - where should the pipeline store the output files (the direcotry is created if not existing)
  - **name of the log directory** - where should the pipeline store the log files
  - **name(s) of your input samples** - please note: If your sample is named `sample1.fq.gz` then `sample1` will be kept as naming scheme throughout the entire run to indicate output files that belong to this input file, e.g. the pipeline will create a file called `sample1.3pSites.noIP.bed.gz`. If you have multiple input files, just follow the given pattern with one sample name per line (and a dash that indicates another list item).
 - **groups of your input samples**, must contain "type","control","treat", makesure "-" in your sample name replaced with "."


# Start a run

Once you set up your configuration file, running the pipeline locally on your computer is as easy as invoking:
`sh step.sh`




