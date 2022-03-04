# Bulk RNAseq processing pipeline
Bulk RNAseq workflow for Biowulf based on STAR read mapper.

## To build conda environment and install required tools run the following command:
```
mamba env create -f=environment_bbtools.yml -n rnaseq
```
## To run the entire snakemake pipeline
```
conda activate rnaseq
snakemake --profile slurm Snakefile
```

## To run a specific rule within the pipeline
First, try a dry run to check that everything works as expected
```
snakemake -R --until $MY_RULE --cores $CPUs -n
```
Then run the pipeline with the command below:
```
snakemake -R --until $MY_RULE --cores $CPUs
```
