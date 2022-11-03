# NGS_pipeline on a Plasmodium dataset

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A520.04.0-brightgreen.svg)](https://www.nextflow.io/)

## Introduction

The pipeline does Next Generation sequence analysis.

### A summary of the steps followed in our analysis include; 

- Pre-processing of the reads
  - Quality Check using `Fastqc` and `multiqc`
- Alignment of samples to the reference using `bwa`
- Variant calling using `octopus`
- Genomic coverage computation using `bedtools`

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. 

## Installation 
[Nextflow](https://www.nextflow.io/docs/latest/getstarted.html)

## Running the pipeline
The pipeline does not require installation as `NextFlow` will automatically fetch it from `GitHub`.

### Test data
To execute the pipeline, run:
   
    nextflow run nanjalaruth/NGS_pipeline -profile slurm -resume --reads "*_{1,2}.fastq.gz" --reference "path to the reference genome"
    
## To run the updated version of this pipeline, run:

    nextflow pull nanjalaruth/NGS_pipeline
    
## Arguments

### Required Arguments
| Argument  | Usage                            | Description                                                          |
|-----------|----------------------------------|----------------------------------------------------------------------|
| -profile  | \slurm\>                    | Configuration profile to use.                                       |
| --input  | \</project/\*\_{1,2}\*.fastq.gz\> | Directory pattern for fastq files.                                   |
| --reference_genome    | \<hg19\>              | Path to the reference genome to which the samples will be mapped |
| -r    | \<revision\>  | Pipeline revision     |
| --singleEnd  |             | Specifies that the input files are not paired reads (default is paired-end). |

## Support
I track open tasks using github's [issues](https://github.com/nanjalaruth/NGS_pipeline/issues)
