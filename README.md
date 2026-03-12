# edna-12s-pipeline
Trim primers, run dada2, classify using Sintax 
# 12S eDNA metabarcoding pipeline

This repository contains scripts to process 12S amplicon data using cutadapt, DADA2, and USEARCH SINTAX.

## Folder structure

- `scripts/` = pipeline scripts
- `RawSeqs/` = input FASTQ files
- `reference_db/` = reference databases for SINTAX
- `results/` = output files

## Requirements

- R
- cutadapt
- USEARCH v11

## Before running

1. Put paired FASTQ files in `RawSeqs/`
2. Put SINTAX databases in `reference_db/`
3. Update paths in `scripts/00_config.R` if needed

## Run

```r
source("scripts/run_pipeline.R")
