# 12S eDNA metabarcoding pipeline

This repository contains a 12S metabarcoding workflow using:

- cutadapt
- DADA2
- USEARCH SINTAX

It supports two usage styles:

1. full pipeline execution
2. DADA2-style stepwise tutorial

## Repository structure

```text
edna-12s-pipeline/
├── README.md
├── .gitignore
├── scripts/
├── RawSeqs/
├── reference_db/
└── results/
```

## Required software

- R
- cutadapt
- usearch

## Required R packages

Run this in R:

```r
install.packages(c("dplyr", "readr", "stringr", "tibble"))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c("dada2", "Biostrings"))
```

## Input files

Put paired-end reads in:

```text
RawSeqs/
```

Example:

```text
Sample1_R1_001.fastq.gz
Sample1_R2_001.fastq.gz
Sample2_R1_001.fastq.gz
Sample2_R2_001.fastq.gz
```

## Reference databases

Put these into:

```text
reference_db/
```

Expected files:

```text
all_seqs_INBO_riaz_amplified.sintax.fasta
MIDORI2_LONGEST_NUC_GB269_srRNA_SINTAX.fasta.gz
12s_verts.trimmed_RiazV5.sintax.fasta
```

The MIDORI file can remain gzipped.

---

## Copy and paste tutorial: run full pipeline in R

Open R in the project directory and run:

```r
source("scripts/run_pipeline.R")
```

---

## Copy and paste tutorial: run full pipeline on Linux server

From the project directory:

```bash
Rscript scripts/run_pipeline.R
```

To save a log:

```bash
Rscript scripts/run_pipeline.R > pipeline.log 2>&1
```

To monitor progress:

```bash
tail -f pipeline.log
```

---

## Copy and paste tutorial: DADA2-style stepwise analysis in R

### Step 1. Run cutadapt

```r
source("scripts/01_cutadapt.R")
```

### Step 2. Run DADA2 ASV inference

```r
source("scripts/02_dada2_asv.R")
```

### Step 3. Run taxonomy assignment

```r
source("scripts/03_sintax_assign.R")
```

### Step 4. Run the tutorial version

```r
source("scripts/04_dada2_tutorial.R")
```

This tutorial script is written in a stepwise style similar to the DADA2 tutorial, with separate sections for:

- reading files
- plotting quality
- filtering
- learning errors
- dereplication
- ASV inference
- merging
- chimera removal
- read tracking
- exporting ASVs

---

## Expected outputs

Results will appear in:

```text
results/
```

Main outputs include:

```text
asvs.nochim.fasta
asv_lookup.tsv
seqtab_asv.csv
asv_taxonomy_INBO.tsv
asv_taxonomy_MIDORI.tsv
asv_taxonomy_CLARE.tsv
species_abundance_INBO.csv
species_abundance_MIDORI.csv
species_abundance_CLARE.csv
```

---

## Notes

This repository does not include:

- raw sequencing files
- reference databases
- usearch binary
- results

These must be added locally before running the pipeline.
