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
---

# Blank correction and contamination cleanup

After taxonomy assignment with the CLARE reference (`12s_verts.trimmed_RiazV5.sintax.fasta`), the repository includes an additional workflow to identify and correct potential contamination using blank controls.

This workflow:

- keeps fish taxa (Actinopteri and Hyperoartia are not removed)
- separates **PCR blanks**, **extraction blanks**, and **site blanks**
- calculates blank-based limits of detection (LOD)
- generates cleaned ASV matrices using:
  - laboratory blanks only
  - site blanks only
  - both blank types
- extracts the **lowest accepted taxonomic level and confidence**
- produces diagnostic tables to evaluate spillover into site blanks

---

# Workflow overview

The cleanup stage occurs after taxonomy assignment:

```
FASTQ
  ↓
cutadapt
  ↓
DADA2 denoising
  ↓
ASV table
  ↓
SINTAX taxonomy (CLARE database)
  ↓
Create cleanup input matrix
  ↓
Blank contamination correction
  ↓
Cleaned ASV tables
```

---

# Required input

These scripts expect the workflow output:

```
results/asv_taxonomy_abundance_CLARE.csv
```

This file is produced by:

```
scripts/03_sintax_assign.R
```

---

# Step 1 — Create cleanup input matrix

Run:

```r
source("scripts/05_make_cleanup_input_from_CLARE.R")
```

This converts the ASV-level CLARE output into a taxon × sample matrix suitable for blank filtering.

Output:

```
results/ncl_matrix_raw.csv
```

---

# Step 2 — Run blank contamination cleanup

Run:

```r
source("scripts/06_blank_cleanup_from_workflow.R")
```

This script:

- identifies blank types
- calculates blank-based LOD thresholds
- applies contamination filtering
- generates diagnostic summaries

---

# Output files

Key outputs are written to:

```
results/
```

### Cleaned ASV matrices

| File | Description |
|-----|-------------|
| `ncl_cleaned_labLOD.csv` | contamination corrected using PCR + extraction blanks |
| `ncl_cleaned_siteLOD.csv` | contamination corrected using site blanks |
| `ncl_cleaned_bothLOD.csv` | contamination corrected using both blank types |

---

### Presence / absence matrices

```
ncl_cleaned_labLOD_pa.csv
ncl_cleaned_siteLOD_pa.csv
ncl_cleaned_bothLOD_pa.csv
```

---

### Diagnostic tables

| File | Purpose |
|----|----|
| `cleanup_blank_diagnostic_table.csv` | evaluates blank spillover patterns |
| `cleanup_sample_classification.csv` | shows blank/sample classification |
| `cleanup_summary.csv` | summary of read counts before and after filtering |

---

### Long format data

```
ncl_cleaned_long.csv
```

Contains:

- ASV / taxon
- taxonomy
- lowest taxonomic level
- confidence score
- sample name
- site
- raw reads
- cleaned reads

---

# Recommended interpretation workflow

Start by comparing:

```
ncl_cleaned_labLOD.csv
ncl_cleaned_bothLOD.csv
```

and inspect:

```
cleanup_blank_diagnostic_table.csv
```

This helps determine whether **site blanks represent local spillover** or true contamination before deciding how aggressively to filter.

---

# Copy-paste tutorial

### Run the full workflow

```r
source("scripts/run_pipeline.R")
```

---

### Run cleanup only

```r
source("scripts/05_make_cleanup_input_from_CLARE.R")
source("scripts/06_blank_cleanup_from_workflow.R")
```

---

# Notes

Large files are intentionally excluded from the repository:

- FASTQ files
- reference databases
- results

These should be added locally before running the pipeline.

---

## Notes

This repository does not include:

- raw sequencing files
- reference databases
- usearch binary
- results

These must be added locally before running the pipeline.
