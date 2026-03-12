# =========================
# 00_config.R
# Central configuration
# =========================

# Primers
FWD <- "ACTGGGATTAGATACCCC"
REV <- "TAGAACAGGCTCCTCTAG"

# Thread count
threads <- max(1, min(10, parallel::detectCores()))

# Directories
raw_dir <- "RawSeqs"
cutadapt_dir <- file.path("results", "cutadapt")
filtered_dir <- file.path("results", "filtered")
reference_dir <- "reference_db"
unzipped_ref_dir <- file.path("results", "reference_db_unzipped")

dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create(cutadapt_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(filtered_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(unzipped_ref_dir, showWarnings = FALSE, recursive = TRUE)

# External tools
# Assumes cutadapt and usearch are already available on PATH
cutadapt_bin <- "cutadapt"
usearch_bin <- "usearch"

# SINTAX settings
sintax_cutoff <- 0.7

# Reference databases
reference_dbs <- data.frame(
  db_name = c("INBO", "MIDORI", "CLARE"),
  db_file = c(
    "all_seqs_INBO_riaz_amplified.sintax.fasta",
    "MIDORI2_LONGEST_NUC_GB269_srRNA_SINTAX.fasta.gz",
    "12s_verts.trimmed_RiazV5.sintax.fasta"
  ),
  stringsAsFactors = FALSE
)
