# =========================
# 02_dada2_asv.R
# DADA2 pipeline and ASV export
# =========================

source("scripts/00_config.R")

library(dada2)
library(Biostrings)
library(tibble)
library(readr)

path <- cutadapt_dir

fnFs <- sort(list.files(
  path,
  pattern = "_R1_001.*\\.(fastq|fq)(\\.gz)?$",
  full.names = TRUE
))
fnRs <- sort(list.files(
  path,
  pattern = "_R2_001.*\\.(fastq|fq)(\\.gz)?$",
  full.names = TRUE
))

if (length(fnFs) == 0 || length(fnRs) == 0) {
  stop("No cutadapt output files found in results/cutadapt/")
}

sample.names <- sub("_R1_001.*$", "", basename(fnFs))

filtFs <- file.path(filtered_dir, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filtered_dir, paste0(sample.names, "_R_filt.fastq.gz"))

cat("Filtering and trimming reads\n")
out <- filterAndTrim(
  fnFs, filtFs,
  fnRs, filtRs,
  truncLen = c(75, 75),
  maxN = 0,
  maxEE = c(2, 2),
  truncQ = 2,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = TRUE
)

write.csv(out, file.path("results", "filtering_summary.csv"))

cat("Learning error rates\n")
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

cat("Dereplicating reads\n")
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)

names(derepFs) <- sample.names
names(derepRs) <- sample.names

cat("Running dada\n")
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

cat("Merging pairs\n")
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)

cat("Making sequence table\n")
seqtab <- makeSequenceTable(mergers)
write.csv(seqtab, file.path("results", "seqtab_raw.csv"))

cat("Removing chimeras\n")
seqtab.nochim <- removeBimeraDenovo(
  seqtab,
  method = "consensus",
  multithread = TRUE,
  verbose = TRUE
)

write.csv(seqtab.nochim, file.path("results", "seqtab_nochim.csv"))

cat("Writing ASV fasta and lookup\n")
asv_seqs <- colnames(seqtab.nochim)
asv_ids <- paste0("ASV", seq_along(asv_seqs))

dna <- DNAStringSet(asv_seqs)
names(dna) <- asv_ids
writeXStringSet(
  dna,
  filepath = file.path("results", "asvs.nochim.fasta"),
  format = "fasta"
)

asv_lookup <- tibble(
  ASV = asv_ids,
  sequence = asv_seqs,
  length = nchar(asv_seqs)
)
write_tsv(asv_lookup, file.path("results", "asv_lookup.tsv"))

seqtab_asv <- seqtab.nochim
colnames(seqtab_asv) <- asv_ids

write.csv(
  data.frame(sample = rownames(seqtab_asv), seqtab_asv, check.names = FALSE),
  file.path("results", "seqtab_asv.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(ASV = colnames(seqtab_asv), t(seqtab_asv), check.names = FALSE),
  file.path("results", "seqtab_asv_transposed.csv"),
  row.names = FALSE
)

saveRDS(seqtab_asv, file.path("results", "seqtab_asv.rds"))

cat("DADA2 complete\n")
