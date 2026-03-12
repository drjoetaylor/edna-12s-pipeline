# ============================================================
# 04_dada2_tutorial.R
# Copy-and-paste DADA2 style tutorial for this 12S pipeline
# ============================================================

# ----------------------------
# SECTION 1. Load packages
# ----------------------------

library(dada2)
library(Biostrings)
library(tibble)
library(readr)

# ----------------------------
# SECTION 2. Define paths
# ----------------------------

path <- "results/cutadapt"

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

sample.names <- sub("_R1_001.*$", "", basename(fnFs))

# Inspect file names
fnFs
fnRs
sample.names

# ----------------------------
# SECTION 3. Visualise quality
# ----------------------------

# Plot a few files to inspect quality profiles
plotQualityProfile(fnFs[1:min(6, length(fnFs))])
plotQualityProfile(fnRs[1:min(6, length(fnRs))])

# ----------------------------
# SECTION 4. Filter and trim
# ----------------------------

filt_path <- file.path("results", "filtered")
dir.create(filt_path, showWarnings = FALSE, recursive = TRUE)

filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

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

out

write.csv(out, file.path("results", "tutorial_filtering_summary.csv"))

# ----------------------------
# SECTION 5. Learn error rates
# ----------------------------

errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)

# ----------------------------
# SECTION 6. Dereplicate reads
# ----------------------------

derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)

names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Example derep object
derepFs[[1]]

# ----------------------------
# SECTION 7. Infer ASVs
# ----------------------------

dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

# Inspect one sample
dadaFs[[1]]

# ----------------------------
# SECTION 8. Merge paired reads
# ----------------------------

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)

# Inspect first sample
head(mergers[[1]])

# ----------------------------
# SECTION 9. Build sequence table
# ----------------------------

seqtab <- makeSequenceTable(mergers)

dim(seqtab)
table(nchar(colnames(seqtab)))

write.csv(seqtab, file.path("results", "tutorial_seqtab_raw.csv"))

# ----------------------------
# SECTION 10. Remove chimeras
# ----------------------------

seqtab.nochim <- removeBimeraDenovo(
  seqtab,
  method = "consensus",
  multithread = TRUE,
  verbose = TRUE
)

dim(seqtab.nochim)
sum(seqtab.nochim) / sum(seqtab)

write.csv(seqtab.nochim, file.path("results", "tutorial_seqtab_nochim.csv"))

# ----------------------------
# SECTION 11. Track reads
# ----------------------------

getN <- function(x) sum(getUniques(x))

track <- cbind(
  input = out[, 1],
  filtered = out[, 2],
  denoisedF = sapply(dadaFs, getN),
  denoisedR = sapply(dadaRs, getN),
  merged = sapply(mergers, nrow),
  nonchim = rowSums(seqtab.nochim)
)

rownames(track) <- sample.names
track

write.csv(track, file.path("results", "tutorial_read_tracking.csv"))

# ----------------------------
# SECTION 12. Export ASVs
# ----------------------------

asv_seqs <- colnames(seqtab.nochim)
asv_ids <- paste0("ASV", seq_along(asv_seqs))

dna <- DNAStringSet(asv_seqs)
names(dna) <- asv_ids

writeXStringSet(
  dna,
  filepath = file.path("results", "tutorial_asvs.nochim.fasta"),
  format = "fasta"
)

asv_lookup <- tibble(
  ASV = asv_ids,
  sequence = asv_seqs,
  length = nchar(asv_seqs)
)

write_tsv(asv_lookup, file.path("results", "tutorial_asv_lookup.tsv"))

seqtab_asv <- seqtab.nochim
colnames(seqtab_asv) <- asv_ids

write.csv(
  data.frame(sample = rownames(seqtab_asv), seqtab_asv, check.names = FALSE),
  file.path("results", "tutorial_seqtab_asv.csv"),
  row.names = FALSE
)

# ----------------------------
# SECTION 13. Save objects
# ----------------------------

saveRDS(seqtab.nochim, file.path("results", "tutorial_seqtab_nochim.rds"))
saveRDS(seqtab_asv, file.path("results", "tutorial_seqtab_asv.rds"))
saveRDS(track, file.path("results", "tutorial_read_tracking.rds"))
