#!/usr/bin/env Rscript
# Builds a small realistic test set from DADA2's bundled V4 example reads (2 x 250 bp, mouse gut).
# Real 515F/806R primers (degenerate bases resolved at random per read) are prepended, so trimLeft
# and primer detection are exercised. 8 samples, 2 groups with different mixing ratios.
#   Rscript tests/make_test_data.R tests/data                  # primers at a fixed position
#   Rscript tests/make_test_data.R tests/data_spacer spacer    # 0-7 nt spacers + half of the pairs reverse-oriented
#   Rscript tests/make_test_data.R tests/data_noprimer none    # primers already removed
suppressPackageStartupMessages(library(ShortRead))
args <- commandArgs(trailingOnly = TRUE)
out  <- if (length(args) >= 1) args[1] else "tests/data"
mode <- if (length(args) >= 2) args[2] else "fixed"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
set.seed(1)
ex <- function(f) system.file("extdata", f, package = "dada2")
s1F <- readFastq(ex("sam1F.fastq.gz")); s1R <- readFastq(ex("sam1R.fastq.gz"))
s2F <- readFastq(ex("sam2F.fastq.gz")); s2R <- readFastq(ex("sam2R.fastq.gz"))

iupac <- list(A="A",C="C",G="G",T="T",R=c("A","G"),Y=c("C","T"),M=c("A","C"),K=c("G","T"),
              S=c("C","G"),W=c("A","T"),H=c("A","C","T"),V=c("A","C","G"),N=c("A","C","G","T"))
resolve <- function(p, n) vapply(seq_len(n), function(i)
  paste(vapply(strsplit(p, "")[[1]], function(b) sample(iupac[[b]], 1), ""), collapse = ""), "")
FWD <- "GTGYCAGCMGCCGCGGTAA"     # 515F, 19 nt
REV <- "GGACTACNVGGGTWTCTAAT"    # 806R, 20 nt

add_primer <- function(fq, primer, spacer = integer(length(fq))) {
  if (mode == "none") return(fq)
  n <- length(fq)
  sp <- vapply(spacer, function(k) paste(sample(c("A","C","G","T"), k, TRUE), collapse = ""), "")
  p <- paste0(sp, resolve(primer, n))
  q <- substr(as.character(quality(quality(fq))), 1, nchar(p))         # reuse the read's own first qualities
  ShortReadQ(DNAStringSet(paste0(p, as.character(sread(fq)))),
             FastqQuality(paste0(q, as.character(quality(quality(fq))))), id(fq))
}
mix <- function(frac1, n = 2500) {
  k1 <- rbinom(1, n, frac1)
  i1 <- sample(length(s1F), k1, replace = TRUE); i2 <- sample(length(s2F), n - k1, replace = TRUE)
  list(F = append(s1F[i1], s2F[i2]), R = append(s1R[i1], s2R[i2]))
}
design <- data.frame(sample = sprintf("S%02d", 1:8), group = rep(c("control", "treatment"), each = 4),
                     frac = c(0.85, 0.80, 0.90, 0.75, 0.20, 0.15, 0.25, 0.10))
for (i in seq_len(nrow(design))) {
  m <- mix(design$frac[i])
  n <- length(m$F)
  spF <- if (mode == "spacer") sample(0:7, n, TRUE) else integer(n)
  spR <- if (mode == "spacer") sample(0:7, n, TRUE) else integer(n)
  F <- add_primer(m$F, FWD, spF); R <- add_primer(m$R, REV, spR)
  if (mode == "spacer") {                       # ligation-type libraries: half of the pairs come out swapped
    sw <- seq_len(n) %% 2 == 0
    F2 <- append(F[!sw], R[sw]); R2 <- append(R[!sw], F[sw]); F <- F2; R <- R2
  }
  writeFastq(F, file.path(out, paste0(design$sample[i], "_R1.fastq.gz")), compress = TRUE)
  writeFastq(R, file.path(out, paste0(design$sample[i], "_R2.fastq.gz")), compress = TRUE)
}
file.copy(ex("example_train_set.fa.gz"), file.path(out, "example_train_set.fa.gz"), overwrite = TRUE)
write.csv(data.frame(sample = design$sample,
                     fastq_1 = paste0(design$sample, "_R1.fastq.gz"),
                     fastq_2 = paste0(design$sample, "_R2.fastq.gz"),
                     group = design$group), file.path(out, "samplesheet.csv"), row.names = FALSE, quote = FALSE)
cat("Test data written to", out, "\n")
