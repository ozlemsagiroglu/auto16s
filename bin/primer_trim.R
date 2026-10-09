#!/usr/bin/env Rscript
# Remove the detected primers from one sample (sequence-based, per read). Reverse-oriented pairs are
# swapped; pairs without both primers are dropped. With status "none" the reads are passed through.
.script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(.script)), "utils.R"))
source(file.path(dirname(normalizePath(.script)), "primer_utils.R"))
opt <- parse_args()

pr <- read.delim(opt$primers, stringsAsFactors = FALSE)
get <- function(k) pr$value[pr$key == k]
status <- get("status")
id <- opt$sample
outF <- paste0(id, "_R1.trim.fastq.gz"); outR <- paste0(id, "_R2.trim.fastq.gz")


if (status == "none") {
  # no primers to remove: re-write compressed (input may be uncompressed .fastq)
  sF <- FastqStreamer(opt$r1, n = 2e5); sR <- FastqStreamer(opt$r2, n = 2e5)
  writeFastq(ShortReadQ(), outF, mode = "w", compress = TRUE)
  writeFastq(ShortReadQ(), outR, mode = "w", compress = TRUE)
  n_in <- 0
  repeat {
    a <- yield(sF); b <- yield(sR)
    if (!length(a)) break
    n_in <- n_in + length(a)
    writeFastq(a, outF, mode = "a", compress = TRUE); writeFastq(b, outR, mode = "a", compress = TRUE)
  }
  close(sF); close(sR)
  n_out <- n_in; n_swap <- 0; offs <- integer(0)
} else {
  fp <- get("fw_primer"); rp <- get("rv_primer")
  sF <- FastqStreamer(opt$r1, n = 2e5); sR <- FastqStreamer(opt$r2, n = 2e5)
  n_in <- 0; n_out <- 0; n_swap <- 0; offs <- integer(0)
  writeFastq(ShortReadQ(), outF, mode = "w", compress = TRUE)
  writeFastq(ShortReadQ(), outR, mode = "w", compress = TRUE)
  repeat {
    a <- yield(sF); b <- yield(sR)
    if (!length(a)) break
    if (length(a) != length(b)) stop("R1 and R2 of ", id, " have different numbers of reads")
    t <- trim_chunk(a, b, fp, rp)
    n_in <- n_in + length(a); n_out <- n_out + length(t$F); n_swap <- n_swap + t$n_swap
    offs <- c(offs, t$offsets)
    if (length(t$F)) {
      writeFastq(t$F, outF, mode = "a", compress = TRUE)
      writeFastq(t$R, outR, mode = "a", compress = TRUE)
    }
  }
  close(sF); close(sR)
}
stats <- data.frame(sample = id, input = n_in, with_primers = n_out,
                    with_primers_pct = round(100 * n_out / max(n_in, 1), 1),
                    swapped_pct = round(100 * n_swap / max(n_out, 1), 1),
                    spacer_max = if (length(offs)) max(0, max(offs)) else 0)
write_tsv(stats, paste0(id, ".primer_stats.tsv"))
message(sprintf("%s: %d read pairs, %d with primers (%.1f%%), %d swapped", id, n_in, n_out, stats$with_primers_pct, n_swap))
