#!/usr/bin/env Rscript
# Remove the detected primers from one sample (sequence-based, per read). Reverse-oriented pairs are
# swapped; pairs without both primers are dropped. With status "none" the reads are passed through.
# At the 3' end, reads that run past the amplicon (read-through) are cut before the opposite primer or
# the Illumina adapter.
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
  n_in <- 0; cut_F <- 0; cut_R <- 0
  repeat {
    a <- yield(sF); b <- yield(sR)
    if (!length(a)) break
    n_in <- n_in + length(a)
    ta <- trim_readthrough(a); tb <- trim_readthrough(b)
    cut_F <- cut_F + sum(ta$cut); cut_R <- cut_R + sum(tb$cut)
    writeFastq(ta$reads, outF, mode = "a", compress = TRUE); writeFastq(tb$reads, outR, mode = "a", compress = TRUE)
  }
  close(sF); close(sR)
  n_out <- n_in; n_swap <- 0; offs <- integer(0); left_F <- NA; left_R <- NA
} else {
  fp <- get("fw_primer"); rp <- get("rv_primer")
  sF <- FastqStreamer(opt$r1, n = 2e5); sR <- FastqStreamer(opt$r2, n = 2e5)
  n_in <- 0; n_out <- 0; n_swap <- 0; offs <- integer(0); left_F <- 0; left_R <- 0; cut_F <- 0; cut_R <- 0
  writeFastq(ShortReadQ(), outF, mode = "w", compress = TRUE)
  writeFastq(ShortReadQ(), outR, mode = "w", compress = TRUE)
  repeat {
    a <- yield(sF); b <- yield(sR)
    if (!length(a)) break
    if (length(a) != length(b)) stop("R1 and R2 of ", id, " have different numbers of reads")
    t <- trim_chunk(a, b, fp, rp)
    n_in <- n_in + length(a); n_out <- n_out + length(t$F); n_swap <- n_swap + t$n_swap
    offs <- c(offs, t$offsets)
    # check: after trimming, the read starts must no longer match the primer
    left_F <- left_F + sum(!is.na(primer_end(sread(t$F), fp)))
    left_R <- left_R + sum(!is.na(primer_end(sread(t$R), rp)))
    # 3' end: R1 runs into the reverse primer (rc) and adapter, R2 into the forward primer (rc)
    tF <- trim_readthrough(t$F, rp); tR <- trim_readthrough(t$R, fp)
    t$F <- tF$reads; t$R <- tR$reads; cut_F <- cut_F + sum(tF$cut); cut_R <- cut_R + sum(tR$cut)
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
                    spacer_max = if (length(offs)) max(0, max(offs)) else 0,
                    primer_left_R1_pct = round(100 * left_F / max(n_out, 1), 2),
                    primer_left_R2_pct = round(100 * left_R / max(n_out, 1), 2),
                    readthrough_cut_R1_pct = round(100 * cut_F / max(n_out, 1), 2),
                    readthrough_cut_R2_pct = round(100 * cut_R / max(n_out, 1), 2))
write_tsv(stats, paste0(id, ".primer_stats.tsv"))
message(sprintf("%s: %d read pairs, %d with primers (%.1f%%), %d swapped; 3' read-through cut in %.1f%% (R1) / %.1f%% (R2)",
                id, n_in, n_out, stats$with_primers_pct, n_swap, stats$readthrough_cut_R1_pct, stats$readthrough_cut_R2_pct))
