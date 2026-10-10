#!/usr/bin/env Rscript
# DADA2 paired-end workflow
#   Input: reads with primers already removed (PRIMER_TRIM)
#   1. diagnostics on a subset of reads: quality per cycle, amplicon length
#   2. choose truncLen (quality drop + guaranteed R1/R2 overlap) unless given
#   3. filterAndTrim (truncLen, maxEE) -> learnErrors F/R -> dada F/R
#      -> mergePairs -> removeBimeraDenovo
.script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(.script)), "utils.R"))
opt <- parse_args()
suppressPackageStartupMessages({ library(dada2); library(ggplot2) })

cpus        <- as.integer(opt$cpus)
# Second guard: memory that is actually free now (other programs, WSL limits). ~2 GB per thread.
mem_avail <- tryCatch({
  l <- grep("^MemAvailable:", readLines("/proc/meminfo"), value = TRUE)
  as.numeric(gsub("[^0-9]", "", l)) / 1024^2
}, error = function(e) NA, warning = function(w) NA)
if (length(mem_avail) == 1 && !is.na(mem_avail)) cpus <- max(1L, min(cpus, as.integer(mem_avail %/% 2)))
message(sprintf("Threads: %d (free memory: %s GB)", cpus, ifelse(is.na(mem_avail), "unknown", sprintf("%.1f", mem_avail))))
fw_len      <- 0L   # primers are removed before this step
rv_len      <- 0L
max_ee      <- as.numeric(opt$max_ee)
qmin        <- as.numeric(opt$trunc_qmin)
min_overlap <- as.integer(opt$min_overlap)
set.seed(as.integer(opt$seed))

man <- read.delim(opt$manifest, stringsAsFactors = FALSE)
# samples without any read left after primer removal cannot be processed
empty_in <- vapply(man$r1, function(f) { s <- ShortRead::FastqStreamer(f, n = 1); on.exit(close(s)); length(ShortRead::yield(s)) == 0 }, TRUE)
if (any(empty_in)) message("WARNING: no reads left after primer removal, sample dropped: ", paste(man$sample[empty_in], collapse = ", "))
dropped_primer <- man$sample[empty_in]
man <- man[!empty_in, , drop = FALSE]
if (!nrow(man)) stop("No sample has reads with the detected primers.")
fwd <- setNames(man$r1, man$sample)
rev <- setNames(man$r2, man$sample)
message("Samples: ", nrow(man))

# ======================================================================
# 1. Diagnostics on the first reads of every sample (R1/R2 stay paired)
# ======================================================================
head_reads <- function(f, n) { s <- ShortRead::FastqStreamer(f, n = n); on.exit(close(s)); ShortRead::yield(s) }
n_each <- max(500, ceiling(30000 / nrow(man)))
fqF <- lapply(fwd, head_reads, n = n_each)
fqR <- lapply(rev, head_reads, n = n_each)

qual_matrix <- function(fqs) {
  ms <- lapply(fqs, function(fq) as(Biostrings::quality(fq), "matrix"))
  w  <- max(vapply(ms, ncol, 1L))
  do.call(rbind, lapply(ms, function(m) if (ncol(m) < w) cbind(m, matrix(NA, nrow(m), w - ncol(m))) else m))
}
qF <- qual_matrix(fqF); qR <- qual_matrix(fqR)
seqF <- unlist(lapply(fqF, function(x) as.character(ShortRead::sread(x))), use.names = FALSE)
seqR <- unlist(lapply(fqR, function(x) as.character(ShortRead::sread(x))), use.names = FALSE)
lenF <- nchar(seqF); lenR <- nchar(seqR)

cycle_stats <- function(m, read) {
  data.frame(read = read, cycle = seq_len(ncol(m)),
             median = apply(m, 2, median, na.rm = TRUE),
             q25 = apply(m, 2, quantile, 0.25, na.rm = TRUE),
             q75 = apply(m, 2, quantile, 0.75, na.rm = TRUE))
}
cs <- rbind(cycle_stats(qF, "R1 (forward)"), cycle_stats(qR, "R2 (reverse)"))

primer_msgs <- character(0)

# --- amplicon length (without primers) from the overlap of R1 and reverse-complemented R2
rc <- function(x) as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(x)))
# only pairs without ambiguous bases (N): nwalign accepts A/C/G/T only
acgt <- which(grepl("^[ACGT]+$", seqF) & grepl("^[ACGT]+$", seqR))
idx <- if (length(acgt)) acgt[unique(round(seq(1, length(acgt), length.out = min(400, length(acgt)))))] else integer(0)
amp_len <- vapply(idx, function(i) {
  a <- substring(seqF[i], fw_len + 1)
  b <- rc(substring(seqR[i], rv_len + 1))
  al <- dada2::nwalign(a, b, endsfree = TRUE)
  s1 <- strsplit(al[1], "")[[1]]; s2 <- strsplit(al[2], "")[[1]]
  both <- s1 != "-" & s2 != "-"
  if (sum(both) < 20 || mean(s1[both] != s2[both]) > 0.25) return(NA_real_)
  max(which(s2 != "-")) - which(s1 != "-")[1] + 1
}, numeric(1))
amp <- if (sum(!is.na(amp_len)) >= 20) median(amp_len, na.rm = TRUE) else NA
message(sprintf("Estimated amplicon length without primers: %s (from %d of %d read pairs)",
                ifelse(is.na(amp), "not determined", amp), sum(!is.na(amp_len)), length(idx)))

# ======================================================================
# 2. truncLen
# ======================================================================
cap <- function(len) as.integer(floor(quantile(len, 0.10)))   # reads shorter than truncLen are discarded
capF <- cap(lenF); capR <- cap(lenR)
qual_cut <- function(med, start, capv) {
  bad <- which(med < qmin & seq_along(med) > start + 20)
  min(if (length(bad)) bad[1] - 1 else length(med), capv)
}
trunc_note <- character(0)
exp_pass <- NA
if (!is.na(opt$trunc_len_f) && !is.na(opt$trunc_len_r)) {
  tF <- as.integer(opt$trunc_len_f); tR <- as.integer(opt$trunc_len_r)
  trunc_note <- "set by user"
} else {
  # Per read: the longest truncation at which it still passes filterAndTrim
  # (expected errors <= maxEE, no base with Q <= truncQ before it, read long enough).
  pass_len <- function(m, sq) {
    m[is.na(m)] <- 0
    npos <- gregexpr("N", sq, fixed = TRUE)                   # maxN = 0: a read fails at its first N
    for (i in which(vapply(npos, function(x) x[1] > 0, TRUE))) m[i, npos[[i]][npos[[i]] <= ncol(m)]] <- 0
    ee <- t(apply(10^(-m / 10), 1, cumsum))
    ok <- ee <= max_ee & m > 2
    apply(ok, 1, function(r) { b <- which(!r); if (length(b)) b[1] - 1L else length(r) })
  }
  LF <- pass_len(qF, seqF); LR <- pass_len(qR, seqR)
  amp_hi <- if (is.na(amp)) NA else as.numeric(quantile(amp_len, 0.99, na.rm = TRUE))   # longest common amplicon
  if (!is.na(amp_hi)) {
    need <- ceiling(amp_hi) + min_overlap + 8          # 8 bp safety margin
    wF <- ncol(qF); wR <- ncol(qR)
    # kept[a, b] = share of read pairs passing with truncLen = (a, b): 2-D reverse cumulative count
    tab <- matrix(0, wF + 1, wR + 1)
    for (k in seq_along(LF)) tab[LF[k] + 1, LR[k] + 1] <- tab[LF[k] + 1, LR[k] + 1] + 1
    kept <- apply(apply(tab[(wF + 1):1, (wR + 1):1, drop = FALSE], 2, cumsum), 1, cumsum)
    kept <- t(kept)[(wF + 1):1, (wR + 1):1, drop = FALSE] / length(LF)     # row a+1 = truncLen a
    grid <- expand.grid(tF = 50:wF, tR = 50:wR)
    grid <- grid[grid$tF + grid$tR >= need, ]
    if (nrow(grid)) {
      grid$kept <- kept[cbind(grid$tF + 1, grid$tR + 1)]
      best <- max(grid$kept)
      # within 1 percentage point of the best: prefer the longest reads (more overlap, more sequence for denoising)
      cand <- grid[grid$kept >= best - 0.01, ]
      cand <- cand[order(-(cand$tF + cand$tR), -cand$kept), ][1, ]
      tF <- cand$tF; tR <- cand$tR
      old <- c(qual_cut(cs$median[cs$read == "R1 (forward)"], fw_len, capF), qual_cut(cs$median[cs$read == "R2 (reverse)"], rv_len, capR))
      old_kept <- if (sum(old) >= need) kept[old[1] + 1, old[2] + 1] else NA
      trunc_note <- sprintf("chosen to keep the most read pairs through maxEE = %g with >= %d bp overlap: %.1f%% expected to pass%s",
                            max_ee, min_overlap + 8, 100 * cand$kept,
                            if (is.na(old_kept)) "" else sprintf(" (median-quality rule %d/%d: %.1f%%)", old[1], old[2], 100 * old_kept))
      exp_pass <- cand$kept
    } else {
      tF <- wF; tR <- wR
      trunc_note <- "reads cannot cover the amplicon with the required overlap even untruncated"
      message("WARNING: R1+R2 are too short to overlap over the amplicon; many pairs will not merge.")
    }
  } else {
    tF <- qual_cut(cs$median[cs$read == "R1 (forward)"], fw_len, capF)
    tR <- qual_cut(cs$median[cs$read == "R2 (reverse)"], rv_len, capR)
    trunc_note <- c(sprintf("quality: first cycle with median Q < %g", qmin), "amplicon length unknown: overlap not checked")
  }
}
exp_overlap <- if (is.na(amp)) NA else (tF - fw_len) + (tR - rv_len) - amp
message(sprintf("truncLen = c(%d, %d); expected overlap = %s bp", tF, tR, exp_overlap))

trunc_report <- data.frame(
  parameter = c("truncLen R1", "truncLen R2", "estimated amplicon length (without primers)",
                "expected overlap (bp)", "read pairs expected to pass filtering", "truncLen decision"),
  value = c(tF, tR, ifelse(is.na(amp), "NA", amp), ifelse(is.na(exp_overlap), "NA", exp_overlap),
            ifelse(is.na(exp_pass), "NA", sprintf("%.1f%%", 100 * exp_pass)),
            paste(trunc_note, collapse = "; ")))

# --- diagnostic plots
vl <- data.frame(read = c("R1 (forward)", "R2 (reverse)"), trim = c(fw_len, rv_len), trunc = c(tF, tR))
p <- ggplot(cs, aes(cycle)) +
  geom_ribbon(aes(ymin = q25, ymax = q75), fill = "grey80") +
  geom_line(aes(y = median), linewidth = 0.7) +
  geom_hline(yintercept = qmin, linetype = 3, colour = "grey40") +
  geom_vline(data = vl, aes(xintercept = trunc), colour = "#0072B2", linewidth = 0.8) +
  facet_wrap(~ read, ncol = 1) + coord_cartesian(ylim = c(0, 42)) +
  labs(title = "Read quality and truncation (after primer removal)",
       subtitle = sprintf("blue = truncLen (%d / %d) | amplicon %s bp, R1/R2 overlap %s bp",
                          tF, tR, ifelse(is.na(amp), "?", amp), ifelse(is.na(exp_overlap), "?", exp_overlap)),
       x = "Cycle (base position)", y = "Quality (median, IQR)") + theme_amp()
save_fig(p, "dada2_quality_truncation", 9, 6.5, mqc = TRUE)



# ======================================================================
# 3. DADA2
# ======================================================================
dir.create("filt", showWarnings = FALSE)
filtF <- setNames(file.path("filt", paste0(man$sample, "_F_filt.fastq.gz")), man$sample)
filtR <- setNames(file.path("filt", paste0(man$sample, "_R_filt.fastq.gz")), man$sample)
out <- filterAndTrim(fwd, filtF, rev, filtR, truncLen = c(tF, tR),
                     maxN = 0, maxEE = c(max_ee, max_ee), truncQ = 2, rm.phix = TRUE,
                     compress = TRUE, multithread = cpus)
rownames(out) <- man$sample
keep <- out[, "reads.out"] > 0 & file.exists(filtF) & file.exists(filtR)
if (any(!keep)) message("WARNING: no reads left after filtering, sample dropped: ", paste(man$sample[!keep], collapse = ", "))
if (sum(keep) == 0) stop("No reads passed filtering. Check primer detection and truncLen in the report.")
fF <- filtF[keep]; fR <- filtR[keep]

errF <- learnErrors(fF, multithread = cpus, randomize = TRUE)
errR <- learnErrors(fR, multithread = cpus, randomize = TRUE)
save_fig(plotErrors(errF, nominalQ = TRUE) + ggtitle("Error model R1"), "dada2_errors_R1", 9, 7)
save_fig(plotErrors(errR, nominalQ = TRUE) + ggtitle("Error model R2"), "dada2_errors_R2", 9, 7)

wrap <- function(x, cls) if (inherits(x, cls)) setNames(list(x), names(fF)) else x
dadaF <- wrap(dada(fF, err = errF, multithread = cpus), "dada"); names(dadaF) <- names(fF)
dadaR <- wrap(dada(fR, err = errR, multithread = cpus), "dada"); names(dadaR) <- names(fR)
mergers <- mergePairs(dadaF, fF, dadaR, fR, minOverlap = min_overlap, trimOverhang = TRUE)
if (is.data.frame(mergers)) mergers <- setNames(list(mergers), names(fF))

seqtab <- makeSequenceTable(mergers)
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = cpus, verbose = TRUE)
saveRDS(seqtab.nochim, "seqtab_nochim.rds")

# ======================================================================
# 4. Read tracking
# ======================================================================
getN <- function(x) sum(getUniques(x))
na_vec <- function() setNames(rep(NA_real_, nrow(man)), man$sample)
dF <- na_vec(); dR <- na_vec(); mg <- na_vec(); nc <- na_vec()
dF[names(dadaF)] <- vapply(dadaF, getN, 1); dR[names(dadaR)] <- vapply(dadaR, getN, 1)
mg[names(mergers)] <- vapply(mergers, function(m) sum(m$abundance), 1)
nc[rownames(seqtab.nochim)] <- rowSums(seqtab.nochim)
ps_files <- list.files(".", pattern = "\\.primer_stats\\.tsv$")
pstat <- if (length(ps_files)) do.call(rbind, lapply(ps_files, read.delim)) else NULL
raw_in <- if (!is.null(pstat)) pstat$input[match(man$sample, pstat$sample)] else out[, "reads.in"]
track <- data.frame(sample = man$sample, input = raw_in, primers_found = out[, "reads.in"], filtered = out[, "reads.out"],
                    denoisedF = dF, denoisedR = dR, merged = mg, nonchim = nc)
track$merge_rate_pct <- round(100 * track$merged / pmin(track$denoisedF, track$denoisedR), 1)
track$retained_pct   <- round(100 * track$nonchim / track$input, 1)
write_tsv(track, "dada2_read_tracking.tsv")
write_mqc_table(track, "dada2_read_tracking_mqc.tsv", "dada2_read_tracking", "DADA2 read tracking",
                "Read pairs per sample through each DADA2 step.")

med_merge <- median(track$merge_rate_pct, na.rm = TRUE)
chim_kept <- sum(seqtab.nochim) / sum(seqtab)
message(sprintf("Median merge rate: %.1f%% | reads kept after chimera removal: %.1f%%", med_merge, 100 * chim_kept))
if (!is.na(med_merge) && med_merge < 70) message("WARNING: median merge rate < 70%. truncLen may be too short for the overlap; see dada2_quality_truncation.png")
if (chim_kept < 0.8) message("WARNING: < 80% of reads survived chimera removal. Most common cause: primers still in the reads (see primer detection).")

# ---------- summary table with all warnings (shown at the top of the MultiQC report) ----------
# per-sample primer removal check (from PRIMER_TRIM)
if (!is.null(pstat)) {
  pcols <- intersect(c("sample", "input", "with_primers", "with_primers_pct", "swapped_pct", "spacer_max",
                       "primer_left_R1_pct", "primer_left_R2_pct"), colnames(pstat))
  write_mqc_table(pstat[order(pstat$sample), pcols], "primer_removal_mqc.tsv", "primer_removal", "Primer removal per sample",
                  "Read pairs with both primers (kept), share reverse-oriented (swapped), and the check after trimming: primer_left = % of trimmed reads that still start with the primer (should be ~0).")
  left <- suppressWarnings(max(c(pstat$primer_left_R1_pct, pstat$primer_left_R2_pct), na.rm = TRUE))
  if (is.finite(left) && left > 1)
    primer_msgs <- c(primer_msgs, sprintf("Up to %.1f%% of trimmed reads still start with a primer: primer removal was incomplete (see 'Primer removal per sample').", left))
}
qc_warn <- c(primer_msgs,
  if (!is.na(med_merge) && med_merge < 70) "Median merge rate < 70%: truncLen may be too short for the overlap.",
  if (chim_kept < 0.8) "< 80% of reads survived chimera removal: primers may still be in the reads (see primer detection).",
  if (any(!keep)) paste("Samples with no reads after filtering (dropped):", paste(man$sample[!keep], collapse = ", ")),
  if (length(dropped_primer)) paste("Samples with no reads after primer removal (dropped):", paste(dropped_primer, collapse = ", ")))
trunc_report <- rbind(trunc_report, data.frame(
  parameter = c("median merge rate", "reads kept after chimera removal", "QC warnings"),
  value = c(sprintf("%.1f%%", med_merge), sprintf("%.1f%%", 100 * chim_kept),
            if (length(qc_warn)) paste(qc_warn, collapse = " | ") else "none")))
write_tsv(trunc_report, "dada2_truncation.tsv")
write_mqc_table(trunc_report, "dada2_truncation_mqc.tsv", "dada2_truncation", "DADA2 summary and warnings",
                "Primer lengths, truncation, overlap check, merge and chimera rates. Read the QC warnings row first.")

steps <- c("input", "primers_found", "filtered", "denoisedF", "merged", "nonchim")
tl <- do.call(rbind, lapply(steps, function(s) data.frame(sample = track$sample, step = s,
                                                         pct = 100 * track[[s]] / track$input)))
tl$step <- factor(tl$step, levels = steps,
                  labels = c("input", "primers found", "filtered", "denoised", "merged", "non-chimeric"))
p <- ggplot(tl, aes(step, pct, group = sample)) +
  geom_line(colour = "grey60", alpha = 0.7) + geom_point(size = 1.5, colour = "#0072B2") +
  scale_y_continuous(limits = c(0, 100)) +
  labs(title = "Reads retained through DADA2",
       subtitle = sprintf("median merge rate %.1f%% | %.1f%% of reads kept after chimera removal", med_merge, 100 * chim_kept),
       x = NULL, y = "% of input read pairs") + theme_amp()
save_fig(p, "dada2_read_tracking_plot", 8, 5, mqc = TRUE)

lens <- data.frame(length = nchar(colnames(seqtab.nochim)), reads = colSums(seqtab.nochim))
main_len <- lens$length[which.max(lens$reads)]
near <- sum(lens$reads[abs(lens$length - main_len) <= 5]) / max(1, sum(lens$reads))
p <- ggplot(aggregate(reads ~ length, lens, sum), aes(length, reads)) + geom_col(fill = "#0072B2", width = 1) +
  scale_y_continuous(labels = scales::comma) +
  scale_x_continuous(limits = c(min(min(lens$length), main_len - 20) - 1, max(max(lens$length), main_len + 20) + 1),
                     breaks = function(l) unique(round(scales::extended_breaks()(l)))) +
  labs(title = "Merged sequences: ASV length (R1 + R2 merged, without primers)",
       subtitle = sprintf("%s ASVs | main peak %d bp | %.1f%% of reads within 5 bp of it\nOne main peak at the region length is expected (V4 ~253 bp, V3-V4 ~400-430 bp); distant peaks are often off-target",
                          format(nrow(lens), big.mark = ","), main_len, 100 * near),
       x = "ASV length (bp)", y = "Reads") + theme_amp()
save_fig(p, "dada2_asv_length", 9, 4.8, mqc = TRUE)

write_session("dada2")
