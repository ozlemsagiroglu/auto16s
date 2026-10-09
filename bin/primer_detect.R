#!/usr/bin/env Rscript
# Detect which 16S primer pair is at the start of the reads (or use the pair given by the user).
# Looks at the first reads of every sample and matches them against assets/primers_16s.tsv.
.script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(.script)), "utils.R"))
source(file.path(dirname(normalizePath(.script)), "primer_utils.R"))
opt <- parse_args()
suppressPackageStartupMessages(library(ggplot2))

man <- read.delim(opt$manifest, stringsAsFactors = FALSE)
lib <- read.delim(opt$primer_db, comment.char = "#", stringsAsFactors = FALSE)
user_fp <- if (is.na(opt$fw_primer)) NA else toupper(opt$fw_primer)
user_rp <- if (is.na(opt$rv_primer)) NA else toupper(opt$rv_primer)
if (xor(is.na(user_fp), is.na(user_rp))) stop("Give both --fw_primer and --rv_primer, or neither.")
if (!is.na(user_fp)) lib <- rbind(data.frame(name = c("user forward", "user reverse"), direction = c("F", "R"),
                                             region = "user", sequence = c(user_fp, user_rp)), lib)

# ---------- first reads of every sample (pairs stay in sync) ----------
n_each <- max(300, min(2000, ceiling(40000 / nrow(man))))
first_reads <- function(f) { s <- FastqStreamer(f, n = n_each); on.exit(close(s)); sread(yield(s)) }
rF <- do.call(c, unname(lapply(man$r1, first_reads)))
rR <- do.call(c, unname(lapply(man$r2, first_reads)))
n  <- min(length(rF), length(rR)); rF <- rF[seq_len(n)]; rR <- rR[seq_len(n)]

H1 <- vapply(lib$sequence, function(p) !is.na(primer_end(rF, p)), logical(n))   # reads x primers, R1 starts
H2 <- vapply(lib$sequence, function(p) !is.na(primer_end(rR, p)), logical(n))   # R2 starts
lib$R1_pct <- round(100 * colMeans(H1), 1); lib$R2_pct <- round(100 * colMeans(H2), 1)

# ---------- best forward/reverse pair, either orientation ----------
fwd_i <- which(lib$direction == "F"); rev_i <- which(lib$direction == "R")
if (!is.na(user_fp)) { fwd_i <- 1L; rev_i <- 2L }
# exact-match rate per primer: tie-break between near-identical primers (515F / 515F-Y, 806R / 806RB)
exact_rate <- setNames(vapply(seq_len(nrow(lib)), function(i) mean(!is.na(primer_end(rF, lib$sequence[i], 0)) |
                                                                  !is.na(primer_end(rR, lib$sequence[i], 0))), 1),
                       seq_len(nrow(lib)))
best <- list(score = 0, exact = 0)
for (f in fwd_i) for (r in rev_i) {
  if (max(lib$R1_pct[c(f, r)], lib$R2_pct[c(f, r)]) < 2) next
  normal  <- H1[, f] & H2[, r]
  swapped <- !normal & H1[, r] & H2[, f]
  score   <- mean(normal | swapped)
  exact <- exact_rate[[as.character(f)]] + exact_rate[[as.character(r)]]
  if (score + 1e-3 * exact > best$score + 1e-3 * best$exact)
    best <- list(f = f, r = r, score = score, exact = exact, swapped = mean(swapped) / max(score, 1e-9))
}

top_prefix <- function(x) { t <- sort(table(as.character(subseq(x[width(x) >= 20], 1, 20))), decreasing = TRUE)
  sprintf("%s (%.0f%% of reads)", names(t)[1], 100 * t[1] / sum(t)) }

warnings <- character(0)
if (best$score >= 0.2) {
  status <- if (!is.na(user_fp)) "user" else "detected"
  fp <- lib$sequence[best$f]; rp <- lib$sequence[best$r]
  fname <- lib$name[best$f]; rname <- lib$name[best$r]
  region <- if (status == "user") "user" else paste(unique(c(lib$region[best$f], lib$region[best$r])), collapse = "-")
  if (best$score < 0.7) warnings <- c(warnings, sprintf(
    "Only %.0f%% of read pairs carry the %s/%s primers. Pairs without them are removed; check whether these are the right primers.",
    100 * best$score, fname, rname))
  message(sprintf("Primers: %s (%s) / %s (%s), %s; in %.1f%% of read pairs, %.1f%% of them in reverse orientation",
                  fname, fp, rname, rp, region, 100 * best$score, 100 * best$swapped))
} else {
  status <- "none"; fp <- rp <- fname <- rname <- region <- NA
  if (!is.na(user_fp)) stop("The given primers were found in only ", round(100 * best$score, 1), "% of read pairs.")
  warnings <- c(warnings, paste0(
    "No known 16S primer pair found at the read starts: assuming primers were already removed. ",
    "If your primers are not in assets/primers_16s.tsv, give them with --fw_primer/--rv_primer. ",
    "Most common read start R1: ", top_prefix(rF), "; R2: ", top_prefix(rR), "."))
}
for (w in warnings) message("WARNING: ", w)

res <- data.frame(key = c("status", "fw_name", "fw_primer", "rv_name", "rv_primer", "region",
                          "pairs_with_primers_pct", "reverse_orientation_pct", "warnings"),
                  value = c(status, fname, fp, rname, rp, region,
                            sprintf("%.1f", 100 * best$score),
                            if (status == "none") NA else sprintf("%.1f", 100 * best$swapped),
                            if (length(warnings)) paste(warnings, collapse = " | ") else "none"))
write_tsv(res, "primers.tsv")
shown <- res; shown$key <- c("status", "forward primer", "forward sequence", "reverse primer", "reverse sequence",
                             "16S region", "read pairs with both primers (%)", "of these in reverse orientation (%)", "warnings")
write_mqc_table(shown, "primer_detection_mqc.tsv", "primer_detection", "Primer detection",
                "Primers identified at the start of the reads and removed by sequence (spacers and reverse-oriented pairs are handled).")

# ---------- plot ----------
lib$label <- paste0(lib$name, " (", lib$region, ")")
lib$chosen <- seq_len(nrow(lib)) %in% c(best$f, best$r)
show <- lib[order(-(pmax(lib$R1_pct, lib$R2_pct))), ][seq_len(min(8, nrow(lib))), ]
long <- rbind(data.frame(label = show$label, read = "R1 start", pct = show$R1_pct, chosen = show$chosen),
              data.frame(label = show$label, read = "R2 start", pct = show$R2_pct, chosen = show$chosen))
long$label <- factor(long$label, levels = rev(show$label))
p <- ggplot(long, aes(pct, label, fill = chosen)) + geom_col(width = 0.7) +
  facet_wrap(~ read) + scale_x_continuous(limits = c(0, 100)) +
  scale_fill_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "grey70"), labels = c(`TRUE` = "selected", `FALSE` = "other"), name = NULL) +
  labs(title = "Primer detection",
       subtitle = if (status == "none") "No known primer pair found: reads used as they are" else
         sprintf("%s / %s (%s): in %.0f%% of read pairs, %.0f%% of them reverse-oriented",
                 fname, rname, region, 100 * best$score, 100 * best$swapped),
       x = "% of reads starting with the primer (up to 12 bases offset, 1-2 mismatches)", y = NULL) +
  theme_amp() + theme(legend.position = "bottom")
save_fig(p, "primer_detection", 9, 5, mqc = TRUE)
