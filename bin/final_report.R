#!/usr/bin/env Rscript
# Final report: one self-contained HTML file that reads like a short study report
# (samples -> methods with the values used in this run -> results -> files -> software and references).
# Only base R is needed for the HTML itself; figures are embedded as base64 PNG.
.script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(.script)), "utils.R"))
opt <- parse_args()

dir  <- opt$inputs
gcol <- opt$group_col
qcut <- as.numeric(opt$qval)
silva_file <- basename(opt$silva_db)
silva_ver  <- regmatches(silva_file, regexpr("(?<=v)[0-9]+(\\.[0-9]+)?", silva_file, perl = TRUE))
silva_txt  <- if (length(silva_ver) && grepl("silva", silva_file, ignore.case = TRUE)) sprintf("SILVA %s (%s)", silva_ver, silva_file) else silva_file
f    <- function(name) file.path(dir, name)
has  <- function(name) file.exists(f(name)) && file.size(f(name)) > 0
rd   <- function(name, ...) if (has(name)) read.delim(f(name), check.names = FALSE, stringsAsFactors = FALSE, ...) else NULL
rdc  <- function(name) if (has(name)) read.csv(f(name), check.names = FALSE, stringsAsFactors = FALSE) else NULL

# ---------------------------------------------------------------- helpers
esc <- function(x) { x <- gsub("&", "&amp;", as.character(x), fixed = TRUE); x <- gsub("<", "&lt;", x, fixed = TRUE); gsub(">", "&gt;", x, fixed = TRUE) }
num <- function(x, d = 0) ifelse(is.na(x), "NA", formatC(as.numeric(x), format = "f", digits = d, big.mark = ","))
pct <- function(x, d = 1) ifelse(is.na(x), "NA", sprintf(paste0("%.", d, "f%%"), as.numeric(x)))
pfmt <- function(p) ifelse(is.na(p), "NA", ifelse(p < 0.001, "&lt; 0.001", sprintf("= %.3f", p)))

b64 <- function(path) {
  r <- readBin(path, "raw", file.size(path)); n <- length(r); pad <- (3 - n %% 3) %% 3
  x <- matrix(as.integer(c(r, as.raw(rep(0L, pad)))), nrow = 3)
  v <- x[1, ] * 65536 + x[2, ] * 256 + x[3, ]
  ch <- c(LETTERS, letters, 0:9, "+", "/")
  out <- ch[rbind(v %/% 262144, (v %/% 4096) %% 64, (v %/% 64) %% 64, v %% 64) + 1]
  if (pad) out[(length(out) - pad + 1):length(out)] <- "="
  paste(out, collapse = "")
}

fig_no <- 0; tab_no <- 0
figure <- function(name, caption) {
  if (!has(name)) return("")
  fig_no <<- fig_no + 1
  # display at the size the figure was designed for (PNG width / 300 dpi, ~100 px per inch), never wider than the page
  hdr <- readBin(f(name), "raw", 24)
  px <- sum(as.integer(hdr[17:20]) * 256^(3:0))
  sprintf('<figure><img src="data:image/png;base64,%s" alt="%s" style="width:%dpx"><figcaption><b>Figure @FIG@.</b> %s</figcaption></figure>\n',
          b64(f(name)), esc(name), as.integer(round(px / 300 * 100)), caption)
}
table_html <- function(df, caption = NULL, max_rows = 50, digits = 3) {
  if (is.null(df) || !nrow(df)) return("")
  note <- ""
  if (nrow(df) > max_rows) { note <- sprintf('<p class="note">First %d of %d rows; the full table is in the output folder.</p>', max_rows, nrow(df)); df <- head(df, max_rows) }
  cells <- lapply(df, function(col) {
    if (is.numeric(col)) {
      ifelse(is.na(col), "", ifelse(abs(col - round(col)) < 1e-9 & abs(col) < 1e9, num(col, 0), formatC(signif(col, digits), format = "g", digits = digits)))
    } else esc(ifelse(is.na(col), "", col))
  })
  right <- vapply(df, is.numeric, TRUE)
  head_row <- paste0("<tr>", paste0("<th", ifelse(right, ' class="r"', ""), ">", esc(names(df)), "</th>", collapse = ""), "</tr>")
  rows <- vapply(seq_len(nrow(df)), function(i)
    paste0("<tr>", paste0("<td", ifelse(right, ' class="r"', ""), ">", vapply(cells, `[`, "", i), "</td>", collapse = ""), "</tr>"), "")
  cap <- ""
  if (!is.null(caption)) { tab_no <<- tab_no + 1; cap <- sprintf("<caption><b>Table @TAB@.</b> %s</caption>", caption) }
  sprintf('<div class="tw"><table>%s<thead>%s</thead><tbody>%s</tbody></table></div>%s\n', cap, head_row, paste(rows, collapse = ""), note)
}
kv <- function(df) if (is.null(df)) list() else setNames(as.list(df[[2]]), df[[1]])

# package versions from the sessionInfo files written by every R step
versions <- list()
for (s in list.files(dir, pattern = "_sessionInfo\\.txt$", full.names = TRUE)) {
  txt <- readLines(s, warn = FALSE)
  r <- grep("^R version", txt, value = TRUE)
  if (length(r)) versions[["R"]] <- sub("^R version ([0-9.]+).*", "\\1", r[1])
  for (m in regmatches(txt, gregexpr("[A-Za-z][A-Za-z0-9.]*_[0-9][0-9.\\-]*", txt))) for (pv in m) {
    k <- sub("_.*", "", pv); versions[[k]] <- sub("^[^_]*_", "", pv)
  }
}
ver <- function(k) if (!is.null(versions[[k]])) versions[[k]] else "?"
fastqc_ver <- "?"
zips <- list.files(dir, pattern = "_fastqc\\.zip$", full.names = TRUE)
if (length(zips)) {
  inner <- grep("fastqc_data.txt$", utils::unzip(zips[1], list = TRUE)$Name, value = TRUE)[1]
  con <- unz(zips[1], inner); fl <- readLines(con, n = 1, warn = FALSE); close(con)
  fastqc_ver <- sub("^##FastQC\\s+", "", fl)
}
multiqc_ver <- "?"
if (has("multiqc_data.json")) {
  j <- readLines(f("multiqc_data.json"), warn = FALSE)
  m <- regmatches(j, regexpr('"config_version": *"[^"]+"', j))
  if (length(m)) multiqc_ver <- sub('.*"([^"]+)"$', "\\1", m[1])
}

# ---------------------------------------------------------------- data
sheet    <- read.csv(opt$samplesheet, check.names = FALSE, stringsAsFactors = FALSE)
primers  <- kv(rd("primers.tsv"))
trunc    <- kv(rd("dada2_truncation.tsv"))
track    <- rd("dada2_read_tracking.tsv")
taxf     <- rd("taxonomic_filtering.tsv")
depths   <- rd("sample_depths.tsv")
alpha_t  <- rd("alpha_diversity_tests.tsv")
beta_t   <- rd("beta_diversity_tests.tsv")
pairw    <- rd("permanova_pairwise.tsv")
interp   <- if (has("beta_diversity_interpretation.txt")) paste(readLines(f("beta_diversity_interpretation.txt"), warn = FALSE), collapse = " ") else ""
m3_all   <- rdc("maaslin3_all_results.csv")
m3_sig   <- rdc("maaslin3_significant.csv")
genus_ra <- rd("genus_relative_abundance.tsv")
asv_tab  <- rd("asv_table.tsv")
rare_tab <- rd("asv_table_rarefied.tsv")
pstats   <- do.call(rbind, lapply(list.files(dir, pattern = "\\.primer_stats\\.tsv$", full.names = TRUE), read.delim))

n_samples <- nrow(sheet)
groups    <- table(sheet[[gcol]])
grp_txt   <- paste(sprintf("%s (%d)", names(groups), as.integer(groups)), collapse = ", ")
rare_depth <- if (!is.null(rare_tab) && ncol(rare_tab) > 1) sum(rare_tab[[2]]) else NA
n_asv     <- if (!is.null(asv_tab)) nrow(asv_tab) else NA
n_genera  <- if (!is.null(genus_ra)) sum(!grepl("unclassified", genus_ra[[1]])) else NA
excluded  <- if (!is.null(depths)) depths$sample[!as.logical(depths$in_alpha_beta)] else character(0)
ref_level <- if (!is.null(m3_all) && nrow(m3_all)) sub("^.* vs ", "", m3_all$comparison[1]) else NA
input_total <- if (!is.null(track)) sum(track$input) else NA

# ---------------------------------------------------------------- summary
warn <- character(0)
if (!is.null(trunc[["QC warnings"]]) && trunc[["QC warnings"]] != "none") warn <- c(warn, trunc[["QC warnings"]])
if (!is.null(primers$warnings) && primers$warnings != "none") warn <- c(warn, primers$warnings)
if (length(excluded)) warn <- c(warn, sprintf("Left out of alpha and beta diversity (too few reads): %s", paste(excluded, collapse = ", ")))
if (any(groups < 3)) warn <- c(warn, "At least one group has fewer than 3 samples; statistical tests have very little power.")

find_lines <- character(0)
if (!is.null(alpha_t) && nrow(alpha_t)) {
  sig_a <- alpha_t[!is.na(alpha_t$p) & alpha_t$p < 0.05, ]
  find_lines <- c(find_lines, if (nrow(sig_a))
    sprintf("Alpha diversity differed between groups for %s (p &lt; 0.05).", paste(unique(sig_a$measure), collapse = ", "))
  else "Alpha diversity (Observed, Shannon, Simpson) did not differ significantly between groups (all p &ge; 0.05).")
}
if (!is.null(beta_t) && nrow(beta_t)) {
  pm <- beta_t[grepl("PERMANOVA", beta_t$test), ][1, ]
  r2 <- suppressWarnings(as.numeric(sub("^R2 = ([0-9.]+).*", "\\1", pm$statistic)))
  find_lines <- c(find_lines, sprintf("Community composition (Bray-Curtis): PERMANOVA R<sup>2</sup> = %.2f, p %s, i.e. %.0f%% of the variation between samples is explained by %s. %s",
                                      r2, pfmt(pm$p), 100 * r2, esc(gcol), esc(interp)))
}
if (!is.null(m3_all)) {
  na_ <- if (is.null(m3_sig)) 0 else sum(m3_sig$model == "abundance"); np_ <- if (is.null(m3_sig)) 0 else sum(m3_sig$model == "prevalence")
  find_lines <- c(find_lines, sprintf("MaAsLin 3 tested %s genera: %d differed in abundance and %d in prevalence (q &lt; %g, reference group: %s).",
                                      num(length(unique(m3_all$genus))), na_, np_, qcut, esc(ref_level)))
}

summary_html <- paste0(
  '<div class="cards">',
  sprintf('<div class="card"><div class="v">%d</div><div class="k">samples</div><div class="s">%s</div></div>', n_samples, esc(grp_txt)),
  sprintf('<div class="card"><div class="v">%s</div><div class="k">read pairs in</div><div class="s">median %s per sample</div></div>', num(input_total), num(median(track$input))),
  sprintf('<div class="card"><div class="v">%s</div><div class="k">kept to the end</div><div class="s">median of samples, after chimera removal</div></div>', pct(median(track$retained_pct))),
  sprintf('<div class="card"><div class="v">%s</div><div class="k">ASVs</div><div class="s">%s named genera</div></div>', num(n_asv), num(n_genera)),
  sprintf('<div class="card"><div class="v">%s</div><div class="k">rarefaction depth</div><div class="s">reads per sample for diversity</div></div>', num(rare_depth)),
  '</div>',
  "<h3>Main results</h3><ul>", paste0("<li>", find_lines, "</li>", collapse = ""), "</ul>",
  if (length(warn)) paste0('<div class="warn"><b>Check:</b><ul>', paste0("<li>", esc(warn), "</li>", collapse = ""), "</ul></div>")
  else '<p class="ok">No QC warnings.</p>',
  '<p class="note">The statements above report the test results only; their biological interpretation is left to the reader.</p>')

# ---------------------------------------------------------------- samples
meta_cols <- setdiff(names(sheet), c("fastq_1", "fastq_2"))
samp <- sheet[, meta_cols, drop = FALSE]
if (!is.null(track)) samp <- merge(samp, track[, c("sample", "input", "nonchim", "retained_pct")], by = "sample", all.x = TRUE, sort = FALSE)
if (!is.null(depths)) samp <- merge(samp, data.frame(sample = depths$sample, reads_after_taxonomic_filter = depths$reads,
                                                     in_diversity = ifelse(as.logical(depths$in_alpha_beta), "yes", "no")), by = "sample", all.x = TRUE, sort = FALSE)
names(samp)[names(samp) == "input"] <- "read_pairs_in"; names(samp)[names(samp) == "nonchim"] <- "read_pairs_final"
names(samp)[names(samp) == "retained_pct"] <- "kept_pct"
samp <- samp[order(samp[[gcol]], samp$sample), ]

# ---------------------------------------------------------------- methods
pr_txt <- if (identical(primers$status, "none")) {
  "No known 16S primer pair was found at the start of the reads; the reads were used as delivered (primers already removed)."
} else sprintf("The primers were identified automatically by matching the read starts against a library of common 16S primers: %s (%s) and %s (%s), 16S region %s%s. They were removed from every read by sequence (IUPAC codes, 1-2 mismatches, variable-length spacers and reverse-oriented pairs handled); read pairs without both primers were discarded (%s of read pairs carried both).",
             esc(primers$fw_name), esc(primers$fw_primer), esc(primers$rv_name), esc(primers$rv_primer), esc(primers$region),
             if (identical(primers$status, "user")) " (given by the user)" else "", pct(as.numeric(primers$pairs_with_primers_pct)))
methods_html <- paste0(
  "<p>", sprintf("Raw paired-end reads (%d samples) were checked with FastQC %s. ", n_samples, esc(fastqc_ver)), pr_txt, "</p>",
  "<p>", sprintf("Amplicon sequence variants (ASVs) were inferred with DADA2 %s in R %s. Reads were truncated at %s (R1) and %s (R2) bases; these positions were chosen from the data as the pair that keeps the most read pairs through the quality filter while leaving R1 and R2 overlapping by at least 20 bp over the amplicon (estimated length %s bp without primers). Reads were filtered with <code>filterAndTrim</code> (maxEE = %s, truncQ = 2, maxN = 0, PhiX removal). Error models were learned separately for R1 and R2, reads were denoised with <code>dada</code>, merged with <code>mergePairs</code> (minimum overlap %s bp) and chimeras were removed with <code>removeBimeraDenovo</code> (consensus). Taxonomy was assigned with <code>assignTaxonomy</code> (minBoot = 50) against the %s reference.",
                 ver("dada2"), ver("R"), esc(trunc[["truncLen R1"]]), esc(trunc[["truncLen R2"]]), esc(trunc[["estimated amplicon length (without primers)"]]),
                 esc(opt$max_ee), esc(opt$min_overlap), esc(silva_txt)), "</p>",
  "<p>", sprintf("The ASV table was combined with the sample metadata in phyloseq %s. ASVs without a kingdom or phylum assignment, eukaryotic ASVs, chloroplasts and mitochondria were removed. For alpha and beta diversity, samples were rarefied to %s reads (random seed %s)%s. Alpha diversity (observed ASVs, Shannon and Simpson indices) was compared between groups with the %s. Beta diversity was computed as Bray-Curtis dissimilarity, ordinated by principal coordinates analysis (PCoA) and tested with PERMANOVA (<code>adonis2</code>, %s permutations) and with <code>betadisper</code> for differences in within-group dispersion (vegan %s).",
                 ver("phyloseq"), num(rare_depth), esc(opt$seed),
                 if (length(excluded)) sprintf("; %d sample(s) below this depth were left out", length(excluded)) else "",
                 if (length(groups) > 2) "Kruskal-Wallis test followed by pairwise Wilcoxon tests (Benjamini-Hochberg)" else "Wilcoxon rank-sum test",
                 esc(opt$permutations), ver("vegan")), "</p>",
  "<p>", sprintf("Differential abundance and prevalence at genus level were tested on the unrarefied counts with MaAsLin 3 %s (total-sum scaling, log transformation, genera present in at least %s%% of samples, reference group: %s). MaAsLin 3 fits an abundance model (samples in which the genus was found) and a prevalence model (presence/absence); associations with a Benjamini-Hochberg q-value &lt; %g were considered significant. Figures were made with ggplot2 %s.",
                 ver("maaslin3"), num(100 * as.numeric(opt$min_prevalence)), esc(ref_level), qcut, ver("ggplot2")), "</p>",
  "<p>", sprintf("The analysis was run with the auto16s pipeline %s (Nextflow %s). The QC summary report was made with MultiQC %s.",
                 esc(opt$pipeline_version), esc(opt$nextflow_version), esc(multiqc_ver)), "</p>")

param_df <- data.frame(
  setting = c("Group column", "Primers", "truncLen R1 / R2", "Expected R1/R2 overlap (bp)", "Read pairs expected to pass filtering",
              "maxEE", "Minimum merge overlap (bp)", "Taxonomy reference", "Rarefaction depth / seed", "PERMANOVA permutations",
              "MaAsLin 3 minimum prevalence / q-value / reference"),
  value = c(gcol, if (identical(primers$status, "none")) "none found" else sprintf("%s / %s", primers$fw_name, primers$rv_name),
            sprintf("%s / %s", trunc[["truncLen R1"]], trunc[["truncLen R2"]]), trunc[["expected overlap (bp)"]],
            if (is.null(trunc[["read pairs expected to pass filtering"]])) "NA" else trunc[["read pairs expected to pass filtering"]],
            opt$max_ee, opt$min_overlap, silva_txt, sprintf("%s / %s", num(rare_depth), opt$seed), opt$permutations,
            sprintf("%s / %g / %s", opt$min_prevalence, qcut, ref_level)), stringsAsFactors = FALSE)

# ---------------------------------------------------------------- results
top_genera <- NULL
if (!is.null(genus_ra)) {
  ra <- as.matrix(genus_ra[, -1, drop = FALSE]); rownames(ra) <- genus_ra[[1]]
  g  <- setNames(sheet[[gcol]], sheet$sample)[colnames(ra)]
  means <- sapply(split(seq_len(ncol(ra)), g), function(ix) rowMeans(ra[, ix, drop = FALSE]))
  if (is.null(dim(means))) means <- matrix(means, ncol = 1, dimnames = list(rownames(ra), unique(g)))
  o <- order(-rowMeans(means))
  top_genera <- data.frame(genus = rownames(means)[o], round(100 * means[o, , drop = FALSE], 2), check.names = FALSE)
  names(top_genera)[-1] <- paste0(names(top_genera)[-1], " (%)")
  top_genera <- head(top_genera, 15)
}
m3_tab <- NULL
if (!is.null(m3_sig) && nrow(m3_sig)) {
  m3_tab <- m3_sig[order(m3_sig$model, m3_sig$qval_individual), intersect(c("genus", "model", "comparison", "coef", "stderr", "qval_individual", "N_not_zero"), names(m3_sig))]
  names(m3_tab)[names(m3_tab) == "qval_individual"] <- "q"
  names(m3_tab)[names(m3_tab) == "N_not_zero"] <- "samples_with_genus"
}
alpha_show <- alpha_t
if (!is.null(alpha_show)) alpha_show <- alpha_show[, colSums(!is.na(alpha_show)) > 0, drop = FALSE]

results_html <- paste0(
  '<h3 id="r-qc">4.1 Read processing</h3>',
  sprintf("<p>Of %s read pairs, a median of %s per sample remained after primer removal, quality filtering, denoising, merging and chimera removal. The median merge rate was %s and %s of the merged reads were kept after chimera removal.</p>",
          num(input_total), pct(median(track$retained_pct)), esc(trunc[["median merge rate"]]), esc(trunc[["reads kept after chimera removal"]])),
  table_html(track, "Read pairs per sample after each step (primers found, filtered, denoised R1/R2, merged, non-chimeric)."),
  if (!is.null(pstats) && "primer_left_R1_pct" %in% names(pstats))
    sprintf("<p>Primer removal check: after removal, at most %s of R1 and %s of R2 reads still contained the primer.</p>",
            pct(max(pstats$primer_left_R1_pct, na.rm = TRUE), 2), pct(max(pstats$primer_left_R2_pct, na.rm = TRUE), 2)) else "",
  figure("dada2_quality_truncation.png", "Read quality after primer removal (median and interquartile range per position) and the chosen truncation positions."),
  figure("dada2_read_tracking_plot.png", "Share of read pairs retained through the DADA2 steps; one line per sample."),
  figure("dada2_asv_length.png", "Length of the merged sequences (ASVs). One main peak at the length of the amplified region is expected."),
  '<h3 id="r-tax">4.2 Taxonomic filtering and sequencing depth</h3>',
  table_html(taxf, "ASVs and reads removed before the analysis."),
  figure("sample_depths_plot.png", "Reads per sample after taxonomic filtering; the dashed line marks the rarefaction depth."),
  figure("rarefaction_curves.png", "Rarefaction curves: observed ASVs against sequencing depth. Curves that level off indicate sufficient depth."),
  '<h3 id="r-comp">4.3 Taxonomic composition</h3>',
  table_html(top_genera, "The 15 most abundant genera: mean relative abundance per group (all reads).", max_rows = 15),
  figure("phylum_barplot_groups.png", "Mean phylum composition per group."),
  figure("genus_barplot_groups.png", "Mean genus composition per group."),
  figure("genus_barplot_samples.png", "Genus composition per sample."),
  figure("genus_heatmap.png", "Relative abundance of the most abundant genera per sample (log10 scale)."),
  '<h3 id="r-alpha">4.4 Alpha diversity</h3>',
  figure("alpha_diversity.png", sprintf("Alpha diversity per group (rarefied to %s reads).", num(rare_depth))),
  table_html(alpha_show, "Alpha diversity tests."),
  '<h3 id="r-beta">4.5 Beta diversity</h3>',
  sprintf("<p>%s</p>", esc(interp)),
  figure("pcoa_bray.png", "Principal coordinates analysis of Bray-Curtis dissimilarities. Diamonds mark the group centroids; lines join each sample to its centroid; dashed outlines show the area covered by each group."),
  figure("betadisper_bray.png", "Distance of each sample to its group centroid (within-group dispersion, betadisper)."),
  table_html(beta_t, "PERMANOVA and betadisper results."),
  table_html(pairw, "Pairwise PERMANOVA between groups."),
  '<h3 id="r-da">4.6 Differential abundance and prevalence (MaAsLin 3)</h3>',
  if (is.null(m3_tab)) sprintf("<p>No genus differed significantly between groups (q &lt; %g).</p>", qcut)
  else table_html(m3_tab, sprintf("Genera with q &lt; %g. abundance: differs in abundance where present; prevalence: differs in the share of samples in which it was found. Positive coefficient: higher than in %s.", qcut, esc(ref_level)), max_rows = 40),
  figure("maaslin3_volcano.png", "MaAsLin 3 coefficients against significance, separately for the abundance and the prevalence model."),
  figure("maaslin3_coefficients.png", "Coefficients (with standard error) of the significant genera."),
  figure("maaslin3_abundance_boxplots.png", "Relative abundance of the significant abundance hits, in the samples where the genus was found."),
  figure("maaslin3_prevalence_bars.png", "Share of samples per group in which the significant prevalence hits were found."))

files_df <- data.frame(
  folder = c("report/", "multiqc/", "phyloseq/", "composition/", "alpha_diversity/", "beta_diversity/", "maaslin3/", "dada2/", "primers/", "taxonomy/", "fastqc/", "pipeline_info/"),
  content = c("This report", "QC summary report (MultiQC)", "phyloseq objects (all reads and rarefied), ASV tables, taxonomy, ASV sequences, metadata",
              "Counts and relative abundances per phylum, family and genus; figures", "Values per sample, tests, figure",
              "Bray-Curtis distances, PCoA coordinates, tests, figures", "All and significant results, figures, raw MaAsLin 3 output",
              "Read tracking, truncation, error models, ASV lengths", "Detected primers and per-sample statistics",
              "SILVA assignment per ASV", "FastQC reports before and after filtering", "Nextflow execution report, timeline, trace"),
  stringsAsFactors = FALSE)
sw_df <- data.frame(software = c("auto16s", "Nextflow", "FastQC", "R", "DADA2", "phyloseq", "vegan", "MaAsLin 3", "ggplot2", "MultiQC"),
                    version = c(opt$pipeline_version, opt$nextflow_version, fastqc_ver, ver("R"), ver("dada2"), ver("phyloseq"),
                                ver("vegan"), ver("maaslin3"), ver("ggplot2"), multiqc_ver), stringsAsFactors = FALSE)
refs <- c(
  "Callahan B.J. et al. (2016). DADA2: High-resolution sample inference from Illumina amplicon data. <i>Nat Methods</i> 13:581-583.",
  "Quast C. et al. (2013). The SILVA ribosomal RNA gene database project. <i>Nucleic Acids Res</i> 41:D590-D596.",
  "McMurdie P.J. &amp; Holmes S. (2013). phyloseq: an R package for reproducible interactive analysis and graphics of microbiome census data. <i>PLoS ONE</i> 8:e61217.",
  "Oksanen J. et al. vegan: Community Ecology Package. https://github.com/vegandevs/vegan",
  "Anderson M.J. (2001). A new method for non-parametric multivariate analysis of variance. <i>Austral Ecol</i> 26:32-46.",
  "Nickols W.A. et al. (2026). MaAsLin 3: refining and extending generalized multivariate linear models for meta-omic association discovery. <i>Nat Methods</i>.",
  "Wickham H. (2016). ggplot2: Elegant Graphics for Data Analysis. Springer.",
  "Andrews S. (2010). FastQC: a quality control tool for high throughput sequence data.",
  "Ewels P. et al. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. <i>Bioinformatics</i> 32:3047-3048.",
  "Di Tommaso P. et al. (2017). Nextflow enables reproducible computational workflows. <i>Nat Biotechnol</i> 35:316-319.",
  "Weinstein M.M. et al. (2019). FIGARO: an efficient and objective tool for optimizing microbiome rRNA gene trimming parameters. <i>bioRxiv</i> 610394.")

# ---------------------------------------------------------------- page
css <- "
:root{--fg:#1d1d1f;--mut:#5f6368;--line:#e3e3e3;--bg:#fff;--soft:#f6f7f9;--acc:#0072B2;--warn:#fff4e5;--warnb:#e69f00}
@media (prefers-color-scheme: dark){:root{--fg:#e8e8e8;--mut:#a0a4a8;--line:#3a3a3a;--bg:#161616;--soft:#202124;--acc:#56B4E9;--warn:#3a2e14;--warnb:#e69f00}
 figure img{background:#fff;border-radius:4px}}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--fg);font:15px/1.6 -apple-system,'Segoe UI',Roboto,Helvetica,Arial,sans-serif}
main{max-width:1020px;margin:0 auto;padding:24px 16px 64px}
h1{font-size:28px;margin:0 0 4px}h2{font-size:21px;margin:44px 0 10px;padding-top:8px;border-top:1px solid var(--line)}h3{font-size:17px;margin:28px 0 8px}
.sub{color:var(--mut);margin:0 0 20px}nav{background:var(--soft);border-radius:8px;padding:10px 16px;margin:16px 0}nav a{color:var(--acc);text-decoration:none;margin-right:14px;white-space:nowrap}
.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(160px,1fr));gap:10px;margin:12px 0}
.card{background:var(--soft);border-radius:8px;padding:12px 14px}.card .v{font-size:24px;font-weight:700}.card .k{font-weight:600}.card .s{color:var(--mut);font-size:13px}
.warn{background:var(--warn);border-left:4px solid var(--warnb);padding:8px 14px;border-radius:4px;margin:12px 0}.warn ul{margin:4px 0}
.ok{color:#009E73;font-weight:600}.note{color:var(--mut);font-size:13px}
figure{margin:18px 0}figure img{max-width:100%;height:auto;display:block;margin:0 auto}figcaption{color:var(--mut);font-size:13px;margin-top:6px}
.tw{overflow-x:auto;margin:12px 0}table{border-collapse:collapse;font-size:13px;width:100%}caption{text-align:left;color:var(--mut);padding-bottom:6px;caption-side:top}
th,td{border-bottom:1px solid var(--line);padding:5px 8px;text-align:left;vertical-align:top}th{background:var(--soft);font-weight:600}td.r,th.r{text-align:right;font-variant-numeric:tabular-nums}
code{background:var(--soft);padding:1px 4px;border-radius:3px;font-size:13px}ol.refs li{margin-bottom:4px}
@media print{nav{display:none}h2{break-before:page}figure,table{break-inside:avoid}}
"
html <- paste0(
  '<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">',
  '<title>auto16s report</title><style>', css, '</style></head><body><main>',
  '<h1>16S rRNA amplicon analysis report</h1>',
  sprintf('<p class="sub">auto16s %s &middot; run %s &middot; %s</p>', esc(opt$pipeline_version), esc(opt$run_name), format(Sys.time(), "%Y-%m-%d %H:%M")),
  '<nav><a href="#summary">1 Summary</a><a href="#samples">2 Samples</a><a href="#methods">3 Methods</a><a href="#results">4 Results</a><a href="#files">5 Files</a><a href="#software">6 Software and references</a></nav>',
  '<h2 id="summary">1 Summary</h2>', summary_html,
  '<h2 id="samples">2 Samples</h2>',
  sprintf("<p>%d samples in %d groups (column <code>%s</code>): %s.</p>", n_samples, length(groups), esc(gcol), esc(grp_txt)),
  table_html(samp, "Samples, metadata from the samplesheet, and read pairs at the start and the end of processing.", max_rows = 500),
  '<h2 id="methods">3 Methods</h2>', methods_html,
  table_html(param_df, "Values used in this run (chosen automatically unless set by the user)."),
  '<h2 id="results">4 Results</h2>', results_html,
  '<h2 id="files">5 Files</h2><p>All results are in the output folder; every figure is also available as PNG and PDF, and every number shown here as a table.</p>',
  table_html(files_df),
  '<h2 id="software">6 Software and references</h2>', table_html(sw_df),
  '<ol class="refs">', paste0("<li>", refs, "</li>", collapse = ""), '</ol>',
  '</main></body></html>')
# number figures and tables in page order
number <- function(x, tag) { parts <- strsplit(x, tag, fixed = TRUE)[[1]]; n <- length(parts) - 1
  if (n < 1) return(x); paste0(paste0(parts[-length(parts)], seq_len(n), collapse = ""), parts[length(parts)]) }
html <- number(number(html, "@FIG@"), "@TAB@")
writeLines(html, opt$out, useBytes = TRUE)
message(sprintf("Report written: %s (%.1f MB, %d figures, %d tables)", opt$out, file.size(opt$out) / 1e6, fig_no, tab_no))
