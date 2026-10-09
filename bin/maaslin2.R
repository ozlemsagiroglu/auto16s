#!/usr/bin/env Rscript
# Genus-level differential abundance with MaAsLin2 (unrarefied counts; TSS + LOG inside MaAsLin2)
.script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(.script)), "utils.R"))
opt <- parse_args()
suppressPackageStartupMessages({ library(phyloseq); library(Maaslin2); library(ggplot2) })

ps   <- readRDS(opt$ps)
gcol <- opt$group_col
qcut <- as.numeric(opt$qval)

# ---------- genus table ----------
ps_genus <- tax_glom(ps, taxrank = "Genus", NArm = FALSE)
tt <- as.data.frame(as(tax_table(ps_genus), "matrix"), stringsAsFactors = FALSE)
fixed <- fix_tax(as(tax_table(ps_genus), "matrix"))
label <- ifelse(is.na(tt$Genus) | tt$Genus == "", fixed[, "Genus"], tt$Genus)   # "Lachnospiraceae (unclassified)"
label <- make.unique(label, sep = " #")
# MaAsLin2 make.names()-mangles ids (Escherichia-Shigella, [Ruminococcus]_torques_group): use syntactic ids
feat <- make.names(label, unique = TRUE)
lookup <- data.frame(feature = feat, genus = label, resolved = !(is.na(tt$Genus) | tt$Genus == ""),
                     family = tt$Family, phylum = tt$Phylum, stringsAsFactors = FALSE)
taxa_names(ps_genus) <- feat

otu <- as(otu_table(ps_genus), "matrix"); if (taxa_are_rows(ps_genus)) otu <- t(otu)   # samples x genera
md  <- data.frame(sample_data(ps_genus), check.names = FALSE, stringsAsFactors = FALSE)

lv <- sort(unique(md[[gcol]]))
ref <- if (!is.na(opt$reference)) opt$reference else lv[1]
if (!ref %in% lv) stop("--maaslin_reference '", ref, "' is not a level of ", gcol, " (", paste(lv, collapse = ", "), ")")
message("Reference level: ", ref)

fit <- Maaslin2(input_data = as.data.frame(otu), input_metadata = md, output = "maaslin2_output",
                fixed_effects = gcol, reference = paste0(gcol, ",", ref),
                normalization = "TSS", transform = "LOG", min_prevalence = as.numeric(opt$min_prevalence),
                cores = as.integer(opt$cpus), plot_heatmap = FALSE, plot_scatter = FALSE)

res <- fit$results
names(res)[grepl("^N\\.not", names(res))] <- "N.not.0"     # "N.not.zero" in the R object, "N.not.0" in the TSV
res <- merge(res, lookup, by = "feature", all.x = TRUE, sort = FALSE)
res <- res[order(res$qval, res$pval), ]
res$comparison <- paste(res$value, "vs", ref)
res$significant <- !is.na(res$qval) & res$qval < qcut
out_cols <- c("genus", "family", "phylum", "resolved", "comparison", "coef", "stderr", "pval", "qval", "N", "N.not.0", "significant")
write.csv(res[, out_cols], "maaslin2_all_results.csv", row.names = FALSE)
sig <- res[res$significant, ]
write.csv(sig[, out_cols], "maaslin2_significant.csv", row.names = FALSE)
message(sprintf("%d genera tested, %d significant (q < %g)", length(unique(res$feature)), length(unique(sig$feature)), qcut))

mqc <- if (nrow(sig)) sig[, c("genus", "comparison", "coef", "qval")] else
  data.frame(genus = "none", comparison = NA, coef = NA, qval = NA)
mqc$coef <- round(mqc$coef, 3); mqc$qval <- signif(mqc$qval, 3)
write_mqc_table(mqc, "maaslin2_significant_mqc.tsv", "maaslin2", "Differential abundance (MaAsLin2, genus)",
                sprintf("Significant genera (q < %g). coef > 0: higher than in the reference group (%s).", qcut, ref))

# ---------- volcano ----------
res$direction <- ifelse(!res$significant, "n.s.", ifelse(res$coef > 0, "higher", "lower"))
res$direction <- factor(res$direction, levels = c("higher", "lower", "n.s."))
lab <- head(res[res$significant, ], 15)
p <- ggplot(res, aes(coef, -log10(qval), colour = direction)) +
  geom_hline(yintercept = -log10(qcut), linetype = 2, colour = "grey50") + geom_vline(xintercept = 0, colour = "grey85") +
  geom_point(size = 2, alpha = 0.8) +
  geom_text(data = lab, aes(label = genus), size = 2.8, vjust = -0.7, show.legend = FALSE, check_overlap = TRUE) +
  scale_colour_manual(values = c(higher = "#D55E00", lower = "#0072B2", `n.s.` = "grey70"), drop = FALSE,
                      labels = c(higher = "higher than reference", lower = "lower than reference", `n.s.` = "not significant")) +
  facet_wrap(~ comparison) +
  labs(title = "Differential abundance (MaAsLin2, genus level)",
       subtitle = sprintf("Reference: %s | dashed line: q = %g", ref, qcut), x = "Coefficient", y = "-log10(q)", colour = NULL) +
  theme_amp() + theme(legend.position = "bottom")
save_fig(p, "maaslin2_volcano", 8, 6)

if (nrow(sig)) {
  # ---------- coefficient plot ----------
  top <- head(sig[order(-abs(sig$coef)), ], 30)
  top$genus <- factor(top$genus, levels = unique(top$genus[order(top$coef)]))
  p <- ggplot(top, aes(coef, genus, colour = coef > 0)) +
    geom_vline(xintercept = 0, colour = "grey70") +
    geom_errorbar(aes(xmin = coef - stderr, xmax = coef + stderr), width = 0.25, orientation = "y") + geom_point(size = 2.5) +
    scale_colour_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "#0072B2"), guide = "none") +
    facet_wrap(~ comparison) +
    labs(title = "Significant genera", subtitle = sprintf("Coefficient +/- SE; positive = higher than %s (q < %g)", ref, qcut),
         x = "MaAsLin2 coefficient", y = NULL) + theme_amp()
  save_fig(p, "maaslin2_coefficients", 8, max(4, 0.28 * nrow(top) + 2), mqc = TRUE)

  # ---------- relative abundance of the top significant genera ----------
  feats <- unique(head(sig$feature[order(sig$qval)], 12))
  rel <- 100 * sweep(otu, 1, rowSums(otu), "/")
  long <- do.call(rbind, lapply(feats, function(f) data.frame(
    genus = lookup$genus[lookup$feature == f], group = md[rownames(rel), gcol], abundance = rel[, f])))
  long$genus <- factor(long$genus, levels = lookup$genus[match(feats, lookup$feature)])
  long$group <- factor(long$group, levels = lv)
  p <- ggplot(long, aes(group, abundance + 0.001, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.35, width = 0.6) +
    geom_point(aes(colour = group), position = position_jitter(width = 0.12, seed = 1), size = 1.4) +
    scale_y_log10(labels = function(x) format(x, scientific = FALSE, drop0trailing = TRUE)) +
    scale_fill_manual(values = group_palette(lv)) + scale_colour_manual(values = group_palette(lv)) +
    facet_wrap(~ genus, scales = "free_y", ncol = 4) +
    labs(title = "Relative abundance of significant genera", subtitle = sprintf("%d genera with the lowest q-value (max. 12); log scale (+0.001%%)", length(feats)),
         x = NULL, y = "Relative abundance (%)") +
    theme_amp() + theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
  save_fig(p, "maaslin2_boxplots", 11, 2.6 * ceiling(length(feats) / 4) + 1.5, mqc = TRUE)
}

write_session("maaslin2")
