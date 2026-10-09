#!/usr/bin/env Rscript
# Beta diversity on rarefied data: Bray-Curtis PCoA, PERMANOVA (adonis2), betadisper,
# pairwise PERMANOVA (BH) when there are > 2 groups
.script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(.script)), "utils.R"))
opt <- parse_args()
suppressPackageStartupMessages({ library(phyloseq); library(vegan); library(ggplot2) })

ps    <- readRDS(opt$ps)
gcol  <- opt$group_col
nperm <- as.integer(opt$permutations)
seed  <- as.integer(opt$seed)

d    <- phyloseq::distance(ps, method = "bray")
grp  <- factor(as.character(sample_data(ps)[[gcol]]))
glev <- levels(grp); gpal <- group_palette(glev)

# ---------- tests ----------
set.seed(seed); ad <- adonis2(d ~ grp, permutations = nperm)
set.seed(seed); bd_fit <- betadisper(d, grp); bd <- permutest(bd_fit, permutations = nperm)
perm_p <- ad$`Pr(>F)`[1]; perm_r2 <- ad$R2[1]; disp_p <- bd$tab$`Pr(>F)`[1]

interpretation <- if (perm_p < 0.05 && disp_p >= 0.05) {
  "Groups differ in community composition; within-group spread is similar, so the PERMANOVA result reflects a shift in location."
} else if (perm_p < 0.05 && disp_p < 0.05) {
  "PERMANOVA is significant but within-group spread also differs (betadisper). The difference may be partly or wholly due to dispersion; interpret with caution and check the PCoA."
} else if (perm_p >= 0.05 && disp_p < 0.05) {
  "No significant difference in composition, but groups differ in within-group spread (one group is more heterogeneous)."
} else {
  "No significant difference in composition or in within-group spread."
}

stats <- data.frame(test = c("PERMANOVA (adonis2)", "betadisper (permutest)"),
                    statistic = c(sprintf("R2 = %.3f, F = %.2f", perm_r2, ad$F[1]), sprintf("F = %.2f", bd$tab$F[1])),
                    p = c(perm_p, disp_p), permutations = nperm)
pairwise <- NULL
if (length(glev) > 2) {
  pr <- combn(glev, 2, simplify = FALSE)
  pairwise <- do.call(rbind, lapply(pr, function(g) {
    s  <- grp %in% g
    dd <- as.dist(as.matrix(d)[s, s])
    set.seed(seed); a <- adonis2(dd ~ droplevels(grp[s]), permutations = nperm)
    data.frame(comparison = paste(g, collapse = " vs "), R2 = round(a$R2[1], 4), F = round(a$F[1], 3), p = a$`Pr(>F)`[1])
  }))
  pairwise$p_adj <- p.adjust(pairwise$p, "BH")
  write_tsv(pairwise, "permanova_pairwise.tsv")
  write_mqc_table(pairwise, "permanova_pairwise_mqc.tsv", "permanova_pairwise", "Pairwise PERMANOVA (Bray-Curtis)",
                  "Pairwise adonis2 between groups, BH-adjusted.")
}
write_tsv(stats, "beta_diversity_tests.tsv")
write_mqc_table(rbind(stats, data.frame(test = "Interpretation", statistic = interpretation, p = NA, permutations = NA)),
                "beta_diversity_tests_mqc.tsv", "beta_tests", "Beta diversity tests (Bray-Curtis)",
                "Rarefied data. PERMANOVA tests location and spread together; betadisper tests spread only.")
writeLines(interpretation, "beta_diversity_interpretation.txt")

# ---------- PCoA ----------
pc <- cmdscale(d, k = 2, eig = TRUE)
ve <- 100 * pc$eig[1:2] / sum(pc$eig[pc$eig > 0])
coords <- data.frame(sample = sample_names(ps), group = grp, PCoA1 = pc$points[, 1], PCoA2 = pc$points[, 2])
write_tsv(coords, "pcoa_bray_coords.tsv")
dm <- as.matrix(d); write_tsv(data.frame(sample = rownames(dm), dm, check.names = FALSE), "bray_curtis_distance.tsv")
cent <- aggregate(cbind(PCoA1, PCoA2) ~ group, coords, mean)

p <- ggplot(coords, aes(PCoA1, PCoA2, colour = group)) +
  geom_hline(yintercept = 0, colour = "grey90") + geom_vline(xintercept = 0, colour = "grey90")
if (all(table(grp) >= 4)) p <- p + stat_ellipse(aes(fill = group), geom = "polygon", alpha = 0.08, level = 0.95, show.legend = FALSE)
p <- p + geom_point(size = 2.6, alpha = 0.9) +
  geom_point(data = cent, shape = 4, size = 5, stroke = 1.5, show.legend = FALSE) +
  scale_colour_manual(values = gpal) + scale_fill_manual(values = gpal) + coord_equal() +
  labs(title = "PCoA, Bray-Curtis",
       subtitle = sprintf("PERMANOVA R2 = %.3f, p %s | betadisper p %s | x = group centroid",
                          perm_r2, fmt_p(perm_p), fmt_p(disp_p)),
       x = sprintf("PCoA1 (%.1f%%)", ve[1]), y = sprintf("PCoA2 (%.1f%%)", ve[2]), colour = gcol) +
  theme_amp()
save_fig(p, "pcoa_bray", 8, 6.5, mqc = TRUE)

# ---------- dispersion ----------
disp <- data.frame(sample = sample_names(ps), group = grp, distance = bd_fit$distances)
write_tsv(disp, "betadisper_distances.tsv")
p <- ggplot(disp, aes(group, distance, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.35, width = 0.6) +
  geom_point(aes(colour = group), position = position_jitter(width = 0.12, seed = 1), size = 2) +
  scale_fill_manual(values = gpal) + scale_colour_manual(values = gpal) +
  labs(title = "Within-group dispersion (betadisper)", subtitle = sprintf("Distance to group centroid; permutest p %s", fmt_p(disp_p)),
       x = NULL, y = "Distance to centroid") + theme_amp() + theme(legend.position = "none")
save_fig(p, "betadisper_bray", 3 + 1.1 * length(glev), 4.8)

write_session("beta")
