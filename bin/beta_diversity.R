#!/usr/bin/env Rscript
# Beta diversity on rarefied data: Bray-Curtis PCoA, PERMANOVA (adonis2), betadisper,
# pairwise PERMANOVA (BH) when there are > 2 groups. With covariates, PERMANOVA tests the group
# adjusted for them (marginal test) and, if requested, the group x covariate interaction.
.script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(.script)), "utils.R"))
opt <- parse_args()
suppressPackageStartupMessages({ library(phyloseq); library(vegan); library(ggplot2) })

gcol  <- opt$group_col
cmp   <- comparison_setup(readRDS(opt$ps), opt)
if (!cmp$ok) skip_comparison(cmp, "beta", gcol)
ps    <- cmp$ps; on <- cmp$on
nperm <- as.integer(opt$permutations)
seed  <- as.integer(opt$seed)

d    <- phyloseq::distance(ps, method = "bray")
grp  <- factor(as.character(sample_data(ps)[[gcol]]))
glev <- levels(grp); gpal <- group_palette(glev)
md   <- data.frame(sample_data(ps), check.names = FALSE)
covs <- cmp$covs[vapply(cmp$covs, function(cv) length(unique(md[[cv]])) > 1, TRUE)]
cdat <- data.frame(grp = grp)
for (cv in covs) cdat[[make.names(cv)]] <- factor(as.character(md[[cv]]))

# ---------- tests ----------
set.seed(seed); ad <- adonis2(d ~ grp, data = cdat, permutations = nperm)
set.seed(seed); bd_fit <- betadisper(d, grp); bd <- permutest(bd_fit, permutations = nperm)
perm_p <- ad$`Pr(>F)`[1]; perm_r2 <- ad$R2[1]; perm_F <- ad$F[1]; disp_p <- bd$tab$`Pr(>F)`[1]
disp_ok <- all(table(grp) >= 3)        # with 2 samples a group's spread cannot be estimated
if (!disp_ok) disp_p <- NA
extra <- NULL; adj_txt <- ""
if (length(covs)) {
  # group effect adjusted for the covariates: marginal test of the group term
  set.seed(seed)
  am <- adonis2(as.formula(paste("d ~", paste(make.names(covs), collapse = " + "), "+ grp")), data = cdat,
                permutations = nperm, by = "margin")
  perm_p <- am["grp", "Pr(>F)"]; perm_r2 <- am["grp", "R2"]; perm_F <- am["grp", "F"]
  adj_txt <- sprintf(", adjusted for %s", paste(covs, collapse = ", "))
  extra <- data.frame(test = sprintf("PERMANOVA (adonis2), %s adjusted for %s", gcol, paste(covs, collapse = ", ")),
                      statistic = sprintf("R2 = %.3f, F = %.2f", perm_r2, perm_F), p = perm_p, permutations = nperm)
  for (cv in covs) {
    extra <- rbind(extra, data.frame(test = sprintf("PERMANOVA (adonis2), %s adjusted for %s", cv, gcol),
                                     statistic = sprintf("R2 = %.3f, F = %.2f", am[make.names(cv), "R2"], am[make.names(cv), "F"]),
                                     p = am[make.names(cv), "Pr(>F)"], permutations = nperm))
  }
  if (cmp$interaction) {
    set.seed(seed)
    ai <- adonis2(as.formula(paste("d ~ grp *", make.names(covs[1]))), data = cdat, permutations = nperm, by = "terms")
    term <- grep(":", rownames(ai), value = TRUE)[1]
    extra <- rbind(extra, data.frame(test = sprintf("PERMANOVA interaction %s x %s (does the %s difference depend on %s?)", gcol, covs[1], gcol, covs[1]),
                                     statistic = sprintf("R2 = %.3f, F = %.2f", ai[term, "R2"], ai[term, "F"]),
                                     p = ai[term, "Pr(>F)"], permutations = nperm))
  }
}

interpretation <- if (!disp_ok) "" else if (perm_p < 0.05 && disp_p >= 0.05) {
  "Groups differ in community composition; within-group spread is similar, so the PERMANOVA result reflects a shift in location."
} else if (perm_p < 0.05 && disp_p < 0.05) {
  "PERMANOVA is significant but within-group spread also differs (betadisper). The difference may be partly or wholly due to dispersion; interpret with caution and check the PCoA."
} else if (perm_p >= 0.05 && disp_p < 0.05) {
  "No significant difference in composition, but groups differ in within-group spread (one group is more heterogeneous)."
} else {
  "No significant difference in composition or in within-group spread."
}
if (!disp_ok) interpretation <- if (!is.na(perm_p) && perm_p < 0.05) "Groups differ in community composition (PERMANOVA). Within-group spread could not be tested (fewer than 3 samples in a group), so a contribution of dispersion cannot be ruled out." else
  "No significant difference in composition (PERMANOVA). Within-group spread could not be tested (fewer than 3 samples in a group)."

stats <- data.frame(test = c(if (length(covs)) "PERMANOVA (adonis2), unadjusted" else "PERMANOVA (adonis2)", "betadisper (permutest)"),
                    statistic = c(sprintf("R2 = %.3f, F = %.2f", ad$R2[1], ad$F[1]),
                                  if (disp_ok) sprintf("F = %.2f", bd$tab$F[1]) else "not tested: fewer than 3 samples in a group"),
                    p = c(ad$`Pr(>F)`[1], disp_p), permutations = nperm)
if (!is.null(extra)) stats <- rbind(stats[1, ], extra, stats[2, ])
stats$p <- signif(stats$p, 3)
pairwise <- NULL
if (length(glev) > 2) {
  pr <- combn(glev, 2, simplify = FALSE)
  pairwise <- do.call(rbind, lapply(pr, function(g) {
    s  <- grp %in% g
    dd <- as.dist(as.matrix(d)[s, s])
    gs <- droplevels(grp[s])
    set.seed(seed); a <- adonis2(dd ~ gs, permutations = nperm)
    data.frame(comparison = paste(g, collapse = " vs "), R2 = round(a$R2[1], 4), F = round(a$F[1], 3), p = a$`Pr(>F)`[1])
  }))
  pairwise$p_adj <- p.adjust(pairwise$p, "BH")
  write_tsv(pairwise, on("permanova_pairwise.tsv"))
  write_mqc_table(pairwise, on("permanova_pairwise_mqc.tsv"), "permanova_pairwise", "Pairwise PERMANOVA (Bray-Curtis)",
                  "Pairwise adonis2 between groups, BH-adjusted.")
}
write_tsv(stats, on("beta_diversity_tests.tsv"))
write_mqc_table(rbind(stats, data.frame(test = "Interpretation", statistic = interpretation, p = NA, permutations = NA)),
                on("beta_diversity_tests_mqc.tsv"), "beta_tests", "Beta diversity tests (Bray-Curtis)",
                "Rarefied data. PERMANOVA tests location and spread together; betadisper tests spread only.")
writeLines(interpretation, on("beta_diversity_interpretation.txt"))

# ---------- PCoA ----------
pc <- cmdscale(d, k = 2, eig = TRUE)
ve <- 100 * pc$eig[1:2] / sum(pc$eig[pc$eig > 0])
coords <- data.frame(sample = sample_names(ps), group = grp, PCoA1 = pc$points[, 1], PCoA2 = pc$points[, 2])
shape_col <- if (length(covs) && length(unique(md[[covs[1]]])) <= 6) covs[1] else NA
if (!is.na(shape_col)) coords$shape <- factor(as.character(md[[shape_col]]))
write_tsv(coords, on("pcoa_bray_coords.tsv"))
dm <- as.matrix(d); write_tsv(data.frame(sample = rownames(dm), dm, check.names = FALSE), on("bray_curtis_distance.tsv"))
cent <- aggregate(cbind(PCoA1, PCoA2) ~ group, coords, mean)

# Each sample is joined to its group centroid (the distance betadisper tests); a dashed line
# outlines the area covered by the samples of a group (convex hull, groups with >= 3 samples).
spokes <- merge(coords, setNames(cent, c("group", "c1", "c2")), by = "group")
hulls <- do.call(rbind, lapply(split(coords, coords$group), function(d) if (nrow(d) >= 3) d[chull(d$PCoA1, d$PCoA2), ] else NULL))
p <- ggplot(coords, aes(PCoA1, PCoA2, colour = group)) +
  geom_hline(yintercept = 0, colour = "grey90") + geom_vline(xintercept = 0, colour = "grey90")
if (!is.null(hulls) && nrow(hulls))
  p <- p + geom_polygon(data = hulls, aes(fill = group), alpha = 0.06, linetype = "dashed", linewidth = 0.5, show.legend = FALSE)
p <- p +
  geom_segment(data = spokes, aes(xend = c1, yend = c2), linewidth = 0.4, alpha = 0.5, show.legend = FALSE) +
  { if (!is.na(shape_col)) geom_point(aes(shape = shape), size = 2.8, alpha = 0.9) else geom_point(size = 2.6, alpha = 0.9) } +
  geom_point(data = cent, aes(fill = group), shape = 23, size = 4.2, colour = "white", stroke = 1, show.legend = FALSE) +
  scale_colour_manual(values = gpal) + scale_fill_manual(values = gpal) +
  labs(title = if (nzchar(cmp$label)) paste("PCoA, Bray-Curtis:", cmp$label) else "PCoA, Bray-Curtis",
       subtitle = sprintf("PERMANOVA%s R2 = %.3f, p %s | %s\nDiamond: group centroid; lines join each sample to it; dashed outline: area covered by the group",
                          adj_txt, perm_r2, fmt_p(perm_p),
                          if (is.na(disp_p)) "betadisper not tested (< 3 samples per group)" else paste("betadisper p", fmt_p(disp_p))),
       x = sprintf("PCoA1 (%.1f%%)", ve[1]), y = sprintf("PCoA2 (%.1f%%)", ve[2]), colour = gcol,
       shape = if (!is.na(shape_col)) shape_col else NULL) +
  theme_amp()
# Equal axis scaling (distances are comparable in both directions); figure size follows the data's shape
rx <- diff(range(coords$PCoA1)); ry <- diff(range(coords$PCoA2))
ratio <- if (rx > 0) min(max(ry / rx, 0.45), 1.4) else 1
p <- p + coord_equal(ratio = if (ry > 0 && rx > 0) ratio / (ry / rx) else 1)
legend_w <- 0.075 * max(nchar(c(gcol, glev, if (!is.na(shape_col)) c(shape_col, levels(coords$shape))))) + 0.9
panel_w <- 5.6
save_fig(p, on("pcoa_bray"), panel_w + legend_w + 0.9, panel_w * ratio + 1.9, mqc = TRUE)

# ---------- dispersion ----------
if (disp_ok) {
  disp <- data.frame(sample = sample_names(ps), group = grp, distance = bd_fit$distances)
  write_tsv(disp, on("betadisper_distances.tsv"))
  p <- ggplot(disp, aes(group, distance, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.35, width = 0.6) +
    geom_point(aes(colour = group), position = position_jitter(width = 0.12, seed = 1), size = 2) +
    scale_fill_manual(values = gpal) + scale_colour_manual(values = gpal) +
    labs(title = "Within-group dispersion (betadisper)", subtitle = if (is.na(disp_p)) "Distance to group centroid; not tested (< 3 samples per group)" else sprintf("Distance to group centroid; permutest p %s", fmt_p(disp_p)),
         x = NULL, y = "Distance to centroid") + theme_amp() + theme(legend.position = "none") +
    x_text_fit(glev, (2 + 1.1 * length(glev)) / length(glev))
  save_fig(p, on("betadisper_bray"), 3 + 1.1 * length(glev), 4.8, mqc = TRUE)
}

write_session(on("beta"))
