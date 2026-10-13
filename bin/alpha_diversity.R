#!/usr/bin/env Rscript
# Alpha diversity on rarefied data: Observed, Shannon, Simpson (1-D)
#   2 groups: Wilcoxon rank-sum | > 2 groups: Kruskal-Wallis + pairwise Wilcoxon (BH)
#   with covariates: additionally a linear model of the group adjusted for the covariates,
#   and (if requested) the test of the group x covariate interaction
.script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(.script)), "utils.R"))
opt <- parse_args()
suppressPackageStartupMessages({ library(phyloseq); library(ggplot2) })

gcol <- opt$group_col
cmp  <- comparison_setup(readRDS(opt$ps), opt)
if (!cmp$ok) skip_comparison(cmp, "alpha", gcol)
ps <- cmp$ps; on <- cmp$on
measures <- c("Observed", "Shannon", "Simpson")

alpha <- estimate_richness(ps, measures = measures)
rownames(alpha) <- sample_names(ps)                    # estimate_richness alters names (make.names)
md <- data.frame(sample_data(ps), check.names = FALSE)
alpha <- data.frame(sample = rownames(alpha), group = as.character(md[[gcol]]),
                    md[, cmp$covs, drop = FALSE], alpha[, measures], check.names = FALSE)
write_tsv(alpha, on("alpha_diversity.tsv"))

glev <- sort(unique(alpha$group)); k <- length(glev)
res <- list(); label <- setNames(character(length(measures)), measures)
for (m in measures) {
  if (k == 2) {
    t <- suppressWarnings(wilcox.test(alpha[[m]] ~ alpha$group))
    res[[length(res) + 1]] <- data.frame(measure = m, test = "Wilcoxon rank-sum", comparison = paste(glev, collapse = " vs "),
                                         statistic = unname(t$statistic), p = t$p.value, p_adj = NA)
    label[m] <- sprintf("%s\nWilcoxon p %s", m, fmt_p(t$p.value))
  } else {
    t <- kruskal.test(alpha[[m]] ~ alpha$group)
    res[[length(res) + 1]] <- data.frame(measure = m, test = "Kruskal-Wallis", comparison = "all groups",
                                         statistic = unname(t$statistic), p = t$p.value, p_adj = NA)
    label[m] <- sprintf("%s\nKruskal-Wallis p %s", m, fmt_p(t$p.value))
    pw <- suppressWarnings(pairwise.wilcox.test(alpha[[m]], alpha$group, p.adjust.method = "none"))$p.value
    pairs <- which(!is.na(pw), arr.ind = TRUE)
    praw  <- pw[pairs]
    res[[length(res) + 1]] <- data.frame(measure = m, test = "pairwise Wilcoxon",
                                         comparison = paste(rownames(pw)[pairs[, 1]], "vs", colnames(pw)[pairs[, 2]]),
                                         statistic = NA, p = praw, p_adj = p.adjust(praw, "BH"))
  }
  # adjusted for covariates: linear model, F test of the group term given the covariates
  covs <- cmp$covs[vapply(cmp$covs, function(cv) length(unique(alpha[[cv]])) > 1, TRUE)]
  if (length(covs)) {
    dat <- data.frame(y = alpha[[m]], group = factor(alpha$group), alpha[, covs, drop = FALSE], check.names = FALSE)
    names(dat)[-(1:2)] <- make.names(covs)
    cv <- paste(make.names(covs), collapse = " + ")
    fit0 <- lm(as.formula(paste("y ~", cv)), dat); fit1 <- lm(as.formula(paste("y ~", cv, "+ group")), dat)
    a <- anova(fit0, fit1)
    res[[length(res) + 1]] <- data.frame(measure = m, test = sprintf("linear model, adjusted for %s", paste(covs, collapse = ", ")),
                                         comparison = paste(glev, collapse = " vs "), statistic = a$F[2], p = a$`Pr(>F)`[2], p_adj = NA)
    if (cmp$interaction) {
      fit2 <- lm(as.formula(paste("y ~ group *", make.names(covs[1]))), dat)
      a2 <- anova(lm(as.formula(paste("y ~ group +", make.names(covs[1]))), dat), fit2)
      res[[length(res) + 1]] <- data.frame(measure = m, test = sprintf("interaction %s x %s (linear model)", gcol, covs[1]),
                                           comparison = sprintf("does the %s difference depend on %s?", gcol, covs[1]),
                                           statistic = a2$F[2], p = a2$`Pr(>F)`[2], p_adj = NA)
    }
  }
}
tests <- do.call(rbind, res)
tests$p <- signif(tests$p, 4); tests$p_adj <- signif(tests$p_adj, 4); tests$statistic <- signif(tests$statistic, 4)
write_tsv(tests, on("alpha_diversity_tests.tsv"))
write_mqc_table(tests, on("alpha_diversity_tests_mqc.tsv"), "alpha_tests", "Alpha diversity tests",
                "Rarefied data. Two groups: Wilcoxon rank-sum; more groups: Kruskal-Wallis and pairwise Wilcoxon with BH correction.")

long <- do.call(rbind, lapply(measures, function(m) data.frame(sample = alpha$sample, group = alpha$group,
                                                               measure = label[m], value = alpha[[m]])))
long$measure <- factor(long$measure, levels = label)
p <- ggplot(long, aes(group, value, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.35, width = 0.6) +
  geom_point(aes(colour = group), position = position_jitter(width = 0.12, seed = 1), size = 2) +
  facet_wrap(~ measure, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = group_palette(glev)) + scale_colour_manual(values = group_palette(glev)) +
  labs(title = if (nzchar(cmp$label)) paste("Alpha diversity:", cmp$label) else "Alpha diversity",
       subtitle = sprintf("Rarefied to %s reads per sample; n = %s",
                          format(unique(sample_sums(ps))[1], big.mark = ","),
                          paste(names(table(alpha$group)), table(alpha$group), sep = ": ", collapse = ", ")),
       x = NULL, y = NULL) +
  theme_amp() + theme(legend.position = "none")
fig_w <- 3 + 2.2 * length(measures) + 0.4 * k
p <- p + x_text_fit(glev, (fig_w - 1) / length(measures) / k)
save_fig(p, on("alpha_diversity"), fig_w, 4.8, mqc = TRUE)

write_session(on("alpha"))
