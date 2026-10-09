#!/usr/bin/env Rscript
# Alpha diversity on rarefied data: Observed, Shannon, Simpson (1-D)
#   2 groups: Wilcoxon rank-sum | > 2 groups: Kruskal-Wallis + pairwise Wilcoxon (BH)
.script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(.script)), "utils.R"))
opt <- parse_args()
suppressPackageStartupMessages({ library(phyloseq); library(ggplot2) })

ps   <- readRDS(opt$ps)
gcol <- opt$group_col
measures <- c("Observed", "Shannon", "Simpson")

alpha <- estimate_richness(ps, measures = measures)
rownames(alpha) <- sample_names(ps)                    # estimate_richness alters names (make.names)
alpha <- data.frame(sample = rownames(alpha), group = as.character(sample_data(ps)[[gcol]]),
                    alpha[, measures], check.names = FALSE)
write_tsv(alpha, "alpha_diversity.tsv")

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
}
tests <- do.call(rbind, res)
tests$p <- signif(tests$p, 4); tests$p_adj <- signif(tests$p_adj, 4)
write_tsv(tests, "alpha_diversity_tests.tsv")
write_mqc_table(tests, "alpha_diversity_tests_mqc.tsv", "alpha_tests", "Alpha diversity tests",
                "Rarefied data. Two groups: Wilcoxon rank-sum; more groups: Kruskal-Wallis and pairwise Wilcoxon with BH correction.")

long <- do.call(rbind, lapply(measures, function(m) data.frame(sample = alpha$sample, group = alpha$group,
                                                               measure = label[m], value = alpha[[m]])))
long$measure <- factor(long$measure, levels = label)
p <- ggplot(long, aes(group, value, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.35, width = 0.6) +
  geom_point(aes(colour = group), position = position_jitter(width = 0.12, seed = 1), size = 2) +
  facet_wrap(~ measure, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = group_palette(glev)) + scale_colour_manual(values = group_palette(glev)) +
  labs(title = "Alpha diversity",
       subtitle = sprintf("Rarefied to %s reads per sample; n = %s",
                          format(unique(sample_sums(ps))[1], big.mark = ","),
                          paste(names(table(alpha$group)), table(alpha$group), sep = ": ", collapse = ", ")),
       x = NULL, y = NULL) +
  theme_amp() + theme(legend.position = "none", axis.text.x = element_text(angle = if (k > 3) 30 else 0, hjust = if (k > 3) 1 else 0.5))
save_fig(p, "alpha_diversity", 3 + 2.2 * length(measures) + 0.4 * k, 4.8, mqc = TRUE)

write_session("alpha")
