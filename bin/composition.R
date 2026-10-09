#!/usr/bin/env Rscript
# Taxonomic composition (relative abundance, unrarefied):
#   per-sample stacked bars, group-mean bars (Phylum/Family/Genus) and a genus heatmap
.script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(.script)), "utils.R"))
opt <- parse_args()
suppressPackageStartupMessages({ library(phyloseq); library(ggplot2) })

ps     <- readRDS(opt$ps)
gcol   <- opt$group_col
top_n  <- as.integer(opt$top_n)
levels <- c("Phylum", "Family", "Genus")
plural <- c(Phylum = "phyla", Family = "families", Genus = "genera")

tax_table(ps) <- tax_table(fix_tax(as(tax_table(ps), "matrix")))
grp  <- setNames(as.character(sample_data(ps)[[gcol]]), sample_names(ps))
glev <- sort(unique(grp))

# taxa x samples relative abundance at one rank
rel_at <- function(rank) {
  m <- as(otu_table(ps), "matrix"); if (!taxa_are_rows(ps)) m <- t(m)
  m <- rowsum(m, as.character(tax_table(ps)[rownames(m), rank]))
  sweep(m, 2, colSums(m), "/")
}

for (lvl in levels) {
  cnt <- as(otu_table(ps), "matrix"); if (!taxa_are_rows(ps)) cnt <- t(cnt)
  cnt <- rowsum(cnt, as.character(tax_table(ps)[rownames(cnt), lvl]))
  write_tsv(data.frame(taxon = rownames(cnt), cnt, check.names = FALSE), paste0(tolower(lvl), "_counts.tsv"))
  rel <- rel_at(lvl)
  write_tsv(data.frame(taxon = rownames(rel), rel, check.names = FALSE), paste0(tolower(lvl), "_relative_abundance.tsv"))

  top  <- names(sort(rowMeans(rel), decreasing = TRUE))[seq_len(min(top_n, nrow(rel)))]
  show <- rel[top, , drop = FALSE]
  if (nrow(rel) > length(top)) show <- rbind(show, Other = colSums(rel[!rownames(rel) %in% top, , drop = FALSE]))
  pal  <- c(taxa_palette(top), Other = "grey80")

  long <- data.frame(sample = rep(colnames(show), each = nrow(show)), taxon = rep(rownames(show), ncol(show)),
                     abundance = as.vector(show))
  long$group <- grp[long$sample]
  long$taxon <- factor(long$taxon, levels = rev(names(pal)))
  # within each group, order samples by the most abundant taxon overall
  sample_order <- colnames(show)[order(grp[colnames(show)], -show[top[1], ])]
  long$sample <- factor(long$sample, levels = sample_order)

  p <- ggplot(long, aes(sample, abundance, fill = taxon)) + geom_col(width = 0.92) +
    facet_grid(~ group, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = pal, breaks = names(pal)) +
    scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
    labs(title = sprintf("%s composition per sample", lvl),
         subtitle = sprintf("Top %d %s by mean relative abundance; the rest is 'Other'", length(top), plural[[lvl]]),
         x = NULL, y = "Relative abundance", fill = lvl) +
    theme_amp() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                                   size = if (ncol(show) > 60) 5 else 8),
                        panel.spacing = unit(0.3, "lines"))
  save_fig(p, paste0(tolower(lvl), "_barplot_samples"), max(8, 0.18 * ncol(show) + 4), 6,
           mqc = lvl == "Genus")

  gm <- aggregate(abundance ~ group + taxon, long, mean)
  p <- ggplot(gm, aes(group, abundance, fill = taxon)) + geom_col(width = 0.7) +
    scale_fill_manual(values = pal, breaks = names(pal)) +
    scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
    labs(title = sprintf("Mean %s composition per group", tolower(lvl)), x = gcol, y = "Mean relative abundance", fill = lvl) +
    theme_amp()
  save_fig(p, paste0(tolower(lvl), "_barplot_groups"), 3 + 1.2 * length(glev), 6, mqc = lvl == "Family")
}

# ---------- genus heatmap (top 25, log10 %) ----------
rel  <- rel_at("Genus")
top  <- names(sort(rowMeans(rel), decreasing = TRUE))[seq_len(min(25, nrow(rel)))]
hm   <- log10(100 * rel[top, , drop = FALSE] + 0.01)
taxon_order  <- if (nrow(hm) > 2) rownames(hm)[hclust(dist(hm))$order] else rownames(hm)
sample_order <- unlist(lapply(glev, function(g) {
  s <- names(grp)[grp == g]
  if (length(s) > 2) s[hclust(dist(t(hm[, s, drop = FALSE])))$order] else s
}))
long <- data.frame(taxon = rep(rownames(hm), ncol(hm)), sample = rep(colnames(hm), each = nrow(hm)), value = as.vector(hm))
long$group  <- grp[long$sample]
long$taxon  <- factor(long$taxon, levels = taxon_order)
long$sample <- factor(long$sample, levels = sample_order)
p <- ggplot(long, aes(sample, taxon, fill = value)) + geom_tile() +
  facet_grid(~ group, scales = "free_x", space = "free_x") +
  scale_fill_viridis_c(name = "Relative\nabundance (%)", breaks = -2:2, labels = c("0.01", "0.1", "1", "10", "100"),
                       limits = c(-2, 2), oob = scales::squish) +
  labs(title = sprintf("Top %d genera", nrow(hm)), subtitle = "log10 scale; genera clustered by abundance profile, samples clustered within group",
       x = NULL, y = NULL) +
  theme_amp() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = if (ncol(hm) > 60) 5 else 8),
                      panel.grid = element_blank(), panel.spacing = unit(0.3, "lines"))
save_fig(p, "genus_heatmap", max(8, 0.18 * ncol(hm) + 5), 7.5, mqc = TRUE)

write_session("composition")
