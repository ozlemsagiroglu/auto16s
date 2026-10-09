#!/usr/bin/env Rscript
# Taxonomic filtering, metadata from the samplesheet, phyloseq object, rarefaction
.script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(.script)), "utils.R"))
opt <- parse_args()
suppressPackageStartupMessages({ library(phyloseq); library(vegan); library(ggplot2) })
seed <- as.integer(opt$seed)
gcol <- opt$group_col

seqtab <- readRDS(opt$seqtab)
taxa   <- readRDS(opt$taxa)
stopifnot(identical(colnames(seqtab), rownames(taxa)))

# ---------- taxonomic filtering ----------
is_na <- function(x) is.na(x) | x == ""
no_kp  <- is_na(taxa[, "Kingdom"]) | is_na(taxa[, "Phylum"])
euk    <- !no_kp & taxa[, "Kingdom"] == "Eukaryota"
chloro <- !no_kp & !euk & !is_na(taxa[, "Order"])  & taxa[, "Order"]  == "Chloroplast"
mito   <- !no_kp & !euk & !is_na(taxa[, "Family"]) & taxa[, "Family"] == "Mitochondria"
keep   <- !(no_kp | euk | chloro | mito)
tot <- colSums(seqtab)
filt_summary <- data.frame(
  category = c("kept", "no Kingdom/Phylum", "Eukaryota", "Chloroplast", "Mitochondria"),
  ASVs  = c(sum(keep), sum(no_kp), sum(euk), sum(chloro), sum(mito)),
  reads = c(sum(tot[keep]), sum(tot[no_kp]), sum(tot[euk]), sum(tot[chloro]), sum(tot[mito])))
filt_summary$reads_pct <- round(100 * filt_summary$reads / sum(tot), 2)
write_tsv(filt_summary, "taxonomic_filtering.tsv")
write_mqc_table(filt_summary, "taxonomic_filtering_mqc.tsv", "taxonomic_filtering", "ASV taxonomic filtering",
                "ASVs and reads removed as unassigned, eukaryotic, chloroplast or mitochondrial.")

seqtab.clean <- as.matrix(seqtab[, keep, drop = FALSE])
seqtab.clean[is.na(seqtab.clean)] <- 0
storage.mode(seqtab.clean) <- "integer"
taxa.clean <- taxa[keep, , drop = FALSE]

# short ids; sequences to FASTA (order: most abundant first)
ord  <- order(colSums(seqtab.clean), decreasing = TRUE)
seqtab.clean <- seqtab.clean[, ord, drop = FALSE]; taxa.clean <- taxa.clean[ord, , drop = FALSE]
seqs <- colnames(seqtab.clean); ids <- paste0("ASV_", seq_along(seqs))
colnames(seqtab.clean) <- ids; rownames(taxa.clean) <- ids
writeLines(paste0(">", ids, "\n", seqs), "asv_sequences.fasta")

# ---------- metadata = samplesheet minus FASTQ paths ----------
md <- read.csv(opt$samplesheet, check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character",
               na.strings = c("", "NA"))
md <- md[, setdiff(colnames(md), c("fastq_1", "fastq_2")), drop = FALSE]
rownames(md) <- md$sample
samples <- rownames(seqtab.clean)
md <- md[samples, , drop = FALSE]

ps <- phyloseq(otu_table(seqtab.clean, taxa_are_rows = FALSE), tax_table(taxa.clean), sample_data(md))
empty <- sample_names(ps)[sample_sums(ps) == 0]
if (length(empty)) {
  message("WARNING: 0 reads after taxonomic filtering, sample removed: ", paste(empty, collapse = ", "))
  ps <- prune_samples(sample_sums(ps) > 0, ps)
}
saveRDS(ps, "phyloseq_object.rds")

otu <- t(as(otu_table(ps), "matrix"))
write_tsv(data.frame(ASV = rownames(otu), otu, check.names = FALSE), "asv_table.tsv")
write_tsv(md, "sample_metadata.tsv")
tt <- as.data.frame(as(tax_table(ps), "matrix"), stringsAsFactors = FALSE)
write_tsv(data.frame(ASV = rownames(tt), tt, sequence = seqs[match(rownames(tt), ids)], check.names = FALSE),
          "asv_taxonomy.tsv")

# ---------- rarefaction depth ----------
# Default: the smallest depth among samples that are not outliers on the low side
# (>= 1,000 reads and >= 10% of the median). Lower samples are left out of alpha/beta only.
depths <- sample_sums(ps)
if (!is.na(opt$rarefy_depth)) {
  depth <- as.integer(opt$rarefy_depth); depth_rule <- "set by user"
} else {
  ok <- depths >= max(1000, 0.1 * median(depths))
  depth <- if (any(ok)) min(depths[ok]) else min(depths)
  depth_rule <- "smallest depth among samples with >= 1,000 reads and >= 10% of the median"
}
excluded <- names(depths)[depths < depth]
message(sprintf("Rarefaction depth %d (%s). Excluded from alpha/beta: %s", depth, depth_rule,
                if (length(excluded)) paste(excluded, collapse = ", ") else "none"))

ps_rare <- rarefy_even_depth(ps, sample.size = depth, rngseed = seed, replace = FALSE, trimOTUs = TRUE, verbose = FALSE)
saveRDS(ps_rare, "phyloseq_object_rarefied.rds")
otu_r <- t(as(otu_table(ps_rare), "matrix"))
write_tsv(data.frame(ASV = rownames(otu_r), otu_r, check.names = FALSE), "asv_table_rarefied.tsv")

grp_n <- table(sample_data(ps_rare)[[gcol]])
if (any(grp_n < 3)) message("WARNING: groups with < 3 samples after rarefaction: ",
                            paste(names(grp_n), grp_n, sep = "=", collapse = ", "))

summ <- data.frame(sample = names(depths), group = md[names(depths), gcol], reads = as.integer(depths),
                   in_alpha_beta = !names(depths) %in% excluded)
write_tsv(summ, "sample_depths.tsv")
write_mqc_table(summ, "sample_depths_mqc.tsv", "sample_depths", "Sample depth and rarefaction",
                sprintf("Reads after taxonomic filtering. Rarefaction depth = %d (%s).", depth, depth_rule))

# ---------- plots ----------
gpal <- group_palette(sort(unique(md[[gcol]])))
summ$sample <- factor(summ$sample, levels = summ$sample[order(summ$reads)])
p <- ggplot(summ, aes(reads, sample, fill = group)) + geom_col(width = 0.75) +
  geom_vline(xintercept = depth, linetype = 2) +
  annotate("text", x = depth, y = 0.6, label = paste0(" depth = ", format(depth, big.mark = ",")),
           hjust = 0, vjust = 0, size = 3.2) +
  scale_fill_manual(values = gpal) + scale_x_continuous(labels = scales::comma, expand = expansion(c(0, 0.05))) +
  labs(title = "Reads per sample after taxonomic filtering",
       subtitle = if (length(excluded)) paste("Below the depth, excluded from alpha/beta:", paste(excluded, collapse = ", ")) else "All samples kept for alpha/beta",
       x = "Reads", y = NULL, fill = gcol) + theme_amp() +
  theme(axis.text.y = element_text(size = if (nrow(summ) > 60) 5 else 8))
save_fig(p, "sample_depths_plot", 8, max(4, 0.16 * nrow(summ) + 2), mqc = TRUE)

rc <- rarecurve(as(otu_table(ps), "matrix"), step = max(1, floor(max(depths) / 100)), tidy = TRUE)
colnames(rc) <- c("sample", "reads", "ASVs")
rc$group <- md[as.character(rc$sample), gcol]
write_tsv(rc, "rarefaction_curves.tsv")
p <- ggplot(rc, aes(reads, ASVs, group = sample, colour = group)) + geom_line(alpha = 0.8) +
  geom_vline(xintercept = depth, linetype = 2) +
  scale_colour_manual(values = gpal) + scale_x_continuous(labels = scales::comma) +
  labs(title = "Rarefaction curves", subtitle = "Curves should flatten before the dashed line (rarefaction depth)",
       x = "Reads sampled", y = "Observed ASVs", colour = gcol) + theme_amp()
save_fig(p, "rarefaction_curves", 8, 5.5, mqc = TRUE)

write_session("phyloseq")
