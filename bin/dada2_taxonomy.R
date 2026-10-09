#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
opt  <- setNames(as.list(args[seq(2, length(args), 2)]), sub("^--", "", args[seq(1, length(args), 2)]))

suppressPackageStartupMessages(library(dada2))
set.seed(as.integer(opt$seed))   # assignTaxonomy bootstraps are random

seqtab <- readRDS(opt$seqtab)
taxa <- assignTaxonomy(seqtab, opt$db, minBoot = as.integer(opt$min_boot),
                       multithread = as.integer(opt$cpus))
saveRDS(taxa, "taxonomy.rds")
write.table(data.frame(sequence = rownames(taxa), taxa, check.names = FALSE),
            "taxonomy.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
