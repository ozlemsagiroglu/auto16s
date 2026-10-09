#!/usr/bin/env Rscript
# Genus-level differential abundance and prevalence with MaAsLin 3 (unrarefied counts).
# MaAsLin 3 fits two models per genus:
#   abundance  - is the genus more/less abundant where it is present? (linear model on log relative abundance)
#   prevalence - is the genus present in more/fewer samples?         (logistic model on presence/absence)
.script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(.script)), "utils.R"))
opt <- parse_args()
suppressPackageStartupMessages({ library(phyloseq); library(maaslin3); library(ggplot2) })

ps   <- readRDS(opt$ps)
gcol <- opt$group_col
qcut <- as.numeric(opt$qval)

# ---------- genus table ----------
ps_genus <- tax_glom(ps, taxrank = "Genus", NArm = FALSE)
tt <- as.data.frame(as(tax_table(ps_genus), "matrix"), stringsAsFactors = FALSE)
fixed <- fix_tax(as(tax_table(ps_genus), "matrix"))
label <- ifelse(is.na(tt$Genus) | tt$Genus == "", fixed[, "Genus"], tt$Genus)   # "Lachnospiraceae (unclassified)"
label <- make.unique(label, sep = " #")
feat <- make.names(label, unique = TRUE)        # syntactic ids, mapped back to genus names below
lookup <- data.frame(feature = feat, genus = label, resolved = !(is.na(tt$Genus) | tt$Genus == ""),
                     family = tt$Family, phylum = tt$Phylum, stringsAsFactors = FALSE)
taxa_names(ps_genus) <- feat

otu <- as(otu_table(ps_genus), "matrix"); if (taxa_are_rows(ps_genus)) otu <- t(otu)   # samples x genera
md  <- data.frame(sample_data(ps_genus), check.names = FALSE, stringsAsFactors = FALSE)

lv  <- sort(unique(md[[gcol]]))
ref <- if (!is.na(opt$reference)) opt$reference else lv[1]
if (!ref %in% lv) stop("--maaslin_reference '", ref, "' is not a level of ", gcol, " (", paste(lv, collapse = ", "), ")")
md[[gcol]] <- factor(md[[gcol]], levels = c(ref, setdiff(lv, ref)))
message("Reference level: ", ref)

# ---------- MaAsLin 3 ----------
# Only arguments known to the installed version are passed (they differ between releases).
args <- list(input_data = as.data.frame(otu), input_metadata = md, output = "maaslin3_output",
             fixed_effects = gcol, reference = paste0(gcol, ",", ref),
             normalization = "TSS", transform = "LOG",
             min_prevalence = as.numeric(opt$min_prevalence), max_significance = qcut,
             plot_summary_plot = FALSE, plot_associations = FALSE, cores = 1, verbosity = "WARN")
args <- args[names(args) %in% names(formals(maaslin3::maaslin3))]
fit  <- do.call(maaslin3::maaslin3, args)

take <- function(x, model) {
  if (is.null(x) || is.null(x$results) || !nrow(x$results)) return(NULL)
  r <- as.data.frame(x$results); r$model <- model; r
}
res <- rbind(take(fit$fit_data_abundance, "abundance"), take(fit$fit_data_prevalence, "prevalence"))
if (is.null(res)) stop("MaAsLin 3 returned no results")
res <- merge(res, lookup, by = "feature", all.x = TRUE, sort = FALSE)
res <- res[res$metadata == gcol, ]
res$comparison  <- paste(res$value, "vs", ref)
res$significant <- !is.na(res$qval_individual) & res$qval_individual < qcut & is.na(res$error)
res <- res[order(res$model, res$qval_individual), ]

cols <- intersect(c("genus", "family", "phylum", "resolved", "model", "comparison", "coef", "stderr",
                    "pval_individual", "qval_individual", "qval_joint", "N", "N_not_zero", "error", "significant"),
                  colnames(res))
write.csv(res[, cols], "maaslin3_all_results.csv", row.names = FALSE)
sig <- res[res$significant, ]
write.csv(sig[, cols], "maaslin3_significant.csv", row.names = FALSE)
message(sprintf("%d genera tested; significant (q < %g): %d abundance, %d prevalence",
                length(unique(res$feature)), qcut, sum(sig$model == "abundance"), sum(sig$model == "prevalence")))

mqc <- if (nrow(sig)) sig[, c("genus", "model", "comparison", "coef", "qval_individual")] else
  data.frame(genus = "none", model = NA, comparison = NA, coef = NA, qval_individual = NA)
mqc$coef <- round(mqc$coef, 3); mqc$qval_individual <- signif(mqc$qval_individual, 3)
write_mqc_table(mqc, "maaslin3_significant_mqc.tsv", "maaslin3", "Differential abundance and prevalence (MaAsLin 3, genus)",
                sprintf("Significant genera (q < %g). abundance: more/less abundant where present; prevalence: present in more/fewer samples. coef > 0: higher than in the reference group (%s).", qcut, ref))

model_lab <- c(abundance = "Abundance (where present)", prevalence = "Prevalence (present / absent)")

# ---------- volcano ----------
ok <- res[is.na(res$error) & !is.na(res$qval_individual), ]
ok$direction <- factor(ifelse(!ok$significant, "n.s.", ifelse(ok$coef > 0, "higher", "lower")),
                       levels = c("higher", "lower", "n.s."))
ok$panel <- factor(model_lab[ok$model], levels = model_lab)
lab <- head(ok[ok$significant, ], 20)
p <- ggplot(ok, aes(coef, -log10(qval_individual), colour = direction)) +
  geom_hline(yintercept = -log10(qcut), linetype = 2, colour = "grey50") + geom_vline(xintercept = 0, colour = "grey85") +
  geom_point(size = 2, alpha = 0.8) +
  geom_text(data = lab, aes(label = genus), size = 2.8, vjust = -0.7, show.legend = FALSE, check_overlap = TRUE) +
  scale_colour_manual(values = c(higher = "#D55E00", lower = "#0072B2", `n.s.` = "grey70"), drop = FALSE,
                      labels = c(higher = "higher than reference", lower = "lower than reference", `n.s.` = "not significant")) +
  scale_x_continuous(expand = expansion(mult = 0.2)) +
  facet_wrap(~ panel, scales = "free_x") +
  labs(title = "Differential abundance and prevalence (MaAsLin 3, genus level)",
       subtitle = sprintf("%s | reference: %s | dashed line: q = %g", unique(ok$comparison)[1], ref, qcut),
       x = "Coefficient", y = "-log10(q)", colour = NULL) +
  theme_amp() + theme(legend.position = "bottom")
save_fig(p, "maaslin3_volcano", 10, 6)

if (nrow(sig)) {
  # ---------- coefficient plot ----------
  top <- head(sig[order(-abs(sig$coef)), ], 40)
  top$panel <- factor(model_lab[top$model], levels = model_lab)
  top$genus <- factor(top$genus, levels = unique(top$genus[order(top$coef)]))
  p <- ggplot(top, aes(coef, genus, colour = coef > 0)) +
    geom_vline(xintercept = 0, colour = "grey70") +
    geom_errorbar(aes(xmin = coef - stderr, xmax = coef + stderr), width = 0.25, orientation = "y") +
    geom_point(size = 2.5) +
    scale_colour_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "#0072B2"), guide = "none") +
    facet_wrap(~ panel, scales = "free_x") +
    labs(title = "Significant genera",
         subtitle = sprintf("Coefficient +/- SE; positive = higher than %s (q < %g)", ref, qcut),
         x = "MaAsLin 3 coefficient", y = NULL) + theme_amp()
  save_fig(p, "maaslin3_coefficients", 10, max(4, 0.28 * length(unique(top$genus)) + 2), mqc = TRUE)

  rel <- 100 * sweep(otu, 1, rowSums(otu), "/")
  grp <- factor(md[rownames(rel), gcol], levels = levels(md[[gcol]]))

  # ---------- abundance hits: relative abundance where present ----------
  fa <- unique(head(sig$feature[sig$model == "abundance"], 12))
  if (length(fa)) {
    long <- do.call(rbind, lapply(fa, function(f) data.frame(
      genus = lookup$genus[lookup$feature == f], group = grp, abundance = rel[, f])))
    long <- long[long$abundance > 0, ]
    long$genus <- factor(long$genus, levels = lookup$genus[match(fa, lookup$feature)])
    p <- ggplot(long, aes(group, abundance, fill = group)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.35, width = 0.6) +
      geom_point(aes(colour = group), position = position_jitter(width = 0.12, seed = 1), size = 1.4) +
      scale_y_log10(labels = function(x) format(x, scientific = FALSE, drop0trailing = TRUE)) +
      scale_fill_manual(values = group_palette(levels(grp))) + scale_colour_manual(values = group_palette(levels(grp))) +
      facet_wrap(~ genus, scales = "free_y", ncol = 4) +
      labs(title = "Abundance hits: relative abundance where present",
           subtitle = sprintf("%d genera with the lowest q-value (max. 12); samples without the genus are not shown; log scale", length(fa)),
           x = NULL, y = "Relative abundance (%)") +
      theme_amp() + theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
    save_fig(p, "maaslin3_abundance_boxplots", 11, 2.6 * ceiling(length(fa) / 4) + 1.5, mqc = TRUE)
  }

  # ---------- prevalence hits: share of samples with the genus ----------
  fp <- unique(head(sig$feature[sig$model == "prevalence"], 12))
  if (length(fp)) {
    prev <- do.call(rbind, lapply(fp, function(f) {
      pr <- tapply(otu[, f] > 0, grp, mean)
      data.frame(genus = lookup$genus[lookup$feature == f], group = names(pr), prevalence = 100 * as.numeric(pr),
                 n = as.integer(table(grp)[names(pr)]))
    }))
    prev$genus <- factor(prev$genus, levels = lookup$genus[match(fp, lookup$feature)])
    prev$group <- factor(prev$group, levels = levels(grp))
    p <- ggplot(prev, aes(group, prevalence, fill = group)) + geom_col(width = 0.65) +
      geom_text(aes(label = sprintf("%.0f%%", prevalence)), vjust = -0.4, size = 3) +
      scale_fill_manual(values = group_palette(levels(grp))) +
      scale_y_continuous(limits = c(0, 110), breaks = seq(0, 100, 25)) +
      facet_wrap(~ genus, ncol = 4) +
      labs(title = "Prevalence hits: share of samples in which the genus is present",
           subtitle = sprintf("%d genera with the lowest q-value (max. 12)", length(fp)),
           x = NULL, y = "Samples with the genus (%)") +
      theme_amp() + theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
    save_fig(p, "maaslin3_prevalence_bars", 11, 2.6 * ceiling(length(fp) / 4) + 1.5, mqc = TRUE)
  }
}

write_session("maaslin3")
