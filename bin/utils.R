# Shared helpers, sourced by every R step:
#   source(file.path(dirname(.script_path()), "utils.R"))

parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) %% 2 != 0) stop("Arguments must be --key value pairs")
  keys <- sub("^--", "", args[seq(1, length(args), 2)])
  vals <- args[seq(2, length(args), 2)]
  vals[vals %in% c("null", "NULL", "NA", "")] <- NA
  setNames(as.list(vals), keys)
}

split_csv <- function(x) if (is.null(x) || is.na(x)) character(0) else trimws(strsplit(x, ",")[[1]])

# ---------- colours ----------
# Okabe-Ito: colour-blind safe, distinct for up to 8 groups
group_palette <- function(levels) {
  oi <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#000000")
  cols <- if (length(levels) <= length(oi)) oi[seq_along(levels)] else grDevices::hcl.colors(length(levels), "Dark 3")
  setNames(cols, levels)
}

taxa_palette <- function(taxa) {
  base <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B", "#E377C2", "#17BECF",
            "#BCBD22", "#AEC7E8", "#FFBB78", "#98DF8A", "#FF9896", "#C5B0D5", "#C49C94", "#F7B6D2",
            "#9EDAE5", "#DBDB8D", "#393B79", "#637939")
  cols <- if (length(taxa) <= length(base)) base[seq_along(taxa)] else grDevices::hcl.colors(length(taxa), "Spectral")
  setNames(cols, taxa)
}

# ---------- ggplot theme / saving ----------
theme_amp <- function(base_size = 11) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   strip.background = ggplot2::element_rect(fill = "grey95", colour = NA),
                   strip.text = ggplot2::element_text(face = "bold"),
                   plot.title = ggplot2::element_text(face = "bold"),
                   plot.subtitle = ggplot2::element_text(colour = "grey30"),
                   legend.key.size = ggplot2::unit(0.45, "cm"))
}

# x-axis labels: rotate them when they would not fit side by side (slot_in = inches available per label)
x_text_fit <- function(labels, slot_in) {
  need <- max(nchar(as.character(labels)), 1) * 0.065 + 0.1
  if (need <= slot_in) ggplot2::theme()
  else ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1, vjust = 1))
}

# Wrap text at spaces and after underscores (long taxon / group names)
wrap_text <- function(x, width) {
  vapply(as.character(x), function(s) {
    if (is.na(s) || nchar(s) <= width) return(s)
    lines <- vapply(strsplit(s, "\n", fixed = TRUE)[[1]], function(l)
      gsub("_ ", "_", paste(strwrap(gsub("_", "_ ", l), width = width), collapse = "\n")), "")
    paste(lines, collapse = "\n")
  }, "", USE.NAMES = FALSE)
}
# characters that fit in the narrowest facet strip when facets get space in proportion to their samples
strip_width_chars <- function(panels_in, n_per_group) {
  max(6, floor(panels_in * min(n_per_group) / sum(n_per_group) / 0.085))
}
# facet labeller that wraps strip labels to `width` characters
wrap_labeller <- function(width) ggplot2::as_labeller(function(x) wrap_text(x, width))

# PNG (300 dpi) + vector PDF; optional low-res copy *_mqc.png that MultiQC embeds as an image section.
# Title and subtitle are wrapped to the figure width so they are never cut off.
save_fig <- function(p, name, width = 8, height = 6, mqc = FALSE) {
  if (!is.null(p$labels$title))    p$labels$title    <- wrap_text(p$labels$title, floor((width - 0.4) / 0.115))
  if (!is.null(p$labels$subtitle)) p$labels$subtitle <- wrap_text(p$labels$subtitle, floor((width - 0.4) / 0.085))
  ggplot2::ggsave(paste0(name, ".png"), p, width = width, height = height, dpi = 300, bg = "white")
  ggplot2::ggsave(paste0(name, ".pdf"), p, width = width, height = height, bg = "white")
  if (mqc) ggplot2::ggsave(paste0(name, "_mqc.png"), p, width = width, height = height, dpi = 110, bg = "white")
  invisible(p)
}

# Table that MultiQC renders as custom content
write_mqc_table <- function(df, file, id, section, description = "", plot_type = "table") {
  writeLines(c(sprintf("# id: '%s'", id), sprintf("# section_name: '%s'", section),
               sprintf("# description: '%s'", gsub("'", "", description)),
               sprintf("# plot_type: '%s'", plot_type)), file)
  suppressWarnings(write.table(df, file, sep = "\t", quote = FALSE, row.names = FALSE, append = TRUE))
}

write_tsv <- function(df, file) write.table(df, file, sep = "\t", quote = FALSE, row.names = FALSE)

fmt_p <- function(p) ifelse(is.na(p), "NA", ifelse(p < 0.001, "< 0.001", sprintf("= %.3f", p)))

# ---------- taxonomy ----------
# Fill unassigned ranks with the last assigned name: NA genus in Lachnospiraceae -> "Lachnospiraceae (unclassified)"
fix_tax <- function(tt) {
  tt <- as.matrix(tt)
  for (i in seq_len(nrow(tt))) {
    last <- NA
    for (r in colnames(tt)) {
      v <- tt[i, r]
      if (is.na(v) || v == "") tt[i, r] <- if (is.na(last)) "Unclassified" else paste0(last, " (unclassified)")
      else last <- v
    }
  }
  tt
}

write_session <- function(step) {
  writeLines(capture.output(sessionInfo()), paste0(step, "_sessionInfo.txt"))
  if (file.exists("Rplots.pdf")) invisible(file.remove("Rplots.pdf"))   # stray default device
}
