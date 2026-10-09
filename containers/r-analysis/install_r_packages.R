# R packages for the phyloseq / vegan / MaAsLin3 / ggplot2 steps.
# Used by the Dockerfile and the Apptainer definition; can also be run on a plain R installation (R >= 4.4, Bioconductor >= 3.20):
#   Rscript containers/r-analysis/install_r_packages.R
repos <- "https://cloud.r-project.org"
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = repos)
install.packages(c("vegan", "ggplot2", "scales"), repos = repos)
BiocManager::install(c("phyloseq", "maaslin3"), ask = FALSE, update = FALSE)
for (p in c("phyloseq", "maaslin3", "vegan", "ggplot2", "scales")) {
  if (!requireNamespace(p, quietly = TRUE)) stop("Package not installed: ", p)
  message(p, " ", as.character(packageVersion(p)))
}
