# Image for PHYLOSEQ_BUILD, COMPOSITION, ALPHA/BETA_DIVERSITY and MAASLIN3.
# Published as ghcr.io/<owner>/auto16s-r by .github/workflows/r-image.yml, so users normally do not build it.
# Local build:  docker build -t auto16s-r:local containers/r-analysis
FROM bioconductor/bioconductor_docker:RELEASE_3_22

LABEL org.opencontainers.image.source="https://github.com/ozlemsagiroglu/auto16s" \
      org.opencontainers.image.description="R packages for the auto16s pipeline (phyloseq, vegan, MaAsLin3, ggplot2)" \
      org.opencontainers.image.licenses="MIT"

RUN apt-get update && apt-get install -y --no-install-recommends procps \
 && rm -rf /var/lib/apt/lists/*

COPY install_r_packages.R /tmp/install_r_packages.R
RUN Rscript /tmp/install_r_packages.R && rm /tmp/install_r_packages.R
