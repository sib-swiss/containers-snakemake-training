# Use rocker image with R v4.1.3 as base image
FROM rocker/r-ver:4.1.3

# Add metadata to image
LABEL org.opencontainers.image.created="2026-06-03" \
    org.opencontainers.image.authors="Antonin Thiébaut" \
    org.opencontainers.image.description="Container with DESeq2 and friends"

# Install system dependencies
RUN apt update && apt install -y \
    fonts-noto \
    git \
    libcurl4-openssl-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libfribidi-dev \
    libgit2-dev \
    libharfbuzz-dev \
    libjpeg-dev \
    libpng-dev \
    libssl-dev \
    libtiff5-dev \
    libxml2-dev \
    tzdata \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Install renv
RUN R -e "install.packages('renv', repos = c(CRAN = 'https://cloud.r-project.org'))"

# Set working directory
WORKDIR /project

# Copy renv file inside container
COPY renv.lock renv.lock

# Set renv and R library paths
ENV RENV_PATHS_LIBRARY=renv/library
ENV R_LIBS_USER=/project/renv/library

# Change default shell to Rscript
SHELL ["Rscript", "-e"]

# Restore all packages
RUN dir.create(Sys.getenv('RENV_PATHS_LIBRARY'), recursive = TRUE, showWarnings = FALSE); \
    .libPaths(c(Sys.getenv('RENV_PATHS_LIBRARY'), .libPaths())); \
    options(repos = c( \
      CRAN     = 'https://cloud.r-project.org', \
      BioCsoft = 'https://bioconductor.org/packages/3.14/bioc', \
      BioCann  = 'https://bioconductor.org/packages/3.14/data/annotation', \
      BioCexp  = 'https://bioconductor.org/packages/3.14/data/experiment', \
      BioCother= 'https://bioconductor.org/packages/3.14/extra' \
    )); \
    options(renv.config.repos.override = getOption('repos')); \
    install.packages('BiocManager'); \
    BiocManager::install(version = '3.14', ask = FALSE, update = FALSE); \
    renv::restore(lockfile = 'renv.lock', prompt = FALSE)
