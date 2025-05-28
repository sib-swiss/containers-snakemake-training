# Use rocker image with specific R version
FROM rocker/r-ver:4.1.3

# Add tags
LABEL org.opencontainers.image.created="2025-05-28" \
    org.opencontainers.image.authors="Antonin Thiébaut" \
    org.opencontainers.image.description="Container with DESeq2 and friends"

# Change to bash shell
SHELL ["/bin/bash", "-c"]

# Update repos
RUN apt update
RUN apt-get update && apt install libz-dev libxml2 && apt clean

# Install R packages with pinned version using remotes package
RUN R --slave -e "install.packages('remotes', repos = 'https://cloud.r-project.org')"
RUN R --slave -e "remotes::install_version('data.table', version = '1.14.6', repos = 'https://cloud.r-project.org')"
RUN R --slave -e "remotes::install_version('ggplot2', version = '3.4.0', repos = 'https://cloud.r-project.org')"
RUN R --slave -e "remotes::install_version('ggrepel', version = '0.9.2', repos = 'https://cloud.r-project.org')"
RUN R --slave -e "remotes::install_version('pheatmap', version = '1.0.12', repos = 'https://cloud.r-project.org')"
RUN R --slave -e "remotes::install_version('RColorBrewer', version = '1.1.3', repos = 'https://cloud.r-project.org')"
RUN R --slave -e "remotes::install_version('scales', version = '1.3.0', repos = 'https://cloud.r-project.org')"
# RUN R --slave -e "remotes::install_version('scales', version = '1.2.1', repos = 'https://cloud.r-project.org')"
RUN R --slave -e "remotes::install_version('stringr', version = '1.5.0', repos = 'https://cloud.r-project.org')"

# Install Bioconductor and set its version
RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install(version = '3.14')"

# Install Bioconductor packages
RUN R -e "BiocManager::install('apeglm')"
RUN R -e "BiocManager::install('DESeq2')"
RUN R -e "BiocManager::install('IHW')"
RUN R -e "BiocManager::install('limma')"

# Set working directory
WORKDIR /opt
