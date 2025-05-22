FROM rocker/r2u:jammy

LABEL org.opencontainers.image.created="2023-04-12" \
    org.opencontainers.image.authors="Geert van Geest" \
    org.opencontainers.image.description="Container with DESeq2 and friends"

RUN install2.r \
    apeglm \
    BiocParallel \
    data.table \
    DESeq2 \
    ggplot2 \
    ggrepel \
    IHW \
    limma \
    optparse \
    pheatmap \
    RColorBrewer \
    scales \
    stringr

WORKDIR /opt

COPY test_deseq2.R .

ENV PATH=/opt:$PATH
