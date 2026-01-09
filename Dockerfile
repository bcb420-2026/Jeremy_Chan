FROM risserlin/bcb420-base-image:winter2026-arm64 AS final

RUN R -q -e "install.packages('BiocManager'); \
             BiocManager::install(c('DESeq2', 'enrichplot'), ask = FALSE, update = FALSE); \
             install.packages('pheatmap')"