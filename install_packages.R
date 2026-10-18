# ==============================================================================
# Script Instalasi Paket R & Bioconductor
# Analisis RNA-seq Diferensial & Pengayaan Fungsional (GSE133053)
# R version 4.5.2 / Bioconductor 3.22
# ==============================================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

# Paket CRAN
cran_packages <- c(
    "ggplot2",       # v4.0.2
    "pheatmap",      # v1.0.13
    "ggrepel",       # v0.9.7
    "patchwork",     # v1.3.2
    "scales",
    "readr"
)

for (pkg in cran_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg, repos = "https://cloud.r-project.org")
    }
}

# Paket Bioconductor
bioc_packages <- c(
    "DESeq2",          # v1.50.2
    "apeglm",          # v1.32.0
    "clusterProfiler", # v4.18.4
    "org.At.tair.db",  # v3.22.0
    "AnnotationDbi",   # v1.72.0
    "KEGGREST"         # v1.50.0
)

for (pkg in bioc_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        BiocManager::install(pkg, ask = FALSE, update = FALSE)
    }
}

message("Instalasi seluruh dependensi R selesai!")
