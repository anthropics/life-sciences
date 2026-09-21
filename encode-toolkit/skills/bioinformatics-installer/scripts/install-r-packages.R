#!/usr/bin/env Rscript
# Install R/Bioconductor packages for ENCODE data analysis
# Usage: Rscript install-r-packages.R [--all | --core | --chipseq | --rnaseq | --singlecell | --methylation | --deconvolution | --visualization | --stats]
#
# Categories:
#   --all            Install all packages (default if no argument)
#   --core           GenomicRanges, rtracklayer, biomaRt, etc.
#   --chipseq        DiffBind, ChIPQC, ChIPseeker, chromVAR
#   --rnaseq         DESeq2, edgeR, limma
#   --singlecell     Seurat, Signac, scater, scran
#   --methylation    DMRcate, bsseq, methylKit
#   --deconvolution  BisqueRNA, DWLS (prints the GitHub install lines for BayesPrism and InstaPrism)
#   --visualization  ComplexHeatmap, EnhancedVolcano, Gviz, ggplot2
#   --stats          sva, WGCNA, ReactomePA

# --- Setup ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) args <- "--all"

cat("============================================\n")
cat("ENCODE Bioinformatics R Package Installer\n")
cat("============================================\n\n")

# Install BiocManager if not present
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  cat("Installing BiocManager...\n")
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

# Lock Bioconductor version based on R version
r_version <- paste(R.version$major, R.version$minor, sep = ".")
cat(sprintf("R version: %s\n", r_version))

bioc_version <- BiocManager::version()
cat(sprintf("Bioconductor version: %s\n\n", bioc_version))

# Packages that could not be installed or loaded; the script exits non-zero when any are left
failed_pkgs <- character(0)

# Helper function
install_pkgs <- function(pkgs, category) {
  cat(sprintf("\n--- Installing %s packages (%d) ---\n", category, length(pkgs)))
  for (pkg in pkgs) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      cat(sprintf("  [OK] %s already installed\n", pkg))
    } else {
      cat(sprintf("  [..] Installing %s...\n", pkg))
      tryCatch({
        BiocManager::install(pkg, ask = FALSE, update = FALSE)
        if (requireNamespace(pkg, quietly = TRUE)) {
          cat(sprintf("  [OK] %s installed successfully\n", pkg))
        } else {
          cat(sprintf("  [!!] %s install completed but package not loadable\n", pkg))
          failed_pkgs <<- c(failed_pkgs, pkg)
        }
      }, error = function(e) {
        cat(sprintf("  [FAIL] %s: %s\n", pkg, conditionMessage(e)))
        failed_pkgs <<- c(failed_pkgs, pkg)
      })
    }
  }
}

# --- Package Definitions ---

core_pkgs <- c(
  "GenomicRanges", "GenomicFeatures", "rtracklayer", "IRanges",
  "GenomeInfoDb", "BiocGenerics", "S4Vectors", "AnnotationDbi", "biomaRt"
)

chipseq_pkgs <- c(
  "DiffBind", "ChIPQC", "ChIPseeker", "chromVAR",
  "annotatr", "clusterProfiler",
  "org.Hs.eg.db", "org.Mm.eg.db",
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "TxDb.Mmusculus.UCSC.mm10.knownGene",
  "TFBSTools", "motifmatchr"
)

rnaseq_pkgs <- c(
  "DESeq2", "edgeR", "limma",
  "tximport", "tximeta"
)

singlecell_pkgs <- c(
  "Seurat", "Signac",
  "SingleCellExperiment", "scater", "scran",
  "celda"
)

deconvolution_pkgs <- c(
  "BisqueRNA", "DWLS"
  # BayesPrism and InstaPrism require GitHub install:
  # devtools::install_github("Danko-Lab/BayesPrism/BayesPrism")
  # devtools::install_github("humengying0907/InstaPrism")
)

methylation_pkgs <- c(
  "DMRcate", "bsseq", "methylKit"
)

visualization_pkgs <- c(
  "ComplexHeatmap", "EnhancedVolcano", "Gviz", "ggplot2",
  "pheatmap", "RColorBrewer", "viridis"
)

stats_pkgs <- c(
  "sva", "WGCNA", "ReactomePA"
)

# --- Execute Based on Arguments ---

install_all <- "--all" %in% args

if (install_all || "--core" %in% args) {
  install_pkgs(core_pkgs, "Core Genomic Infrastructure")
}

if (install_all || "--chipseq" %in% args) {
  install_pkgs(chipseq_pkgs, "ChIP-seq / Peak Analysis")
}

if (install_all || "--rnaseq" %in% args) {
  install_pkgs(rnaseq_pkgs, "RNA-seq Differential Expression")
}

if (install_all || "--singlecell" %in% args) {
  install_pkgs(singlecell_pkgs, "Single-Cell Analysis")
}

if (install_all || "--deconvolution" %in% args) {
  install_pkgs(deconvolution_pkgs, "Deconvolution (CRAN/Bioc)")
  cat("\n--- GitHub-only deconvolution packages ---\n")
  cat("  Run manually:\n")
  cat("    devtools::install_github('Danko-Lab/BayesPrism/BayesPrism')\n")
  cat("    devtools::install_github('humengying0907/InstaPrism')\n")
  cat("    devtools::install_github('dtsoucas/DWLS')   # if CRAN version unavailable\n")
}

if (install_all || "--methylation" %in% args) {
  install_pkgs(methylation_pkgs, "DNA Methylation Analysis")
}

if (install_all || "--visualization" %in% args) {
  install_pkgs(visualization_pkgs, "Visualization")
}

if (install_all || "--stats" %in% args) {
  install_pkgs(stats_pkgs, "Statistics and Batch Correction")
}

# --- Summary ---
cat("\n============================================\n")
if (length(failed_pkgs) > 0) {
  cat(sprintf("Installation INCOMPLETE: %d package(s) failed:\n  %s\n",
              length(failed_pkgs), paste(failed_pkgs, collapse = ", ")))
  cat("Fix the errors above and run the script again.\n")
  cat("============================================\n")
  quit(status = 1)
}
cat("Installation complete.\n")
cat("Run BiocManager::valid() to check for version mismatches.\n")
cat("============================================\n")
