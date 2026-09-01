# Packages required by the notebooks in analysis_code/ and models_code/.
#
# Run once from the repository root:
#
#     Rscript install_packages.R
#
# Only packages that are not already available are installed, so it is safe to
# re-run. The list covers every package the notebooks load directly;
# dependencies are resolved by the installers.

cran <- c(
  "broom",
  "caret",
  "circlize",
  "colorspace",
  "ComplexUpset",
  "compositions",
  "cowplot",
  "data.table",
  "devtools",
  "dplyr",
  "FSelectorRcpp",
  "ggnewscale",
  "ggplot2",
  "ggpubr",
  "ggrepel",
  "ggsci",
  "ggsignif",
  "glmnet",       # needed to work with the saved caret/glmnet classifier
  "gridExtra",
  "here",
  "HGNChelper",
  "lemon",
  "magrittr",
  "patchwork",
  "pheatmap",
  "plyr",
  "ppcor",
  "pROC",
  "RColorBrewer",
  "readr",
  "readxl",
  "recipes",
  "reshape2",
  "reticulate",
  "scales",
  "stringr",
  "tidyr",
  "tidyverse",
  "zoo"
)

bioc <- c(
  "biomaRt",
  "Biostrings",
  "BSgenome.Hsapiens.UCSC.hg19",
  "clusterProfiler",
  "ComplexHeatmap",
  "dorothea",
  "enrichplot",
  "GenomeInfoDb",
  "GenomicRanges",
  "IRanges",
  "org.Hs.eg.db",
  "preprocessCore",
  "rtracklayer",
  "SummarizedExperiment"
)

not_installed <- function(pkgs) {
  pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
}

for (pkg in c("BiocManager", "remotes")) {
  if (length(not_installed(pkg))) install.packages(pkg)
}

if (length(not_installed(cran))) install.packages(not_installed(cran))
if (length(not_installed(bioc))) BiocManager::install(not_installed(bioc))

# rlucas is not on CRAN or Bioconductor; it ships inside the LUCAS workflow repo.
# A copy is also vendored at data/reproduce_lucas_wflow/code/rlucas.
if (length(not_installed("rlucas"))) {
  remotes::install_github(
    "cancer-genomics/reproduce_lucas_wflow",
    subdir = "code/rlucas"
  )
}

still_missing <- not_installed(c(cran, bioc, "rlucas"))
if (length(still_missing)) {
  warning("Not installed: ", paste(still_missing, collapse = ", "))
} else {
  message("All required packages are available.")
}
