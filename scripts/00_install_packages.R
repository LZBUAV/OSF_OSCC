# =============================================================================
# Install the R packages used by the analysis (versions used in the study: environment/R_packages.csv)
# Run from scripts/:  Rscript 00_install_packages.R      (safe to re-run)
# Python (pySCENIC) environment: see environment/python_pyscenic_requirements.txt and README.md
# =============================================================================
options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 3600)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
installed <- function() rownames(installed.packages())

cran_pkgs <- c("Seurat", "SeuratObject", "harmony", "Matrix", "data.table", "dplyr", "tidyr", "reshape2", "jsonlite", "png",
               "ggplot2", "ggpubr", "ggrepel", "patchwork", "pheatmap", "ggVennDiagram", "msigdbr", "igraph",
               "glmnet", "pROC", "randomForest", "xgboost", "gbm", "e1071", "plsRglm", "mboost", "naivebayes", "MASS", "shapviz")
bioc_pkgs <- c("limma", "GSVA", "UCell", "fgsea", "clusterProfiler", "ComplexHeatmap", "Biobase", "monocle")

miss <- setdiff(cran_pkgs, installed()); if (length(miss)) install.packages(miss)
miss <- setdiff(bioc_pkgs, installed()); if (length(miss)) BiocManager::install(miss, update = FALSE, ask = FALSE)
if (!"monocle" %in% installed()) remotes::install_github("cole-trapnell-lab/monocle-release", upgrade = "never")
if (!"CellChat" %in% installed()) remotes::install_github("jinworks/CellChat", upgrade = "never")
if (!"spacexr" %in% installed()) remotes::install_github("dmcable/spacexr", upgrade = "never")
if (!"CytoTRACE2" %in% installed())
  remotes::install_github("digitalcytometry/cytotrace2", subdir = "cytotrace2_r", upgrade = "never")

all_pkgs <- c(cran_pkgs, bioc_pkgs, "CellChat", "spacexr", "CytoTRACE2")
status <- setNames(all_pkgs %in% installed(), all_pkgs)
if (all(status)) cat("All packages installed.\n") else cat("Missing:", paste(names(status)[!status], collapse = ", "), "\n")
