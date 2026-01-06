## =========================================================
## Script: MOFA_SETUP.R
## Purpose:
##   - Install/check MOFA2
##   - Prefer basilisk backend (recommended)
##   - Optional: reticulate/miniconda fallback check
## Author: Havva Ünal
## =========================================================

## ---------------------------
## 0) R / Bioconductor setup
## ---------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

## ---------------------------
## 1) Install MOFA2 (if needed)
## ---------------------------
if (!requireNamespace("MOFA2", quietly = TRUE)) {
  BiocManager::install("MOFA2")
}

suppressPackageStartupMessages({
  library(MOFA2)
})

message("MOFA2 version: ", as.character(packageVersion("MOFA2")))

## ---------------------------
## 2) Backend choice (recommended)
## ---------------------------
## We will run MOFA using basilisk=TRUE in the model script.
## basilisk manages a private Python env automatically.
backend_mode <- "basilisk"
message("Backend mode set to: ", backend_mode)

## ---------------------------
## 3) Quick sanity check (no model fitting here)
## ---------------------------
## Just verify the package loads and core functions exist.
stopifnot(exists("create_mofa", where = asNamespace("MOFA2")))
stopifnot(exists("run_mofa",    where = asNamespace("MOFA2")))

message("MOFA2 setup OK")
