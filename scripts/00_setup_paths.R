## ============================================================
## 00_setup_paths.R — Global setup (GitHub-safe)
## Purpose:
##   - Define project root (repo root)
##   - Define ALL input/output directories
##   - Create required folder structure
## ============================================================

## ---------------------------
## 1) PROJECT ROOT
## ---------------------------
# IMPORTANT:
# All scripts must be run from the repository root
# (gbm-multiomics-thesis)

ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

is_repo_root <- function(p) {
  dir.exists(file.path(p, "scripts")) &&
    dir.exists(file.path(p, "results")) &&
    dir.exists(file.path(p, "data_source"))
}

max_up <- 6
i <- 0
while (!is_repo_root(ROOT) && i < max_up) {
  ROOT <- normalizePath(dirname(ROOT), winslash = "/", mustWork = TRUE)
  i <- i + 1
}

stopifnot(is_repo_root(ROOT))

## 2) INPUT DATA
data_source_dir <- file.path(ROOT, "data_source")
RNA_INPUT_FILE   <- file.path(data_source_dir, "rna_input.rds")
METH_INPUT_FILE  <- file.path(data_source_dir, "meth_input.rds")
SAMPLE_META_FILE <- file.path(data_source_dir, "sample_meta.rds")

stopifnot(file.exists(RNA_INPUT_FILE),
          file.exists(METH_INPUT_FILE),
          file.exists(SAMPLE_META_FILE))

## 3) RESULTS (top-level)
results_dir <- file.path(ROOT, "results")
mofa_dir    <- file.path(results_dir, "MOFA")
gfa_dir     <- file.path(results_dir, "GFA")
compare_dir <- file.path(results_dir, "COMPARE")

dir.create(mofa_dir,    recursive = TRUE, showWarnings = FALSE)
dir.create(gfa_dir,     recursive = TRUE, showWarnings = FALSE)
dir.create(compare_dir, recursive = TRUE, showWarnings = FALSE)

## 4) STANDARD SUBFOLDERS (MOFA & GFA)
make_model_folders <- function(base_dir) {
  dirs <- c(
    file.path(base_dir, "rds"),
    file.path(base_dir, "GSEA", "plots"),
    file.path(base_dir, "GSEA", "tables"),
    file.path(base_dir, "ORA",  "plots"),
    file.path(base_dir, "ORA",  "tables")
  )
  for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  invisible(dirs)
}

make_model_folders(mofa_dir)
make_model_folders(gfa_dir)

## 5) QUICK PRINT (optional but useful)
message("ROOT: ", ROOT)
message("data_source_dir: ", data_source_dir)
message("results_dir: ", results_dir)