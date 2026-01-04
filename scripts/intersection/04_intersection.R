## ==============================================================
## 04_intersection.R — Paired RNA/METH intersection (PROJECT-SAFE)
## --------------------------------------------------------------
## IMPORTANT NOTE (DO NOT REMOVE):
## - This script MUST use the central path config:
##     source("scripts/00_setup_paths.R")
## - Do NOT hard-code paths elsewhere.
## - All outputs of this script MUST be written ONLY to:
##     results/intersection/
## - Inputs are read from:
##     results/transcriptomics_paired/   (RNA tumour SYMBOL matrix)
##     results/paired_inputs/            (paired metadata + METH paired promoter matrix)
## - This script creates the FINAL paired matrices used by MOFA/GFA:
##     rna_paired_final.rds
##     meth_paired_final.rds
##     paired_metadata_final.rds
##     paired_qc_summary.txt
## ==============================================================
## Steps:
##   STEP 0) Setup paths + output folder (results/intersection)
##   STEP 1) Create RNA paired matrix in SYMBOL space (sample subset + order)
##   STEP 2) Intersect genes (RNA vs METH) and save FINAL paired inputs
## ==============================================================
##
## Author: Havva Ünal
## ==============================================================


## =============================================================
## STEP 0) SETUP PATHS (DEFINE ONCE)
## =============================================================
source("scripts/00_setup_paths.R")

# Canonicalize names if setup uses aliases
if (!exists("RESULTS_DIR") && exists("results_dir")) RESULTS_DIR <- results_dir
if (!exists("RESULTS_DIR")) stop("RESULTS_DIR/results_dir not found. Check scripts/00_setup_paths.R")

INTERSECT_DIR <- file.path(RESULTS_DIR, "intersection")
dir.create(INTERSECT_DIR, recursive = TRUE, showWarnings = FALSE)

# Fixed input roots (repo-safe)
TX_PAIRED_DIR <- file.path(RESULTS_DIR, "transcriptomics_paired")
PAIRED_DIR    <- file.path(RESULTS_DIR, "paired_inputs")

if (!dir.exists(TX_PAIRED_DIR)) stop("Missing folder: ", TX_PAIRED_DIR)
if (!dir.exists(PAIRED_DIR))    stop("Missing folder: ", PAIRED_DIR)

message("TX_PAIRED_DIR : ", TX_PAIRED_DIR)
message("PAIRED_DIR    : ", PAIRED_DIR)
message("INTERSECT_DIR : ", INTERSECT_DIR)


## =============================================================
## Helper: pick first existing file from candidate paths
## =============================================================
first_existing <- function(paths) {
  p <- paths[file.exists(paths)][1]
  if (length(p) == 0 || is.na(p)) return(NA_character_)
  p
}


## =============================================================
## STEP 1) CREATE RNA PAIRED MATRIX IN SYMBOL SPACE
## Input:
##   results/transcriptomics_paired/rna_logcounts_tumor_SYMBOL.rds
##   results/paired_inputs/rna_metadata_paired.rds
## Output:
##   results/intersection/rna_logcounts_paired_SYMBOL.rds
## =============================================================
RNA_SYM_TUMOR <- file.path(TX_PAIRED_DIR, "rna_logcounts_tumor_SYMBOL.rds")
RNA_META_PAIRED <- file.path(PAIRED_DIR, "rna_metadata_paired.rds")
OUT_RNA_PAIRED_SYM <- file.path(INTERSECT_DIR, "rna_logcounts_paired_SYMBOL.rds")

if (!file.exists(RNA_SYM_TUMOR)) stop("Missing file: ", RNA_SYM_TUMOR)
if (!file.exists(RNA_META_PAIRED)) stop("Missing file: ", RNA_META_PAIRED)

rna_sym_tumor  <- readRDS(RNA_SYM_TUMOR)                 # genes x tumour samples (SYMBOL)
rna_meta_paired <- as.data.frame(readRDS(RNA_META_PAIRED))

if (!("barcode" %in% colnames(rna_meta_paired))) {
  stop("rna_metadata_paired.rds must contain a column named 'barcode'.")
}
rna_meta_paired$barcode <- as.character(rna_meta_paired$barcode)

barcodes <- rna_meta_paired$barcode

missing <- setdiff(barcodes, colnames(rna_sym_tumor))
if (length(missing) > 0) {
  stop("Some paired barcodes are missing from RNA SYMBOL tumour matrix. Example: ", missing[1])
}

rna_sym_paired <- rna_sym_tumor[, barcodes, drop = FALSE]
stopifnot(identical(colnames(rna_sym_paired), barcodes))

saveRDS(rna_sym_paired, OUT_RNA_PAIRED_SYM)

message("STEP 1 DONE")
message("Saved: ", OUT_RNA_PAIRED_SYM)
message("Dim  : ", paste(dim(rna_sym_paired), collapse = " x "))


## =============================================================
## STEP 2) INTERSECT GENES (RNA vs METH) + SAVE FINAL PAIRED INPUTS
## Inputs:
##   results/intersection/rna_logcounts_paired_SYMBOL.rds
##   results/paired_inputs/meth_promoter_paired_raw.rds  (or meth_promoter_paired.rds)
##   results/paired_inputs/rna_metadata_paired.rds
##   results/paired_inputs/meth_metadata_paired.rds
## Outputs (ONLY to results/intersection/):
##   rna_paired_final.rds
##   meth_paired_final.rds
##   paired_metadata_final.rds
##   paired_qc_summary.txt
## =============================================================
suppressPackageStartupMessages({
  library(dplyr)
})

rna_mat_path <- file.path(INTERSECT_DIR, "rna_logcounts_paired_SYMBOL.rds")

meth_mat_path <- first_existing(c(
  file.path(PAIRED_DIR, "meth_promoter_paired_raw.rds"),
  file.path(PAIRED_DIR, "meth_promoter_paired.rds")
))

meth_meta_path <- file.path(PAIRED_DIR, "meth_metadata_paired.rds")

if (!file.exists(rna_mat_path)) stop("Missing RNA paired SYMBOL matrix: ", rna_mat_path)
if (is.na(meth_mat_path)) stop("Missing METH paired promoter matrix in: ", PAIRED_DIR)
if (!file.exists(RNA_META_PAIRED)) stop("Missing: ", RNA_META_PAIRED)
if (!file.exists(meth_meta_path)) stop("Missing: ", meth_meta_path)

message("RNA  input: ", rna_mat_path)
message("METH input: ", meth_mat_path)

rna_mat  <- readRDS(rna_mat_path)     # genes x samples (SYMBOL)
meth_mat <- readRDS(meth_mat_path)    # genes x samples (should match SYMBOL gene space)

rna_meta  <- as.data.frame(readRDS(RNA_META_PAIRED))
meth_meta <- as.data.frame(readRDS(meth_meta_path))

# --- Align samples (RNA side)
if (!("barcode" %in% colnames(rna_meta))) stop("rna_metadata_paired.rds must contain 'barcode'")
rna_meta$barcode <- as.character(rna_meta$barcode)

if (!identical(colnames(rna_mat), rna_meta$barcode)) {
  rna_meta <- rna_meta[match(colnames(rna_mat), rna_meta$barcode), , drop = FALSE]
  stopifnot(identical(colnames(rna_mat), rna_meta$barcode))
}

# --- Align samples (METH side): meth_meta rownames should match meth_mat colnames
if (!identical(colnames(meth_mat), rownames(meth_meta))) {
  samp_col <- c("Sample ID", "sample_id", "barcode")[c("Sample ID","sample_id","barcode") %in% colnames(meth_meta)][1]
  if (!is.na(samp_col)) {
    meth_meta[[samp_col]] <- as.character(meth_meta[[samp_col]])
    meth_meta <- meth_meta[match(colnames(meth_mat), meth_meta[[samp_col]]), , drop = FALSE]
    rownames(meth_meta) <- colnames(meth_mat)
  }
  stopifnot(identical(colnames(meth_mat), rownames(meth_meta)))
}

# --- Patient ID consistency
if (!("patient_id" %in% colnames(rna_meta)))  rna_meta$patient_id  <- substr(rna_meta$barcode, 1, 12)
if (!("patient_id" %in% colnames(meth_meta))) meth_meta$patient_id <- substr(rownames(meth_meta), 1, 12)

stopifnot(
  ncol(rna_mat) == ncol(meth_mat),
  identical(as.character(rna_meta$patient_id), as.character(meth_meta$patient_id))
)

message("Samples aligned. Paired n = ", ncol(rna_mat))

# --- Gene intersection
common_genes <- intersect(rownames(rna_mat), rownames(meth_mat))
common_genes <- sort(unique(common_genes))

if (length(common_genes) < 1000) {
  stop("Too few common genes (", length(common_genes), "). Check gene IDs (both must be SYMBOL).")
}

rna_final  <- rna_mat[common_genes, , drop = FALSE]
meth_final <- meth_mat[common_genes, , drop = FALSE]
stopifnot(identical(rownames(rna_final), rownames(meth_final)))

message("Common genes (n): ", length(common_genes))
message("Final RNA dim : ", paste(dim(rna_final), collapse = " x "))
message("Final METH dim: ", paste(dim(meth_final), collapse = " x "))

# --- Save FINAL outputs (ONLY to INTERSECT_DIR)
saveRDS(rna_final,  file = file.path(INTERSECT_DIR, "rna_paired_final.rds"))
saveRDS(meth_final, file = file.path(INTERSECT_DIR, "meth_paired_final.rds"))

meta_final <- list(
  rna_meta     = rna_meta,
  meth_meta    = meth_meta,
  common_genes = common_genes,
  rna_source   = basename(rna_mat_path),
  meth_source  = basename(meth_mat_path)
)
saveRDS(meta_final, file = file.path(INTERSECT_DIR, "paired_metadata_final.rds"))

qc_lines <- c(
  paste0("paired_samples\t", ncol(rna_final)),
  paste0("common_genes\t", length(common_genes)),
  paste0("rna_matrix\t", basename(rna_mat_path)),
  paste0("meth_matrix\t", basename(meth_mat_path))
)
writeLines(qc_lines, con = file.path(INTERSECT_DIR, "paired_qc_summary.txt"))

message("STEP 2 DONE")
message("Saved FINAL paired inputs in: ", INTERSECT_DIR)
message("DONE.")
