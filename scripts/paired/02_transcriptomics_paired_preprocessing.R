# ==============================================================
# Script: 02_transcriptomics_paired_tumor_prepare.R
# Purpose:
#   Paired transcriptomics (tumour-only) preprocessing for GitHub:
#     - Read existing tumour-only RNA RDS objects (no downloading)
#     - Prefilter genes
#     - DESeq2 size-factor normalization (design ~ 1)
#     - Save log-normalized counts + VST matrices
#     - Clean ENSEMBL IDs (remove version suffix)
#     - Map ENSEMBL -> SYMBOL/GENENAME and save mapping table
#     - Save tumour-only outputs for downstream paired integration
#
# Inputs (repo):
#   data/transcriptomics/paired/rna_counts_tumor_raw.rds
#   data/transcriptomics/paired/rna_metadata_tumor.rds
#
# Outputs (repo):
#   results/transcriptomics_paired/tumor/*
#
# Author: Havva Ünal
# ==============================================================

source("scripts/00_setup_paths.R")

## =============================================================
## 0) PATHS (GitHub-stable; no mode dependence)
## =============================================================
# Expect these from setup
if (!exists("ROOT"))        stop("ROOT not found. Check scripts/00_setup_paths.R")
if (!exists("DATA_DIR"))    stop("DATA_DIR not found. Check scripts/00_setup_paths.R")
if (!exists("RESULTS_DIR")) stop("RESULTS_DIR not found. Check scripts/00_setup_paths.R")

TX_DIR  <- file.path(DATA_DIR, "transcriptomics", "paired")
OUT_DIR <- file.path(RESULTS_DIR, "transcriptomics_paired")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

message("TX input dir : ", TX_DIR)
message("TX output dir: ", OUT_DIR)

## =============================================================
## 1) LIBRARIES
## =============================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(DESeq2)
  library(matrixStats)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

## =============================================================
## 2) LOAD TUMOUR-ONLY RNA OBJECTS
## =============================================================
rna_counts_path <- file.path(TX_DIR, "rna_counts_tumor_raw.rds")
rna_meta_path   <- file.path(TX_DIR, "rna_metadata_tumor.rds")

if (!file.exists(rna_counts_path)) stop("Missing RNA counts file: ", rna_counts_path)
if (!file.exists(rna_meta_path))   stop("Missing RNA metadata file: ", rna_meta_path)

rna_counts <- readRDS(rna_counts_path)
rna_meta   <- readRDS(rna_meta_path)

if (ncol(rna_counts) != nrow(rna_meta)) {
  stop("Counts columns and metadata rows do not match: ",
       ncol(rna_counts), " vs ", nrow(rna_meta))
}

if (!("barcode" %in% colnames(rna_meta))) {
  stop("RNA metadata must contain a 'barcode' column.")
}

# Align metadata to counts columns (strict)
if (!identical(colnames(rna_counts), rna_meta$barcode)) {
  idx <- match(colnames(rna_counts), rna_meta$barcode)
  if (any(is.na(idx))) stop("Some count sample barcodes are missing in metadata.")
  rna_meta <- rna_meta[idx, , drop = FALSE]
  stopifnot(identical(colnames(rna_counts), rna_meta$barcode))
}

# Patient ID
if (!("patient_id" %in% colnames(rna_meta))) {
  rna_meta$patient_id <- substr(rna_meta$barcode, 1, 12)
}

message("Loaded tumour RNA counts: ", paste(dim(rna_counts), collapse = " x "))
message("Unique patients: ", length(unique(rna_meta$patient_id)))

# Save metadata copy into OUT_DIR
saveRDS(rna_meta, file = file.path(OUT_DIR, "rna_metadata_tumor.rds"))

## =============================================================
## 3) QC: LIBRARY SIZES
## =============================================================
lib_sizes <- colSums(rna_counts)
lib_df <- data.frame(
  barcode      = colnames(rna_counts),
  patient_id   = rna_meta$patient_id,
  library_size = as.numeric(lib_sizes)
)
write.csv(lib_df, file.path(OUT_DIR, "sample_library_sizes.csv"), row.names = FALSE)

## =============================================================
## 4) GENE PREFILTERING
##    Keep genes with >= 10 counts in >= 3 samples
## =============================================================
keep_genes <- rowSums(rna_counts >= 10) >= 3
rna_counts_filt <- rna_counts[keep_genes, , drop = FALSE]

message("Prefilter kept genes: ", sum(keep_genes), " / ", nrow(rna_counts))

saveRDS(rna_counts_filt, file.path(OUT_DIR, "rna_counts_tumor_prefiltered.rds"))

qc_feat <- data.frame(
  metric = c("n_genes_raw", "n_genes_prefiltered", "n_samples"),
  value  = c(nrow(rna_counts), nrow(rna_counts_filt), ncol(rna_counts))
)
write.csv(qc_feat, file.path(OUT_DIR, "rna_feature_qc_summary.csv"), row.names = FALSE)

## =============================================================
## 5) NORMALIZATION (DESeq2 SIZE FACTORS; design ~ 1)
## =============================================================
dds <- DESeqDataSetFromMatrix(
  countData = round(rna_counts_filt),
  colData   = rna_meta,
  design    = ~ 1
)

dds <- estimateSizeFactors(dds)

# log2 normalized counts
rna_log <- log2(counts(dds, normalized = TRUE) + 1)
stopifnot(identical(colnames(rna_log), rna_meta$barcode))
saveRDS(rna_log, file.path(OUT_DIR, "rna_logcounts_tumor.rds"))

# VST (for unsupervised integration)
rna_vst <- assay(vst(dds, blind = TRUE))
stopifnot(identical(colnames(rna_vst), rna_meta$barcode))
saveRDS(rna_vst, file.path(OUT_DIR, "rna_vst_tumor.rds"))

## =============================================================
## 6) ENSEMBL CLEANUP + SYMBOL MAPPING
## =============================================================
ens_raw   <- rownames(rna_log)
ens_clean <- sub("\\..*$", "", ens_raw)

rownames(rna_log) <- ens_clean
rownames(rna_vst) <- ens_clean

map_df <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = unique(ens_clean),
  keytype = "ENSEMBL",
  columns = c("SYMBOL", "GENENAME")
)

map_df <- map_df %>%
  dplyr::filter(!is.na(ENSEMBL)) %>%
  dplyr::arrange(ENSEMBL, SYMBOL, GENENAME) %>%
  dplyr::distinct(ENSEMBL, .keep_all = TRUE)

write.csv(map_df, file.path(OUT_DIR, "gene_id_mapping_ensembl_symbol.csv"), row.names = FALSE)

# SYMBOL matrix versions (drop duplicated SYMBOLs deterministically)
symbol_vec <- map_df$SYMBOL[match(ens_clean, map_df$ENSEMBL)]
sym_names  <- ifelse(is.na(symbol_vec) | symbol_vec == "", ens_clean, symbol_vec)

rna_log_sym <- rna_log
rownames(rna_log_sym) <- sym_names
rna_log_sym <- rna_log_sym[!duplicated(rownames(rna_log_sym)), , drop = FALSE]
saveRDS(rna_log_sym, file.path(OUT_DIR, "rna_logcounts_tumor_SYMBOL.rds"))

rna_vst_sym <- rna_vst
rownames(rna_vst_sym) <- sym_names
rna_vst_sym <- rna_vst_sym[!duplicated(rownames(rna_vst_sym)), , drop = FALSE]
saveRDS(rna_vst_sym, file.path(OUT_DIR, "rna_vst_tumor_SYMBOL.rds"))

## =============================================================
## 7) TOP VARIABLE GENES (VST, SYMBOL)
## =============================================================
gene_var <- matrixStats::rowVars(rna_vst_sym)
names(gene_var) <- rownames(rna_vst_sym)

topN <- 2000
top_genes <- names(sort(gene_var, decreasing = TRUE))[seq_len(min(topN, length(gene_var)))]
writeLines(top_genes, con = file.path(OUT_DIR, "top_variable_genes_2000.txt"))

message("DONE. Paired tumour transcriptomics outputs written to: ", OUT_DIR)


