# ==============================================================
# Script: 03_methylation_processing.R
# Purpose:
#   Build methylation inputs for TCGA-GBM (Illumina 450k):
#     (A) CpG-level beta + M-value matrices + phenotype table
#     (B) Promoter-level gene matrix (TSS200/TSS1500) for TUMOUR samples only
#
# Inputs (repo):
#   data/methylation/
#     - GDC_methylation/                         (GDC client download folder)
#     - gdc_sample_sheet.2025-09-22.tsv          (tumour sample sheet)
#     - gdc_sample_sheet_normals.tsv             (optional)
#
# Outputs (repo):
#   results/methylation/cpg/
#     - methyl_beta_matrix.rds
#     - methyl_M_matrix.rds
#     - methyl_pheno.rds
#   results/methylation/promoter/
#     - meth_gene_matrix_promoter_tumor.rds
#     - meth_promoter_pheno_tumor.rds
#
# Author: Havva Ünal
# ==============================================================

source("scripts/00_setup_paths.R")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(matrixStats)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
})

# -----------------------------
# 0) Paths
# -----------------------------
METH_DIR <- file.path(data_dir, "methylation")
GDC_DIR  <- file.path(METH_DIR, "GDC_methylation")

SS_TUMOR_TSV  <- file.path(METH_DIR, "gdc_sample_sheet.2025-09-22.tsv")
SS_NORMAL_TSV <- file.path(METH_DIR, "gdc_sample_sheet_normals.tsv")  

OUT_CPG  <- file.path(results_dir, "methylation", "cpg")
OUT_PROM <- file.path(results_dir, "methylation", "promoter")
dir.create(OUT_CPG,  recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_PROM, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 1) Parameters / helpers
# -----------------------------
ROW_COVERAGE <- 0.95

short_barcode <- function(x) substr(as.character(x), 1, 15)

# -----------------------------
# 2) Sanity checks
# -----------------------------
if (!dir.exists(GDC_DIR)) stop("Missing folder: ", GDC_DIR)
if (!file.exists(SS_TUMOR_TSV)) stop("Missing tumour sample sheet: ", SS_TUMOR_TSV)

beta_files <- list.files(GDC_DIR, pattern = "level3betas\\.txt$", recursive = TRUE, full.names = TRUE)
message("Found level3betas files: ", length(beta_files))
if (length(beta_files) == 0) {
  stop("No level3betas.txt found under: ", GDC_DIR,
       "\nCheck your GDC download path / manifest.")
}

# -----------------------------
# 3) Load sample sheets (tumour + optional normals)
# -----------------------------
ss_t <- readr::read_tsv(SS_TUMOR_TSV, show_col_types = FALSE)
if ("File ID" %in% colnames(ss_t) && !"file_id" %in% colnames(ss_t)) {
  ss_t <- dplyr::rename(ss_t, file_id = `File ID`)
}
if (!"file_id" %in% colnames(ss_t)) stop("Tumour sample sheet must contain 'File ID' or 'file_id'.")
if (!"Sample ID" %in% colnames(ss_t)) stop("Tumour sample sheet must contain 'Sample ID'.")

meta_all <- ss_t

if (file.exists(SS_NORMAL_TSV)) {
  ss_n <- readr::read_tsv(SS_NORMAL_TSV, show_col_types = FALSE)
  if ("File ID" %in% colnames(ss_n) && !"file_id" %in% colnames(ss_n)) {
    ss_n <- dplyr::rename(ss_n, file_id = `File ID`)
  }
  meta_all <- dplyr::bind_rows(meta_all, ss_n)
}

keep_cols <- c("file_id", "File Name", "Sample ID", "Case ID",
               "Tissue Type", "Tumor Descriptor", "Project ID")
meta_all <- meta_all[, intersect(keep_cols, colnames(meta_all)), drop = FALSE]

message("Sample sheet rows (total): ", nrow(meta_all))
message("Unique file_id in sheet   : ", length(unique(meta_all$file_id)))
message("Unique Sample ID in sheet : ", length(unique(meta_all$`Sample ID`)))

# -----------------------------
# 4) Match level3betas files to Sample IDs
#    (parent folder name = file_id)
# -----------------------------
all_files <- beta_files
fid <- basename(dirname(all_files))

idx <- match(fid, meta_all$file_id)
sids <- meta_all$`Sample ID`[idx]

keep <- !is.na(sids)
all_files <- all_files[keep]
sids <- sids[keep]

# Standardise sample ids to 15-char TCGA sample barcode
sids <- short_barcode(sids)

message("Matched methylation files: ", length(all_files))
message("Unique samples (15-char): ", length(unique(sids)))

if (length(all_files) < 2) stop("Too few matched methylation files. Check file_id matching.")

# -----------------------------
# 5) Build CpG beta matrix (CpG x sample)
# -----------------------------
ref <- read.delim(all_files[1], header = FALSE, stringsAsFactors = FALSE)[, 1]

beta_all <- do.call(cbind, lapply(all_files, function(f) {
  df <- read.delim(f, header = FALSE, stringsAsFactors = FALSE)
  v  <- setNames(df[, 2], df[, 1])
  v[ref]
}))

rownames(beta_all) <- ref
colnames(beta_all) <- sids

# Remove columns with all NA
beta_ok <- beta_all[, colMeans(!is.na(beta_all)) > 0, drop = FALSE]

# CpG coverage filter
keep_rows <- rowMeans(!is.na(beta_ok)) >= ROW_COVERAGE
beta_ok <- beta_ok[keep_rows, , drop = FALSE]

# Median impute per sample
for (j in seq_len(ncol(beta_ok))) {
  nas <- is.na(beta_ok[, j])
  if (any(nas)) beta_ok[nas, j] <- median(beta_ok[, j], na.rm = TRUE)
}

# Convert beta -> M (clamp)
beta_clamped <- pmin(pmax(beta_ok, 1e-6), 1 - 1e-6)
M <- log2(beta_clamped / (1 - beta_clamped))

message("CpG beta matrix dim: ", paste(dim(beta_ok), collapse = " x "))
message("CpG M matrix dim   : ", paste(dim(M), collapse = " x "))

saveRDS(beta_ok, file = file.path(OUT_CPG, "methyl_beta_matrix.rds"))
saveRDS(M,       file = file.path(OUT_CPG, "methyl_M_matrix.rds"))

# -----------------------------
# 6) Build phenotype aligned to CpG matrix columns
# -----------------------------
pheno <- meta_all %>%
  dplyr::mutate(`Sample ID` = short_barcode(`Sample ID`)) %>%
  dplyr::filter(`Sample ID` %in% colnames(M)) %>%
  dplyr::distinct(`Sample ID`, .keep_all = TRUE)

# Make sure pheno is a base data.frame (avoid tibble rowname issues)
pheno <- as.data.frame(pheno)

# Clean Sample ID (character + trim)
pheno$`Sample ID` <- trimws(as.character(pheno$`Sample ID`))

# If Tissue Type missing, infer from sample code (01 tumour; 10/11 normal)
if (!("Tissue Type" %in% colnames(pheno)) || any(is.na(pheno$`Tissue Type`))) {
  codes <- substr(pheno$`Sample ID`, 14, 15)
  pheno$`Tissue Type` <- ifelse(codes %in% c("10", "11"), "Normal", "Tumor")
}

# Align phenotype rows to M matrix columns
pheno <- pheno[match(colnames(M), pheno$`Sample ID`), , drop = FALSE]
stopifnot(all(pheno$`Sample ID` == colnames(M)))

# Set rownames to Sample ID (must match matrix colnames exactly)
rownames(pheno) <- pheno$`Sample ID`

# Add patient_id (12 chars)
pheno$patient_id <- substr(pheno$`Sample ID`, 1, 12)

message("Tissue Type table:")
print(table(pheno$`Tissue Type`, useNA = "ifany"))

saveRDS(pheno, file = file.path(OUT_CPG, "methyl_pheno.rds"))
# -----------------------------
# 7) CpG -> promoter gene matrix (TSS200/TSS1500)
# -----------------------------

anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

common_cpgs <- intersect(rownames(M), rownames(anno))
message("CpGs with 450k annotation: ", length(common_cpgs))

anno_sub <- anno[common_cpgs, c("Name", "UCSC_RefGene_Name", "UCSC_RefGene_Group")]
is_promoter <- grepl("TSS200|TSS1500", anno_sub$UCSC_RefGene_Group)
anno_prom <- anno_sub[is_promoter, , drop = FALSE]

split_genes <- strsplit(as.character(anno_prom$UCSC_RefGene_Name), ";")
mapping <- data.frame(
  CpG  = rep(anno_prom$Name, times = lengths(split_genes)),
  GENE = unlist(split_genes),
  stringsAsFactors = FALSE
)
mapping$GENE <- trimws(mapping$GENE)
mapping <- mapping[mapping$GENE != "" & !is.na(mapping$GENE), , drop = FALSE]

genes <- unique(mapping$GENE)
message("Unique promoter genes: ", length(genes))

prom_mat <- vapply(
  genes,
  function(g) {
    cpgs_g <- mapping$CpG[mapping$GENE == g]
    cpgs_g <- intersect(cpgs_g, rownames(M))
    if (length(cpgs_g) == 0) {
      rep(NA_real_, ncol(M))
    } else {
      colMeans(M[cpgs_g, , drop = FALSE], na.rm = TRUE)
    }
  },
  numeric(ncol(M))
)

prom_mat <- t(prom_mat)
rownames(prom_mat) <- genes
colnames(prom_mat) <- colnames(M)

# Remove zero-variance genes
v <- matrixStats::rowVars(prom_mat, na.rm = TRUE)
prom_mat <- prom_mat[v > 0, , drop = FALSE]

message("Promoter matrix (all samples) dim: ", paste(dim(prom_mat), collapse = " x "))

# -----------------------------
# 8) Tumour-only promoter matrix + tumour pheno (integration-ready)
# -----------------------------
tumor_ids <- rownames(pheno)[pheno$`Tissue Type` == "Tumor"]
tumor_ids <- intersect(tumor_ids, colnames(prom_mat))

message("Tumour IDs after intersect: ", length(tumor_ids))
if (length(tumor_ids) == 0) {
  message("Example pheno rownames: "); print(head(rownames(pheno), 5))
  message("Example prom_mat colnames: "); print(head(colnames(prom_mat), 5))
  stop("No tumour samples found after phenotype filtering (ID mismatch).")
}

meth_tumor  <- prom_mat[, tumor_ids, drop = FALSE]
pheno_tumor <- pheno[tumor_ids, , drop = FALSE]
stopifnot(all(rownames(pheno_tumor) == colnames(meth_tumor)))

saveRDS(meth_tumor,  file = file.path(OUT_PROM, "meth_gene_matrix_promoter_tumor.rds"))
saveRDS(pheno_tumor, file = file.path(OUT_PROM, "meth_promoter_pheno_tumor.rds"))

message("Tumour-only promoter matrix dim: ", paste(dim(meth_tumor), collapse = " x "))
message("Tumour-only pheno rows         : ", nrow(pheno_tumor))
message("Saved CpG outputs     -> ", OUT_CPG)
message("Saved promoter outputs-> ", OUT_PROM)
message("DONE.")

