# ==============================================================
# Script: 01_transcriptomics_deseq2_unpaired.R
# Purpose:
#   Differential expression (Tumor vs Normal) using DESeq2 on the
#   UNPAIRED TCGA-GBM RNA-seq dataset.
#
# NOTES (DO NOT REMOVE):
# - This script relies on a single central path config:
#     source("scripts/00_setup_paths.R")
# - All input/output folders must be defined there (no hard-coded paths here).
# - Output figures/tables are written to `tx_run_dir`.
# ==============================================================

source("scripts/00_setup_paths.R")

## =============================================================
## 0) SAFETY CHECKS (MODE + REQUIRED PATHS)
## =============================================================
if (!exists("TX_MODE")) stop("TX_MODE not found. Check scripts/00_setup_paths.R")
if (TX_MODE != "unpaired") stop("TX_MODE must be 'unpaired' for this script. Current: ", TX_MODE)

if (!exists("tx_dir"))       stop("tx_dir not found. Check scripts/00_setup_paths.R")
if (!exists("tx_unpaired"))  stop("tx_unpaired not found. Check scripts/00_setup_paths.R")
if (!exists("tx_run_dir"))   stop("tx_run_dir not found. Check scripts/00_setup_paths.R")

# Count folder is NOT tx_unpaired in your structure.
# Counts are in: data/transcriptomics/GDC_GBM_counts
tx_countsdir <- file.path(tx_dir, "GDC_GBM_counts")

if (!dir.exists(tx_countsdir)) stop("Counts folder not found: ", tx_countsdir)
if (!dir.exists(tx_unpaired))  stop("Unpaired folder (sample sheet) not found: ", tx_unpaired)

dir.create(tx_run_dir, recursive = TRUE, showWarnings = FALSE)

message("TX_MODE      : ", TX_MODE)
message("tx_dir       : ", tx_dir)
message("tx_unpaired  : ", tx_unpaired, "  (sample sheet folder)")
message("tx_countsdir : ", tx_countsdir, "  (count TSV folder)")
message("tx_run_dir   : ", tx_run_dir)

## =============================================================
## 1) LOAD LIBRARIES
## =============================================================
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(purrr)
  library(stringr)
  library(readr)
  library(tibble)
  library(DESeq2)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(ggplot2)
  library(pheatmap)
  library(grid)
})

## =============================================================
## 2) LOAD COUNT TSV FILES (UUID-named)
## =============================================================
count_files <- list.files(
  tx_countsdir,
  pattern = "\\.tsv$",
  full.names = TRUE
)

stopifnot(length(count_files) > 0)
message("Number of count TSV files: ", length(count_files))

read_one <- function(f) {
  dt <- fread(f, showProgress = FALSE)
  
  if (!("gene_id" %in% names(dt))) stop("No gene_id column found in: ", f)
  
  count_col <- intersect(
    c("unstranded", "raw_counts", "expected_count", "counts"),
    names(dt)
  )[1]
  
  if (is.na(count_col)) stop("No count column found in: ", f)
  
  tibble(
    ENSEMBL = sub("\\.\\d+$", "", dt$gene_id),
    !!basename(f) := as.integer(dt[[count_col]])
  )
}

counts_long <- purrr::reduce(
  lapply(count_files, read_one),
  dplyr::full_join,
  by = "ENSEMBL"
) %>%
  dplyr::arrange(ENSEMBL)

counts <- as.matrix(counts_long[, -1])
rownames(counts) <- counts_long$ENSEMBL

# full_join NA -> 0
if (anyNA(counts)) {
  message("WARNING: NA in merged counts. Replacing with 0.")
  counts[is.na(counts)] <- 0
}

message("Counts matrix: ", nrow(counts), " genes x ", ncol(counts), " samples/files")

## =============================================================
## 3) LOAD SAMPLE SHEET (from tx_unpaired)
## =============================================================
ss_path <- list.files(
  tx_unpaired,
  pattern = "^gdc_sample_sheet.*\\.tsv$",
  full.names = TRUE
)

stopifnot(length(ss_path) == 1)

ss <- read_tsv(ss_path, show_col_types = FALSE)
message("Sample sheet rows: ", nrow(ss))

## =============================================================
## 4) BUILD COLData BY UUID MATCH
## =============================================================
# sample sheet has "File Name" containing UUIDs
key <- ss %>%
  mutate(
    UUID = stringr::str_extract(
      `File Name`,
      "[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}"
    )
  ) %>%
  transmute(
    UUID,
    Sample_ID        = `Sample ID`,
    Tissue_Type      = `Tissue Type`,
    Tumor_Descriptor = `Tumor Descriptor`
  )

# count filenames are UUID.tsv
col_map <- tibble(file = colnames(counts)) %>%
  mutate(
    UUID = stringr::str_extract(
      file,
      "[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}"
    )
  )

j1 <- left_join(col_map, key, by = "UUID")

if (any(is.na(j1$Tissue_Type))) {
  message("ERROR: Some count files do not match sample sheet UUIDs.")
  message("Unmatched (first 20):")
  print(head(j1$file[is.na(j1$Tissue_Type)], 20))
  stop("STOP: Fix mismatch between tx_countsdir files and sample sheet 'File Name'.")
}

message("Tissue types:")
print(table(j1$Tissue_Type))

coldata <- j1 %>%
  transmute(
    sample = make.unique(ifelse(!is.na(Sample_ID), Sample_ID, UUID)),
    Tissue_Type,
    Tumor_Descriptor
  ) %>%
  as.data.frame()

rownames(coldata) <- coldata$sample
coldata$sample <- NULL

coldata$condition <- ifelse(
  grepl("Normal", coldata$Tissue_Type, ignore.case = TRUE),
  "Normal", "Tumor"
)
coldata$condition <- factor(coldata$condition, levels = c("Normal", "Tumor"))

message("Condition table:")
print(table(coldata$condition))

## =============================================================
## 5) ALIGN COUNTS <-> COLData
## =============================================================
stopifnot(ncol(counts) == nrow(coldata))

colnames(counts) <- rownames(coldata)
coldata <- coldata[colnames(counts), , drop = FALSE]

stopifnot(identical(colnames(counts), rownames(coldata)))
stopifnot(all(!is.na(coldata$condition)))

## =============================================================
## 6) DESEQ2 (Tumor vs Normal)
## =============================================================
dds <- DESeqDataSetFromMatrix(
  countData = round(counts),
  colData   = coldata,
  design    = ~ condition
)

keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
message("After prefilter: ", nrow(dds), " genes retained")

dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "Tumor", "Normal"))
print(summary(res))

res_df <- as.data.frame(res) %>%
  tibble::rownames_to_column("ENSEMBL") %>%
  dplyr::mutate(
    SYMBOL = AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys      = ENSEMBL,
      column    = "SYMBOL",
      keytype   = "ENSEMBL",
      multiVals = "first"
    )
  )

# FDR<0.05 için up/down sayımı
deg_counts <- res_df %>%
  dplyr::filter(!is.na(padj), padj < 0.05) %>%
  dplyr::summarise(
    DEGs_total = dplyr::n(),
    Up_in_Tumor = sum(log2FoldChange > 0, na.rm = TRUE),
    Down_in_Tumor = sum(log2FoldChange < 0, na.rm = TRUE)
  )

print(deg_counts)
# Save results
write.csv(res_df, file.path(tx_run_dir, "DESeq2_results_unpaired.csv"), row.names = FALSE)

## =============================================================
## 7) QC TRANSFORM (VST)
## =============================================================
vsd <- vst(dds, blind = FALSE)

## ==============================================================
## 8) PLOTS (show + save)
## ==============================================================
# Volcano
volc <- res_df %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::mutate(class = dplyr::case_when(
    padj < 0.05 & log2FoldChange >=  1 ~ "Up",
    padj < 0.05 & log2FoldChange <= -1 ~ "Down",
    TRUE                               ~ "NS"
  ))

p_volc <- ggplot(volc, aes(x = log2FoldChange, y = -log10(padj), color = class)) +
  geom_point(alpha = 0.7, size = 1.3) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  scale_color_manual(values = c(Up = "red", Down = "blue", NS = "lightgreen")) +
  labs(
    x = "log2 fold-change (Tumor vs Normal)",
    y = "-log10(FDR-adjusted p-value)",
    color = "",
    title = "Volcano plot: GBM Tumor vs Normal (unpaired)"
  ) +
  theme_minimal(base_size = 12)

print(p_volc)
ggsave(file.path(tx_run_dir, "volcano_plot_unpaired.png"), p_volc, width = 7, height = 5, dpi = 300)

# PCA
p_pca <- plotPCA(vsd, intgroup = "condition") +
  ggtitle("PCA (VST): GBM Tumor vs Normal (unpaired)") +
  theme_minimal()

print(p_pca)
ggsave(file.path(tx_run_dir, "PCA_VST_unpaired.png"), p_pca, width = 6, height = 5, dpi = 300)

# Heatmap: Top 50 by padj (unique SYMBOL)
top_sym <- res_df %>%
  dplyr::filter(!is.na(padj), !is.na(SYMBOL), SYMBOL != "") %>%
  dplyr::arrange(padj) %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE) %>%
  dplyr::slice_head(n = 50)

mat <- assay(vsd)[top_sym$ENSEMBL, , drop = FALSE]
ord <- order(coldata$condition)
mat <- mat[, ord, drop = FALSE]

ann <- data.frame(
  condition = factor(
    coldata$condition[ord],
    levels = c("Normal", "Tumor")
  )
)
rownames(ann) <- colnames(mat)

ann_colors <- list(
  condition = c(
    Normal = "blue",
    Tumor  = "red"
  )
)


matz <- t(scale(t(mat)))

heat <- pheatmap(
  matz,
  annotation_col = ann,
  annotation_colors = ann_colors,
  show_rownames = TRUE,
  show_colnames = TRUE,
  labels_row = top_sym$SYMBOL,
  fontsize_row = 5,
  fontsize_col = 7,
  angle_col = 45,
  border_color = NA,
  main = "Top 50 DE genes (VST z-score)"
)

png(file.path(tx_run_dir, "heatmap_top50_DE_unpaired.png"), width = 1400, height = 1100, res = 150)
grid::grid.draw(heat$gtable)
dev.off()

message("DONE. Outputs written to: ", tx_run_dir)
