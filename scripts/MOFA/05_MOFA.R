## =========================================================
## MOFA.R — Full run (clean paths + consistent filenames)
## Inputs  (auto-detected) : results/paired_inputs/*.rds
## Outputs (fixed names)   : results/paired_inputs/{rna_mofa_input.rds, meth_mofa_input.rds, sample_meta_mofa.rds}
##                         : results/MOFA/run_paired01/*
## =========================================================

## =========================
## STEP 0) PATHS (project-safe)
## =========================

## ---------------------------
## 0A) PATHS (project-safe; ONLY from 00_setup_paths.R)
## ---------------------------
source("scripts/00_setup_paths.R")

# Canonicalize (use setup aliases if needed)
if (!exists("ROOT") && exists("project_dir")) ROOT <- project_dir
if (!exists("RESULTS_DIR") && exists("results_dir")) RESULTS_DIR <- results_dir
if (!exists("RESULTS_DIR")) stop("RESULTS_DIR/results_dir not found. Check scripts/00_setup_paths.R")

RUN_ID <- "run_paired01"

## ---- Inputs: MUST be setup-defined paired_inputs_dir ----
if (!exists("paired_inputs_dir")) stop("paired_inputs_dir not found. Check scripts/00_setup_paths.R")
IN_DIR <- paired_inputs_dir
stopifnot(dir.exists(IN_DIR))

## ---- Outputs: MUST be setup-defined mofa_dir (results/MOFA) ----
if (!exists("mofa_dir")) stop("mofa_dir not found. Check scripts/00_setup_paths.R")

OUT_DIR   <- file.path(mofa_dir, RUN_ID)
PLOTS_DIR <- file.path(OUT_DIR, "plots")

dir.create(OUT_DIR,   recursive = TRUE, showWarnings = FALSE)
dir.create(PLOTS_DIR, recursive = TRUE, showWarnings = FALSE)

run_dir   <- OUT_DIR
plots_dir <- PLOTS_DIR

message("ROOT:         ", ROOT)
message("RESULTS_DIR:  ", RESULTS_DIR)
message("IN_DIR:       ", IN_DIR)
message("OUT_DIR:      ", OUT_DIR)
message("PLOTS_DIR:    ", PLOTS_DIR)

stopifnot(dir.exists(IN_DIR))

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

## Helper: find first existing file by patterns (in a given dir)
pick_file <- function(dir, patterns) {
  files <- list.files(dir, full.names = TRUE)
  for (p in patterns) {
    hit <- files[grepl(p, basename(files), ignore.case = TRUE)]
    if (length(hit) > 0) return(hit[1])
  }
  return(NA_character_)
}

to_num_matrix <- function(x) {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x)) stop("Input is not a matrix/data.frame.")
  storage.mode(x) <- "double"
  x
}

## =========================================================
## STEP 1) Load + align paired RNA/METH inputs
## Save MOFA-ready inputs into results/paired_inputs (fixed names)
## =========================================================

## 1A) Auto-detect your existing paired matrices/meta INSIDE IN_DIR
rna_mat_path <- pick_file(IN_DIR, c(
  "^rna_.*mofa.*\\.rds$",
  "^rna_logcounts_paired_SYMBOL\\.rds$",
  "^rna_logcounts_paired\\.rds$",
  "^rna_.*paired.*\\.rds$"
))

meth_mat_path <- pick_file(IN_DIR, c(
  "^meth_.*mofa.*\\.rds$",
  "^meth_.*promoter_.*paired.*\\.rds$",
  "^meth_.*paired.*\\.rds$"
))

rna_meta_path <- pick_file(IN_DIR, c(
  "^rna_metadata_paired\\.rds$",
  "^rna_meta.*paired.*\\.rds$"
))

meth_meta_path <- pick_file(IN_DIR, c(
  "^meth_metadata_paired\\.rds$",
  "^meth_meta.*paired.*\\.rds$"
))

message("RNA matrix : ", rna_mat_path)
message("METH matrix: ", meth_mat_path)
message("RNA meta   : ", rna_meta_path)
message("METH meta  : ", meth_meta_path)

if (is.na(rna_mat_path))  stop("RNA matrix not found in: ", IN_DIR)
if (is.na(meth_mat_path)) stop("METH matrix not found in: ", IN_DIR)

## 1B) Load matrices
rna_mat  <- to_num_matrix(readRDS(rna_mat_path))
meth_mat <- to_num_matrix(readRDS(meth_mat_path))

message("Loaded RNA  dim: ", paste(dim(rna_mat),  collapse = " x "))
message("Loaded METH dim: ", paste(dim(meth_mat), collapse = " x "))

## 1C) Load metadata (optional)
rna_meta  <- if (!is.na(rna_meta_path))  readRDS(rna_meta_path)  else NULL
meth_meta <- if (!is.na(meth_meta_path)) readRDS(meth_meta_path) else NULL

## Prefer patient_id if exists; otherwise use colnames
get_sample_key <- function(meta, mat) {
  if (!is.null(meta) && "patient_id" %in% colnames(meta)) return(as.character(meta$patient_id))
  return(colnames(mat))
}

rna_key  <- get_sample_key(rna_meta,  rna_mat)
meth_key <- get_sample_key(meth_meta, meth_mat)

# Force matrix colnames to keys (only if lengths match)
if (!is.null(rna_key)  && length(rna_key)  == ncol(rna_mat))  colnames(rna_mat)  <- rna_key
if (!is.null(meth_key) && length(meth_key) == ncol(meth_mat)) colnames(meth_mat) <- meth_key

## 1D) Align samples (intersection) + reorder
common_samples <- intersect(colnames(rna_mat), colnames(meth_mat))
if (length(common_samples) == 0) {
  stop("No overlapping samples between RNA and METH (check colnames/patient_id).")
}

common_samples <- sort(common_samples)
rna_mat  <- rna_mat[,  common_samples, drop = FALSE]
meth_mat <- meth_mat[, common_samples, drop = FALSE]

message("Common samples (paired n): ", length(common_samples))

## 1E) Align features only if gene IDs overlap (optional)
common_genes <- intersect(rownames(rna_mat), rownames(meth_mat))
message("Common genes: ", length(common_genes))

if (length(common_genes) > 0) {
  rna_mat  <- rna_mat[ common_genes, , drop = FALSE]
  meth_mat <- meth_mat[common_genes, , drop = FALSE]
  stopifnot(identical(rownames(rna_mat), rownames(meth_mat)))
} else {
  warning("0 common genes. OK only if your views intentionally have different feature spaces. ",
          "If unintended, fix gene ID mapping (SYMBOL vs ENSG etc.).")
}

## 1F) Build sample meta (minimal)
sample_meta <- data.frame(sample_id = common_samples, stringsAsFactors = FALSE)

if (!is.null(rna_meta) && all(c("patient_id","barcode") %in% colnames(rna_meta))) {
  tmp <- rna_meta %>%
    mutate(patient_id = as.character(patient_id)) %>%
    distinct(patient_id, .keep_all = TRUE)
  sample_meta$barcode <- tmp$barcode[match(sample_meta$sample_id, tmp$patient_id)]
}

## 1G) Save MOFA-ready inputs (FIXED NAMES) into IN_DIR
OUT_RNA  <- file.path(IN_DIR, "rna_mofa_input.rds")
OUT_METH <- file.path(IN_DIR, "meth_mofa_input.rds")
OUT_META <- file.path(IN_DIR, "sample_meta_mofa.rds")

saveRDS(rna_mat,      OUT_RNA)
saveRDS(meth_mat,     OUT_METH)
saveRDS(sample_meta,  OUT_META)

message("Saved: ", OUT_RNA,  "  dim = ", paste(dim(rna_mat),  collapse = " x "))
message("Saved: ", OUT_METH, "  dim = ", paste(dim(meth_mat), collapse = " x "))
message("Saved: ", OUT_META, "  n = ", nrow(sample_meta))
message("STEP 1 DONE" )

## =========================================================
## STEP 2) Fit MOFA model
## Reads fixed inputs from IN_DIR; saves mofa_fit.rds to OUT_DIR
## =========================================================

suppressPackageStartupMessages({
  library(MOFA2)
})

rna_mat  <- readRDS(file.path(IN_DIR, "rna_mofa_input.rds"))
meth_mat <- readRDS(file.path(IN_DIR, "meth_mofa_input.rds"))
meta     <- readRDS(file.path(IN_DIR, "sample_meta_mofa.rds"))

stopifnot(
  is.matrix(rna_mat),
  is.matrix(meth_mat),
  ncol(rna_mat) == ncol(meth_mat)
)

message("RNA dim  : ", paste(dim(rna_mat), collapse = " x "))
message("METH dim : ", paste(dim(meth_mat), collapse = " x "))

data_list <- list(RNA = rna_mat, METH = meth_mat)
MOFAobject <- create_mofa(data_list)

# Options (same as yours)
data_opts <- get_default_data_options(MOFAobject)
data_opts$scale_views  <- TRUE
data_opts$scale_groups <- FALSE

model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 15

train_opts <- get_default_training_options(MOFAobject)
train_opts$seed <- 42
train_opts$maxiter <- 1500
train_opts$convergence_mode <- "medium"
train_opts$drop_factor_threshold <- 0.01

MOFAobject <- prepare_mofa(
  MOFAobject,
  data_options     = data_opts,
  model_options    = model_opts,
  training_options = train_opts
)

MOFA_fit <- run_mofa(MOFAobject, use_basilisk = TRUE)

saveRDS(MOFA_fit, file.path(OUT_DIR, "mofa_fit.rds"))
message("STEP 2 DONE  Saved: ", file.path(OUT_DIR, "mofa_fit.rds"))

## Quick check: factor count
facs <- MOFA2::get_factors(MOFA_fit, groups = "all", as.data.frame = FALSE)[[1]]
cat("Number of factors in MOFA_fit:", ncol(facs), "\n")
print(colnames(facs))

## =========================================================
## STEP 3) Extract Z (scores) and W (weights) + save (fixed names)
## =========================================================

# Z (samples x K)
Z_list <- MOFA2::get_factors(MOFA_fit, groups = "all", factors = "all", as.data.frame = FALSE)
Z_all <- Z_list[[1]]

# W (features x K) per view
W_list <- MOFA2::get_weights(MOFA_fit, views = "all", factors = "all", as.data.frame = FALSE)
W_RNA  <- as.matrix(W_list[["RNA"]])
W_METH <- as.matrix(W_list[["METH"]])

cat("Z_all dim (samples x K):", paste(dim(Z_all), collapse = " x "), "\n")
cat("W_RNA  dim (features x K):", paste(dim(W_RNA),  collapse = " x "), "\n")
cat("W_METH dim (features x K):", paste(dim(W_METH), collapse = " x "), "\n")

stopifnot(ncol(Z_all) == ncol(W_RNA), ncol(Z_all) == ncol(W_METH))

saveRDS(Z_all,  file.path(OUT_DIR, "Z_all.rds"))
saveRDS(W_RNA,  file.path(OUT_DIR, "W_RNA_all.rds"))
saveRDS(W_METH, file.path(OUT_DIR, "W_METH_all.rds"))

message("STEP 3 DONE Saved Z/W to OUT_DIR")

## =========================================================
## STEP 4) Variance explained + rank factors + save ranked objects (fixed names)
## =========================================================

ve_raw <- MOFA2::calculate_variance_explained(MOFA_fit)$r2_per_factor[[1]]
stopifnot(all(c("RNA","METH") %in% colnames(ve_raw)))

# Convert to percent if in [0,1]
ve_use <- if (max(ve_raw, na.rm = TRUE) <= 1) ve_raw * 100 else ve_raw

ve_df <- data.frame(
  factor_id = rownames(ve_use),
  VE_RNA    = ve_use[, "RNA"],
  VE_METH   = ve_use[, "METH"],
  stringsAsFactors = FALSE
)
ve_df$VE_TOTAL <- ve_df$VE_RNA + ve_df$VE_METH

# Rank by total VE
ve_df <- ve_df[order(ve_df$VE_TOTAL, decreasing = TRUE), , drop = FALSE]
ve_df$Rank   <- seq_len(nrow(ve_df))
ve_df$Factor <- paste0("F", ve_df$Rank)

# Reorder Z/W to ranked factor order
top_factor_ids <- ve_df$factor_id

Z_ranked      <- Z_all[,  top_factor_ids, drop = FALSE]
W_RNA_ranked  <- W_RNA[,  top_factor_ids, drop = FALSE]
W_METH_ranked <- W_METH[, top_factor_ids, drop = FALSE]

colnames(Z_ranked)      <- ve_df$Factor
colnames(W_RNA_ranked)  <- ve_df$Factor
colnames(W_METH_ranked) <- ve_df$Factor

stopifnot(
  ncol(Z_ranked) == ncol(W_RNA_ranked),
  ncol(Z_ranked) == ncol(W_METH_ranked)
)

saveRDS(ve_df,         file.path(OUT_DIR, "VE_per_factor_ranked.rds"))
write.csv(ve_df,       file.path(OUT_DIR, "VE_per_factor_ranked.csv"), row.names = FALSE)

saveRDS(Z_ranked,      file.path(OUT_DIR, "Z_ranked_F1toFK.rds"))
saveRDS(W_RNA_ranked,  file.path(OUT_DIR, "W_RNA_ranked_F1toFK.rds"))
saveRDS(W_METH_ranked, file.path(OUT_DIR, "W_METH_ranked_F1toFK.rds"))

message("STEP 4 DONE Saved ranked outputs to OUT_DIR")

## =========================================================
## STEP 5) PLOTS: VE barplot + Z heatmap (fixed filenames)
## =========================================================

suppressPackageStartupMessages({
  library(tidyr)
  library(ggplot2)
  library(forcats)
  library(pheatmap)
})

ve_df    <- readRDS(file.path(OUT_DIR, "VE_per_factor_ranked.rds"))
Z_ranked <- readRDS(file.path(OUT_DIR, "Z_ranked_F1toFK.rds"))
Z_ranked <- as.matrix(Z_ranked)

# 5A) Variance explained stacked barplot
ve_long <- ve_df |>
  select(Factor, VE_RNA, VE_METH) |>
  pivot_longer(cols = c(VE_RNA, VE_METH), names_to = "View", values_to = "Percent") |>
  mutate(
    View = recode(View, VE_RNA = "RNA expression", VE_METH = "DNA methylation"),
    Factor = fct_relevel(Factor, ve_df$Factor)
  )

p_ve <- ggplot(ve_long, aes(x = Factor, y = Percent, fill = View)) +
  geom_col() +
  scale_fill_manual(
    values = c(
      "RNA expression"= "#1B4F72" ,
      "DNA methylation" = "#A7C7E7"
    ),
    name = "View"
  ) +
  labs(
    title = "MOFA2: % variance explained per factor (by view)",
    y = "% variance explained",
    x = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title  = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )

print(p_ve)

ggsave(
  filename = file.path(PLOTS_DIR, "MOFA_VE_stacked_bluepalette.png"),
  plot     = p_ve,
  width    = 8,
  height   = 4,
  dpi      = 300,
  bg       = "white"
)

## ========== 5B) FACTOR CLASSES & VE PLOTS (RESULTS, Top-20) ==========

# Define TOPN here (Top-20 as in your old workflow)
TOPN <- min(20L, nrow(ve_df))
ve_top <- ve_df[seq_len(TOPN), , drop = FALSE]

# --- Thresholds (unchanged logic) ---
weak_thr <- 2.0   # <2% total VE
spec_thr <- 0.80  # >=80% RNA share => RNA-specific

# --- RNA share (safe division) ---
ve_top$RNA_share <- ve_top$VE_RNA / pmax(ve_top$VE_TOTAL, 1e-9)

# --- Class assignment (robust, no recycling issues) ---
ve_top$Class <- "Shared"
ve_top$Class[ve_top$VE_TOTAL < weak_thr] <- "Weak (<2%)"
ve_top$Class[ve_top$RNA_share >= spec_thr] <- "RNA-specific"
ve_top$Class[ve_top$RNA_share <= (1 - spec_thr)] <- "METH-specific"

ve_top$Class <- factor(
  ve_top$Class,
  levels = c("METH-specific","RNA-specific","Shared","Weak (<2%)")
)

# --- Save table ---
out_csv <- file.path(run_dir, "VE_classes_top15_MOFA.csv")
write_csv(ve_top, out_csv)

# --- Console summary ---
message("Saved factor class table: ", out_csv)
print(table(ve_top$Class))

## ====== VE by factor class (stacked RNA vs METH) ======

agg_class <- ve_top |>
  group_by(Class) |>
  summarise(
    RNA  = sum(VE_RNA,  na.rm = TRUE),
    METH = sum(VE_METH, na.rm = TRUE),
    .groups = "drop"
  ) |>
  pivot_longer(
    cols = c(RNA, METH),
    names_to = "View",
    values_to = "Percent"
  ) |>
  mutate(View = factor(View, levels = c("RNA", "METH")))

p_class <- ggplot(agg_class, aes(x = Class, y = Percent, fill = View)) +
  geom_col(position = "stack") +
  scale_fill_manual(
    values = c(RNA = "#0b3c5d", METH = "#cfe1f2"),
    name   = "View"
  ) +
  labs(
    title = "MOFA2: Explained variance by factor class (by-view)",
    y     = "Total % variance explained",
    x     = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title      = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom",
    legend.title    = element_text(face = "bold")
  )

print(p_class)

ggsave(
  filename = file.path(plots_dir, "MOFA_VE_by_class_top20.png"),
  plot     = p_class,
  width    = 7.2,
  height   = 4.2,
  dpi      = 300,
  bg       = "white"
)

## ====== VE stacked barplot per factor (Top-20) ======

ve_long <- ve_top |>
  transmute(
    Factor,
    `RNA expression`  = VE_RNA,
    `DNA methylation` = VE_METH
  ) |>
  pivot_longer(-Factor, names_to = "View", values_to = "Percent") |>
  mutate(Factor = fct_relevel(Factor, ve_top$Factor))

p_ve <- ggplot(ve_long, aes(Factor, Percent, fill = View)) +
  geom_col() +
  scale_fill_manual(
    values = c(`RNA expression`  = "#0b3c5d",
               `DNA methylation` = "#cfe1f2"),
    name = "View"
  ) +
  labs(
    title = "MOFA2: % variance explained per factor (by view)",
    y = "% variance explained",
    x = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title  = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )

print(p_ve)

ggsave(
  filename = file.path(plots_dir, "MOFA_VE_stacked_top20.png"),
  plot     = p_ve,
  width    = 8,
  height   = 4,
  dpi      = 300,
  bg       = "white"
)

## =========================================================
## STEP 6) FACTOR SCORE HEATMAP & DOMINANCE (Top-20) — continue
## =========================================================

suppressPackageStartupMessages({
  library(pheatmap)
  library(grid)
})

# ---- safety: output dir ----
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

# ---- align columns by factor names (prevents silent mismatch) ----
common_factors <- intersect(colnames(Z_topN), colnames(W_RNA_topN))
common_factors <- intersect(common_factors, colnames(W_METH_topN))

stopifnot(length(common_factors) > 0)

Z_topN_use      <- Z_topN[, common_factors, drop = FALSE]
W_RNA_topN_use  <- W_RNA_topN[, common_factors, drop = FALSE]
W_METH_topN_use <- W_METH_topN[, common_factors, drop = FALSE]

# ---- dominance scores ----
rna_act  <- colSums(W_RNA_topN_use^2,  na.rm = TRUE)
meth_act <- colSums(W_METH_topN_use^2, na.rm = TRUE)
tot_act  <- rna_act + meth_act
ratio    <- rna_act / pmax(tot_act, 1e-12)

dom_class <- ifelse(
  ratio >= 0.65, "RNA-dominant",
  ifelse(ratio <= 0.35, "METH-dominant", "Balanced")
)

# ---- order by total activity (same as your code) ----
ord <- order(tot_act, decreasing = TRUE)

Z_plot    <- Z_topN_use[, ord, drop = FALSE]
rna_act   <- rna_act[ord]
meth_act  <- meth_act[ord]
tot_act   <- tot_act[ord]
dom_class <- factor(dom_class[ord],
                    levels = c("RNA-dominant","Balanced","METH-dominant"))
colnames(Z_plot) <- colnames(Z_topN_use)[ord]

# ---- z-score + clip ----
Zz <- scale(Z_plot)

clip <- function(m, lim = 2) {
  m[m >  lim] <-  lim
  m[m < -lim] <- -lim
  m
}
Zz <- clip(Zz, lim = 2)

# ---- annotation (same content) ----
ann_col <- data.frame(
  Dominance    = dom_class,
  RNA_contrib  = round(rna_act  / pmax(tot_act, 1e-12), 2),
  METH_contrib = round(meth_act / pmax(tot_act, 1e-12), 2)
)
rownames(ann_col) <- colnames(Zz)

# ---- build heatmap (silent=TRUE gives gtable) ----
hp <- pheatmap(
  Zz,
  color = colorRampPalette(c("#2c7fb8", "#ffffbf", "#d95f0e"))(100),
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  clustering_method = "ward.D2",
  show_rownames = FALSE,
  show_colnames = TRUE,
  angle_col = 90,
  main = "MOFA2 factor scores (Top-15, z-scored)",
  annotation_col = ann_col,
  fontsize = 10,
  border_color = NA,
  silent = TRUE
)

# ---- show in RStudio (correct draw) ----
grid::grid.newpage()
grid::grid.draw(hp$gtable)

# ---- save PNG (bigger so it doesn't look squashed) ----
out_png <- file.path(plots_dir, "MOFA_factor_scores_15.png")
png(out_png, width = 3200, height = 1800, res = 300, bg = "white")
grid::grid.newpage()
grid::grid.draw(hp$gtable)
dev.off()

message("Saved heatmap: ", out_png)

## =========================================================
## STEP 7) Cross-omics correlation:
##        RNA-specific vs METH-specific MOFA2 factors
## =========================================================

## 1) Extract factor class vector (Factor -> Class)
class_vec <- as.character(ve_top$Class)
names(class_vec) <- ve_top$Factor

## 2) Keep only factors present in Z_topN
common_factors <- intersect(names(class_vec), colnames(Z_topN))
class_vec <- class_vec[common_factors]

## 3) Identify RNA- and METH-specific factors
rna_factors  <- names(class_vec)[class_vec == "RNA-specific"]
meth_factors <- names(class_vec)[class_vec == "METH-specific"]

## Safety check
stopifnot(length(rna_factors) > 0)
stopifnot(length(meth_factors) > 0)

## 4) Extract factor score matrices
Z_rna  <- Z_topN[, rna_factors,  drop = FALSE]
Z_meth <- Z_topN[, meth_factors, drop = FALSE]

## 5) Compute Pearson correlation matrix
cor_mat <- cor(
  Z_rna,
  Z_meth,
  method = "pearson",
  use    = "pairwise.complete.obs"
)

## 6) Add informative labels
rownames(cor_mat) <- paste0("RNA_",  colnames(Z_rna))
colnames(cor_mat) <- paste0("METH_", colnames(Z_meth))

## ---- Safety: close any open graphics devices ----
while (!is.null(dev.list())) dev.off()

## ---- Title (shorter avoids clipping) ----
ttl <- "Cross-omics correlation: RNA- vs METH-specific MOFA2 factors"

## ---- Build heatmap object (silent=TRUE is critical) ----
hp <- pheatmap::pheatmap(
  cor_mat,
  color = colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(100),
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "ward.D2",
  border_color = NA,
  fontsize = 11,
  main = ttl,
  silent = TRUE
)
hp

## ---- 1) Force draw in RStudio Plot panel ----
grid::grid.newpage()
grid::grid.draw(hp$gtable)

## ---- 2) Save the exact same plot to file ----
out_png <- file.path(plots_dir, "MOFA_cross_omics_RNA_vs_METH.png")

png(filename = out_png, width = 2600, height = 1700, res = 300, bg = "white")
grid::grid.newpage()
grid::grid.draw(hp$gtable)
dev.off()

message("Saved Figure 14: ", out_png)

#====================================================================
# STEP 8: GSEA/ ORA
#=====================================================================

# libraries:
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(readr)
})

#8.1 Helper: SYMBOL → ENTREZ + ranked vector (GSEA)

sym2ent <- function(sym) {
  sym <- unique(na.omit(as.character(sym)))
  if (length(sym) == 0) return(character(0))
  
  res <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = sym,
    keytype = "SYMBOL",
    column = "ENTREZID",
    multiVals = "first"
  )
  res[!is.na(res)]
}

prep_gsea_rank <- function(w_named) {
  ent <- sym2ent(names(w_named))
  keep <- names(w_named) %in% names(ent)
  ranks <- w_named[keep]
  names(ranks) <- ent[names(ranks)]
  ranks <- ranks[!duplicated(names(ranks))]
  sort(ranks, decreasing = TRUE)
}

#8.2 GSEA – RNA factors

gsea_rna_dir <- file.path(run_dir, "GSEA_RNA")
dir.create(gsea_rna_dir, showWarnings = FALSE)

for (k in seq_len(ncol(W_RNA_topN))) {
  
  w <- W_RNA_topN[, k]
  names(w) <- rownames(W_RNA_topN)
  
  ## --- FIX 1: remove '_RNA' suffix ---
  gene_sym <- sub("_RNA$", "", names(w))
  
  ## --- FIX 2: map SYMBOL -> ENTREZ ---
  entrez <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys     = gene_sym,
    keytype  = "SYMBOL",
    column   = "ENTREZID",
    multiVals = "first"
  )
  
  ## --- build ranked list ---
  keep <- !is.na(entrez)
  ranks <- w[keep]
  names(ranks) <- entrez[keep]
  
  ## remove duplicated ENTREZ IDs
  ranks <- ranks[!duplicated(names(ranks))]
  ranks <- sort(ranks, decreasing = TRUE)
  
  if (length(ranks) < 50) next
  
  gse <- tryCatch(
    gseGO(
      geneList      = ranks,
      OrgDb         = org.Hs.eg.db,
      keyType       = "ENTREZID",
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      verbose       = FALSE
    ),
    error = function(e) NULL
  )
  
  if (!is.null(gse) && nrow(as.data.frame(gse)) > 0) {
    df <- as.data.frame(gse)
    df$Factor <- colnames(W_RNA_topN)[k]
    
    write_csv(
      df,
      file.path(
        gsea_rna_dir,
        paste0("MOFA_RNA_GSEA_", df$Factor[1], ".csv")
      )
    )
  }
}

message("RNA GSEA finished successfully.")

## =========================================================
## STEP 9) Functional enrichment (MOFA2) — GSEA (RNA) + ORA (METH)
## Outputs:
##   OUT_DIR/GSEA_RNA/   (per-factor + combined tables)
##   OUT_DIR/ORA_METH/   (per-factor + combined tables)
##   OUT_DIR/results/enrichment/ (summary count tables)
## =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(clusterProfiler)
})

## ---------------------------
## 9.0) Safety checks
## ---------------------------
RNA_GSEA_DIR  <- file.path(OUT_DIR, "GSEA_RNA")
METH_GSEA_DIR <- file.path(OUT_DIR, "GSEA_METH")
METH_ORA_DIR  <- file.path(OUT_DIR, "ORA_METH")

dir.create(RNA_GSEA_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(METH_GSEA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(METH_ORA_DIR,  recursive = TRUE, showWarnings = FALSE)

message("RNA_GSEA_DIR:  ", RNA_GSEA_DIR)
message("METH_GSEA_DIR: ", METH_GSEA_DIR)
message("METH_ORA_DIR:  ", METH_ORA_DIR)

## ---------------------------
## 9.1) Helpers (robust SYMBOL -> ENTREZ)
## ---------------------------

# Clean feature names to valid HGNC symbols (remove suffixes like "_RNA", "_METH", etc.)
clean_symbol <- function(x) {
  x <- as.character(x)
  
  # Remove common suffixes you created earlier (e.g., TSPAN6_RNA -> TSPAN6)
  x <- sub("_(RNA|METH)$", "", x, perl = TRUE)
  
  # If ENSEMBL IDs appear, remove version (ENSG... .12 -> ENSG...)
  x <- sub("\\.\\d+$", "", x, perl = TRUE)
  
  # Trim whitespace
  x <- trimws(x)
  
  x
}

# Map SYMBOLs to ENTREZ IDs (case-insensitive), returns named vector (original -> ENTREZ)
sym2ent <- function(sym) {
  sym <- unique(na.omit(as.character(sym)))
  if (length(sym) == 0L) return(character(0))
  
  sym_clean <- clean_symbol(sym)
  sym_up <- toupper(sym_clean)
  
  # Fast validity check against OrgDb keys
  valid_keys <- AnnotationDbi::keys(org.Hs.eg.db, keytype = "SYMBOL")
  keep <- sym_up %in% valid_keys
  
  if (!any(keep)) {
    warning("sym2ent: no valid SYMBOL keys after cleaning (check suffixes / ID types).")
    out <- rep(NA_character_, length(sym))
    names(out) <- sym
    return(out)
  }
  
  # Map only valid symbols to reduce warnings
  mapped <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys      = unique(sym_up[keep]),
    keytype   = "SYMBOL",
    column    = "ENTREZID",
    multiVals = "first"
  )
  
  # Re-expand to original order (named by original input)
  out <- mapped[match(sym_up, names(mapped))]
  names(out) <- sym
  out
}

# Prepare ranked gene list for GSEA (named numeric vector: ENTREZ -> weight)
prep_gsea_rank <- function(w_named) {
  stopifnot(!is.null(names(w_named)))
  
  ent <- sym2ent(names(w_named))
  keep <- !is.na(ent) & !duplicated(ent)
  
  if (!any(keep)) {
    warning("prep_gsea_rank: no valid ENTREZ IDs after mapping.")
    return(numeric(0))
  }
  
  ranks <- stats::setNames(as.numeric(w_named[keep]), ent[keep])
  ranks <- ranks[order(ranks, decreasing = TRUE)]
  ranks
}

# Ensure weight matrices are numeric and have rownames
to_num_with_names <- function(W) {
  W <- as.matrix(W)
  storage.mode(W) <- "double"
  if (is.null(rownames(W))) stop("Weights matrix has no rownames.")
  rn <- rownames(W)
  rn[!nzchar(rn)] <- paste0("feat_", which(!nzchar(rn)))
  rownames(W) <- rn
  W
}

W_RNA_topN  <- to_num_with_names(W_RNA_topN)
W_METH_topN <- to_num_with_names(W_METH_topN)

## ---------------------------
## 9.2) GSEA (RNA view) — GO:BP
## ---------------------------

rna_gsea_summary <- list()

for (k in seq_len(ncol(W_RNA_topN))) {
  
  w <- W_RNA_topN[, k]
  names(w) <- rownames(W_RNA_topN)
  
  ranks <- prep_gsea_rank(w)
  
  # Skip if too few mapped genes
  if (length(ranks) < 50) {
    message("MOFA RNA GSEA ", colnames(W_RNA_topN)[k], ": <50 ranked genes after mapping, skipped.")
    next
  }
  
  gse <- tryCatch(
    clusterProfiler::gseGO(
      geneList      = ranks,
      OrgDb         = org.Hs.eg.db,
      ont           = "BP",
      keyType       = "ENTREZID",
      minGSSize     = 10,
      maxGSSize     = 500,
      pvalueCutoff  = 0.05,
      pAdjustMethod = "BH",
      verbose       = FALSE
    ),
    error = function(e) NULL
  )
  
  if (!is.null(gse) && nrow(as.data.frame(gse)) > 0) {
    df <- as.data.frame(gse)
    df$Factor <- colnames(W_RNA_topN)[k]
    
    # Save per-factor table
    readr::write_csv(
      df,
      file.path(RNA_GSEA_DIR, sprintf("MOFA_RNA_GSEA_%s.csv", colnames(W_RNA_topN)[k]))
    )
    
    rna_gsea_summary[[colnames(W_RNA_topN)[k]]] <- df
  } else {
    message("MOFA RNA GSEA ", colnames(W_RNA_topN)[k], ": no significant terms at FDR <= 0.05.")
  }
}

names(rna_gsea_summary)

# Combine all RNA GSEA results
if (length(rna_gsea_summary) > 0) {
  
  ## 9A) Combine all factor tables (keep Factor name from list)
  rna_gsea_all <- dplyr::bind_rows(rna_gsea_summary, .id = "Factor")
  
  # Save ALL results
  readr::write_csv(
    rna_gsea_all,
    file.path(RNA_GSEA_DIR, "MOFA_RNA_GSEA_ALL.csv")
  )
  
  ## 9B) Filter significant terms at FDR <= 0.05
  rna_gsea_sig <- rna_gsea_all |>
    dplyr::mutate(p.adjust = suppressWarnings(as.numeric(p.adjust))) |>
    dplyr::filter(!is.na(p.adjust), p.adjust <= 0.05)
  
  readr::write_csv(
    rna_gsea_sig,
    file.path(RNA_GSEA_DIR, "MOFA_RNA_GSEA_FDR05.csv")
  )
  
  ## 9C) Count significant terms per factor
  rna_gsea_counts <- rna_gsea_sig |>
    dplyr::count(Factor, name = "sig_terms") |>
    dplyr::arrange(dplyr::desc(sig_terms))
  
  readr::write_csv(
    rna_gsea_counts,
    file.path(RNA_GSEA_DIR, "MOFA_RNA_GSEA_sig_counts.csv")
  )
  
  message("RNA GSEA done Factors with significant terms: ", nrow(rna_gsea_counts))
  
} else {
  message("RNA GSEA done No factors produced results.")
}

# plot===

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(forcats)
  library(readr)
})

## Safety: need rna_gsea_counts from Step 9C
stopifnot(exists("rna_gsea_counts"))

## Make factor order nice (F1, F2, ..., F9)
rna_gsea_counts_plot <- rna_gsea_counts |>
  dplyr::mutate(
    Factor_num = suppressWarnings(as.integer(gsub("\\D+", "", Factor)))
  ) |>
  dplyr::arrange(Factor_num) |>
  dplyr::mutate(Factor = forcats::fct_inorder(Factor))

## Plot (blue tones, like your previous figure)
p_rna_gsea_counts <- ggplot(
  rna_gsea_counts_plot,
  aes(x = Factor, y = sig_terms, fill = sig_terms)
) +
  geom_col(width = 0.8) +
  coord_flip() +
  labs(
    title = "MOFA2 RNA factors (GSEA, FDR \u2264 0.05): significant GO:BP terms",
    x = NULL,
    y = "Number of enriched GO:BP terms",
    fill = "GO:BP term count"
  ) +
  scale_fill_gradient(low = "#cfe1f2", high = "#0b3c5d") +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.margin = margin(t = 15, r = 40, b = 15, l = 15),
    legend.position = "right",
    legend.title = element_text(size = 11),
    legend.text  = element_text(size = 10)
  )

print(p_rna_gsea_counts)

## Save (prefer PLOTS_DIR if defined; otherwise save into RNA_GSEA_DIR)
out_dir_plot <- if (exists("PLOTS_DIR")) PLOTS_DIR else RNA_GSEA_DIR
out_png <- file.path(out_dir_plot, "MOFA_RNA_GSEA_sig_term_counts.png")

ggsave(
  filename = out_png,
  plot     = p_rna_gsea_counts,
  width    = 9,
  height   = 4.8,
  dpi      = 300,
  bg       = "white"
)

message("Saved plot: ", out_png)

## ---------------------------
## 9.3) ORA (METH view) — GO:BP
## ---------------------------

TOPN_ORA <- 200
pCut <- 0.05
qCut <- 0.05

# Universe = all features in methylation weights mapped to ENTREZ
meth_universe_entrez <- sym2ent(rownames(W_METH_topN)) %>%
  na.omit() %>%
  unique()

message("MOFA METH ORA universe size: ", length(meth_universe_entrez), " ENTREZ IDs")

meth_ora_summary <- list()

for (k in seq_len(ncol(W_METH_topN))) {
  
  w <- W_METH_topN[, k]
  names(w) <- rownames(W_METH_topN)
  
  ord <- order(abs(w), decreasing = TRUE)
  sel <- names(w)[ord][seq_len(min(TOPN_ORA, length(w)))]
  
  ent <- sym2ent(sel) %>%
    na.omit() %>%
    unique()
  
  if (length(ent) < 10L) {
    message("MOFA METH ORA ", colnames(W_METH_topN)[k], ": <10 ENTREZ IDs, skipped.")
    next
  }
  
  eg <- tryCatch(
    clusterProfiler::enrichGO(
      gene          = ent,
      universe      = meth_universe_entrez,
      OrgDb         = org.Hs.eg.db,
      keyType       = "ENTREZID",
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = pCut,
      qvalueCutoff  = qCut,
      readable      = TRUE
    ),
    error = function(e) NULL
  )
  
  if (!is.null(eg) && nrow(as.data.frame(eg)) > 0) {
    df <- as.data.frame(eg)
    df$Factor <- colnames(W_METH_topN)[k]
    
    # Save per-factor table
    readr::write_csv(
      df,
      file.path(METH_ORA_DIR, sprintf("MOFA_METH_ORA_%s.csv", colnames(W_METH_topN)[k]))
    )
    
    meth_ora_summary[[colnames(W_METH_topN)[k]]] <- df
  } else {
    message("MOFA METH ORA ", colnames(W_METH_topN)[k], ": no significant terms at FDR <= 0.05.")
  }
}

# Combine all METH ORA results
if (length(meth_ora_summary) > 0) {
  
  ## 9A) Combine all factor ORA tables (preserve Factor from list names)
  meth_ora_all <- dplyr::bind_rows(meth_ora_summary, .id = "Factor")
  
  # Save ALL ORA results
  readr::write_csv(
    meth_ora_all,
    file.path(METH_ORA_DIR, "MOFA_METH_ORA_ALL.csv")
  )
  
  ## 9B) Filter significant terms at FDR <= qCut
  meth_ora_sig <- meth_ora_all |>
    dplyr::mutate(p.adjust = suppressWarnings(as.numeric(p.adjust))) |>
    dplyr::filter(!is.na(p.adjust), p.adjust <= qCut)
  
  readr::write_csv(
    meth_ora_sig,
    file.path(METH_ORA_DIR, sprintf("MOFA_METH_ORA_FDR%.2f.csv", qCut))
  )
  
  ## 9C) Count significant terms per factor
  meth_ora_counts <- meth_ora_sig |>
    dplyr::count(Factor, name = "sig_terms") |>
    dplyr::arrange(dplyr::desc(sig_terms))
  
  readr::write_csv(
    meth_ora_counts,
    file.path(METH_ORA_DIR, "MOFA_METH_ORA_sig_counts.csv")
  )
  
  message("METH ORA done Factors with significant terms: ", nrow(meth_ora_counts))
  
} else {
  message("METH ORA done  No factors produced results.")
}

message("STEP 9 DONE Enrichment outputs saved under OUT_DIR")

# plot:==========

# Safety: only plot if we have counts
if (exists("meth_ora_counts") && is.data.frame(meth_ora_counts) && nrow(meth_ora_counts) > 0) {
  
  # Order factors by count (ascending) for a clean plot
  meth_ora_counts <- meth_ora_counts |>
    dplyr::arrange(sig_terms) |>
    dplyr::mutate(Factor = factor(Factor, levels = Factor))
  
  p_meth_ora <- ggplot2::ggplot(meth_ora_counts, ggplot2::aes(x = sig_terms, y = Factor, fill = sig_terms)) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::scale_fill_gradient(low = "#cfe1f2", high = "#0b3c5d") +
    ggplot2::labs(
      title = sprintf("MOFA2 methylation factors (ORA, FDR \u2264 %.2f): significant GO:BP terms", qCut),
      x = "Number of enriched GO:BP terms",
      y = NULL,
      fill = "GO:BP term count"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "right",
      plot.margin = grid::unit(c(10, 18, 10, 12), "pt")
    )
  
  # Show in Plots panel
  print(p_meth_ora)
  
  # Save
  out_png <- file.path(PLOTS_DIR, sprintf("MOFA_METH_ORA_sig_counts_FDR%.2f.png", qCut))
  ggplot2::ggsave(
    filename = out_png,
    plot = p_meth_ora,
    width = 10,
    height = 4.8,
    units = "in",
    dpi = 300,
    bg = "white"
  )
  
  message("Saved ORA plot: ", out_png)
  
} else {
  message("No meth_ora_counts available to plot (no significant ORA terms).")
}

#============== THE END ===========================================================
