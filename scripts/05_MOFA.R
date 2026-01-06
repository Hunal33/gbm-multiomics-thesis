
## ============================================================
## 05_MOFA.R — MOFA2 main run
## Purpose:
##   - Load global setup
##   - Read MOFA input data
##   - Prepare data list for model fitting
## ============================================================

## ---------------------------
## 0) LOAD SETUP (MANDATORY)
## ---------------------------
stopifnot(file.exists("scripts/00_setup_paths.R"))
source("scripts/00_setup_paths.R")

## ---------------------------
## 1) LOAD LIBRARIES
## ---------------------------
suppressPackageStartupMessages({
  library(MOFA2)
})

## ---------------------------
## 2) READ INPUT DATA
## ---------------------------
rna_mat   <- readRDS(RNA_INPUT_FILE)
meth_mat  <- readRDS(METH_INPUT_FILE)
sample_md <- readRDS(SAMPLE_META_FILE)

## ---------------------------
## 3) BASIC SANITY CHECKS
## ---------------------------
stopifnot(
  is.matrix(rna_mat) || is.data.frame(rna_mat),
  is.matrix(meth_mat) || is.data.frame(meth_mat),
  identical(colnames(rna_mat), colnames(meth_mat))
)

message("MOFA inputs loaded.")
message("RNA:  ", paste(dim(rna_mat), collapse = " x "))
message("METH: ", paste(dim(meth_mat), collapse = " x "))


## ---------------------------
## 4) PREPARE MOFA OBJECT
## ---------------------------
rna_mat  <- as.matrix(rna_mat)
meth_mat <- as.matrix(meth_mat)

stopifnot(
  all(is.finite(rna_mat)),
  all(is.finite(meth_mat))
)

data_list <- list(
  RNA  = rna_mat,
  METH = meth_mat
)

mofa_obj <- create_mofa(data_list)

## ---------------------------
## 5) DATA OPTIONS
## ---------------------------
data_opts <- get_default_data_options(mofa_obj)
data_opts$scale_views <- TRUE

## ---------------------------
## 6) MODEL OPTIONS
## ---------------------------
model_opts <- get_default_model_options(mofa_obj)
model_opts$num_factors <- 15

## ---------------------------
## 7) TRAINING OPTIONS
## ---------------------------
train_opts <- get_default_training_options(mofa_obj)
train_opts$seed <- 42
train_opts$maxiter <- 1500
train_opts$convergence_mode <- "medium"
train_opts$drop_factor_threshold <- 0.01

## ---------------------------
## 8) PREPARE + RUN
## ---------------------------
mofa_obj <- prepare_mofa(
  object = mofa_obj,
  data_options     = data_opts,
  model_options    = model_opts,
  training_options = train_opts
)

mofa_fit <- run_mofa(mofa_obj, use_basilisk = TRUE)

## ---------------------------
## 9) SAVE OUTPUTS (RDS folder)
## ---------------------------
saveRDS(mofa_fit, file.path(mofa_dir, "rds", "mofa_fit.rds"))

message("MOFA fit saved: ", file.path(mofa_dir, "rds", "mofa_fit.rds"))

## ---------------------------
## 10) EXTRACT + SAVE CORE OUTPUTS (RDS)
## ---------------------------
# Z: factors (samples x factors)
Z <- get_factors(mofa_fit, factors = "all")$group1
# W: loadings (features x factors) per view
W_list <- get_weights(mofa_fit, factors = "all")

# Save
saveRDS(Z,      file.path(mofa_dir, "rds", "Z_all.rds"))
saveRDS(W_list, file.path(mofa_dir, "rds", "W_all.rds"))

message("Saved: ", file.path(mofa_dir, "rds", "Z_all.rds"))
message("Saved: ", file.path(mofa_dir, "rds", "W_all.rds"))

## ---------------------------
## 11) FACTOR CLASS + PLOT (by-view) 
##     (Shows METH-specific / RNA-specific / Shared / Weak)
## ---------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# Load ranked VE table (Step 12 output)
ve_df <- readRDS(file.path(mofa_dir, "rds", "VE_per_factor_ranked.rds"))

## ---- 11A) FACTOR CLASSIFICATION
weak_thr <- 2.0    # <2% total VE
rna_thr  <- 0.80   # RNA share >= 80%
meth_thr <- 0.20   # RNA share <= 20%

ve_df <- ve_df %>%
  mutate(
    RNA_share = VE_RNA / pmax(VE_TOTAL, 1e-9),
    Class = case_when(
      VE_TOTAL < weak_thr        ~ "Weak",
      RNA_share >= rna_thr       ~ "RNA-specific",
      RNA_share <= meth_thr      ~ "METH-specific",
      TRUE                       ~ "Shared"
    )
  )

# Save + quick summary
saveRDS(ve_df, file.path(mofa_dir, "rds", "MOFA_factor_classes.rds"))
write.csv(ve_df, file.path(mofa_dir, "rds", "MOFA_factor_classes.csv"), row.names = FALSE)

print(table(ve_df$Class))

## ---- 11B) SUM VE BY CLASS (stacked by view)
ve_class_sum <- ve_df %>%
  group_by(Class) %>%
  summarise(
    RNA  = sum(VE_RNA,  na.rm = TRUE),
    METH = sum(VE_METH, na.rm = TRUE),
    n_factors = dplyr::n(),
    .groups = "drop"
  ) %>%
  mutate(Class = factor(Class, levels = c("METH-specific", "RNA-specific", "Shared", "Weak"))) %>%
  pivot_longer(cols = c(RNA, METH), names_to = "View", values_to = "Percent") %>%
  mutate(
    View = factor(View, levels = c("METH", "RNA")) # bottom->top with reverse stack
  )

# counts for subtitle
counts_txt <- ve_df %>%
  count(Class) %>%
  mutate(Class = factor(Class, levels = c("METH-specific", "RNA-specific", "Shared", "Weak"))) %>%
  arrange(Class)

subtitle_txt <- paste0(
  paste0(as.character(counts_txt$Class), "=", counts_txt$n),
  collapse = " | "
)

p_class <- ggplot(ve_class_sum, aes(x = Class, y = Percent, fill = View)) +
  geom_col(position = position_stack(reverse = TRUE), width = 0.75) +
  scale_fill_manual(
    limits = c("METH", "RNA"),
    values = c("RNA" = "#0B3C5D", "METH" = "#CFE1F2"),
    name = "View"
  ) +
  labs(
    title = "MOFA2: Explained variance by factor class (by-view)",
    subtitle = subtitle_txt,
    y = "Total % variance explained",
    x = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "right"
  )

print(p_class)

out_png <- file.path(mofa_dir, "plots", "MOFA_VE_by_factor_class.png")
ggsave(out_png, p_class, width = 7, height = 4, dpi = 300, bg = "white")

stopifnot(file.exists(out_png))

## ---------------------------
## 13) SAVE TOP FACTORS (Z/W)  
## ---------------------------

# Load classified VE table (must contain Rank, VE_RNA, VE_METH, VE_TOTAL, Class)
ve_df <- readRDS(file.path(mofa_dir, "rds", "MOFA_factor_classes.rds"))

# Load Z/W (already saved in Step 10)
Z_all <- readRDS(file.path(mofa_dir, "rds", "Z_all.rds"))
W_all <- readRDS(file.path(mofa_dir, "rds", "W_all.rds"))

# Top K by Rank (SAFE: does not depend on colnames)
K <- 15
top_idx <- ve_df$Rank[ve_df$Rank <= K]

Z_top <- Z_all[, top_idx, drop = FALSE]
W_top <- lapply(W_all, function(Wv) Wv[, top_idx, drop = FALSE])

# Rename columns to F1..FK for consistency
colnames(Z_top) <- paste0("F", seq_len(ncol(Z_top)))
W_top <- lapply(W_top, function(Wv) { colnames(Wv) <- paste0("F", seq_len(ncol(Wv))); Wv })

# Save
saveRDS(Z_top, file.path(mofa_dir, "rds", "Z_top15.rds"))
saveRDS(W_top, file.path(mofa_dir, "rds", "W_top15.rds"))
write.csv(ve_df, file.path(mofa_dir, "tables", "MOFA_factor_classes.csv"), row.names = FALSE)

message("Saved: Z_top15.rds, W_top15.rds, MOFA_factor_classes.csv")

## ---------------------------
## 14) PLOT — % variance explained per factor (by view)  (ranked)
## ---------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(forcats)
})

ve_df <- readRDS(file.path(mofa_dir, "rds", "MOFA_factor_classes.rds"))

ve_long <- ve_df %>%
  select(Factor, VE_RNA, VE_METH) %>%
  pivot_longer(
    cols = c(VE_RNA, VE_METH),
    names_to = "View",
    values_to = "Percent"
  ) %>%
  mutate(
    View = recode(View, VE_RNA = "RNA", VE_METH = "METH"),
    View = factor(View, levels = c("RNA", "METH")),
    Factor = fct_relevel(Factor, ve_df$Factor)
  )

p_ve <- ggplot(ve_long, aes(x = Factor, y = Percent, fill = View)) +
  geom_col(position = position_stack(reverse = TRUE), width = 0.75) +
  scale_fill_manual(
    limits = c("RNA", "METH"),
    values = c("RNA" = "#0B3C5D", "METH" = "#CFE1F2"),
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

out_png <- file.path(mofa_dir, "plots", "MOFA_VE_stacked_per_factor.png")
ggsave(out_png, p_ve, width = 8, height = 4, dpi = 300, bg = "white")
stopifnot(file.exists(out_png))

message("Saved plot: MOFA_VE_stacked_per_factor.png")

## ---------------------------
## 15) MOFA RNA — PREPARE RANKED LISTS FOR GSEA (per factor)
## ---------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

# Load W (RNA) from saved object
W_all <- readRDS(file.path(mofa_dir, "rds", "W_all.rds"))
W_RNA <- W_all[["RNA"]]

stopifnot(!is.null(W_RNA))

# Output folder
gsea_rna_dir <- file.path(mofa_dir, "GSEA", "tables")
dir.create(gsea_rna_dir, recursive = TRUE, showWarnings = FALSE)

# Ensure numeric matrix
W_RNA <- as.matrix(W_RNA)
storage.mode(W_RNA) <- "double"

# For each factor: create ranked vector (named) and save as CSV
for (j in seq_len(ncol(W_RNA))) {
  scores <- W_RNA[, j]
  names(scores) <- rownames(W_RNA)
  
  scores <- sort(scores, decreasing = TRUE)
  
  out_df <- data.frame(
    gene  = names(scores),
    score = as.numeric(scores),
    stringsAsFactors = FALSE
  )
  
  out_file <- file.path(gsea_rna_dir, sprintf("MOFA_RNA_ranked_F%d.csv", j))
  write_csv(out_df, out_file)
}

message("Saved RNA ranked lists for GSEA in: ", gsea_rna_dir)

## ---------------------------
## 16) MOFA RNA — GSEA (GO:BP) across ALL factors (robust mapping)
## Inputs:
##   - results/MOFA/rds/W_all.rds  (expects W_all[["RNA"]] with gene symbols as rownames)
## Outputs:
##   - results/MOFA/GSEA/tables/MOFA_RNA_GSEA_F#.csv
##   - results/MOFA/tables/MOFA_RNA_GSEA_sig_counts.csv
##   - results/MOFA/plots/MOFA_RNA_GSEA_sig_terms_horizontal.png
## ---------------------------

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(dplyr)
  library(readr)
  library(ggplot2)
})

stopifnot(dir.exists(mofa_dir))

GSEA_RNA_DIR <- file.path(mofa_dir, "GSEA", "tables")
dir.create(GSEA_RNA_DIR, recursive = TRUE, showWarnings = FALSE)

# For summary + plot (NOT under GSEA)
GSEA_TABLES_DIR <- file.path(mofa_dir, "GSEA", "tables")
GSEA_PLOTS_DIR  <- file.path(mofa_dir, "GSEA", "plots")

dir.create(GSEA_TABLES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(GSEA_PLOTS_DIR,  recursive = TRUE, showWarnings = FALSE)

## ---------------------------
## Load RNA weights (features x K)
## ---------------------------
W_all <- readRDS(file.path(mofa_dir, "rds", "W_all.rds"))
W_RNA <- as.matrix(W_all[["RNA"]])
storage.mode(W_RNA) <- "double"

stopifnot(!is.null(colnames(W_RNA)), !is.null(rownames(W_RNA)))
## ---- CLEAN RNA gene symbols (remove _RNA suffix)
rownames(W_RNA) <- gsub("_RNA$", "", rownames(W_RNA))

K <- ncol(W_RNA)
message("Loaded MOFA W_RNA: ", paste(dim(W_RNA), collapse = " x "), " (K=", K, ")")

## ---- helpers: robust SYMBOL -> ENTREZ
sym2ent <- function(sym) {
  sym <- unique(na.omit(as.character(sym)))
  if (length(sym) == 0L) return(character(0))
  
  sym_up <- toupper(sym)
  
  # key validation (prevents "None of the keys..." error)
  all_keys <- AnnotationDbi::keys(org.Hs.eg.db, keytype = "SYMBOL")
  sym_ok   <- intersect(sym_up, all_keys)
  if (length(sym_ok) == 0L) {
    out <- rep(NA_character_, length(sym_up))
    names(out) <- sym
    return(out)
  }
  
  res <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys      = sym_ok,
    keytype   = "SYMBOL",
    column    = "ENTREZID",
    multiVals = "first"
  )
  names(res) <- toupper(names(res))
  
  out <- res[match(sym_up, names(res))]
  names(out) <- sym
  out
}

prep_gsea_rank <- function(w_named) {
  stopifnot(!is.null(names(w_named)))
  ent <- sym2ent(names(w_named))
  keep <- !is.na(ent) & !duplicated(ent)
  if (!any(keep)) return(numeric(0))
  r <- setNames(as.numeric(w_named[keep]), ent[keep])
  r[order(r, decreasing = TRUE)]
}

## ---- run GSEA per factor
gsea_results_list <- list()

for (f in colnames(W_RNA)) {
  
  w <- W_RNA[, f]
  names(w) <- rownames(W_RNA)
  
  ranks <- prep_gsea_rank(w)
  message("=== MOFA RNA ", f, " | ranked genes: ", length(ranks), " ===")
  
  if (length(ranks) < 50) {
    message("Skip ", f, ": <50 ranked genes after ENTREZ mapping.")
    next
  }
  
  gse <- tryCatch(
    gseGO(
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
    error = function(e) {
      message("ERROR in ", f, ": ", e$message)
      NULL
    }
  )
  
  if (!is.null(gse) && nrow(as.data.frame(gse)) > 0) {
    df <- as.data.frame(gse)
    df$Factor <- f
    
    out_csv <- file.path(GSEA_RNA_DIR, paste0("MOFA_RNA_GSEA_", f, ".csv"))
    readr::write_csv(df, out_csv)
    
    gsea_results_list[[f]] <- df
    message("Saved: ", out_csv)
  } else {
    message("No significant GO:BP terms for ", f, " at FDR <= 0.05.")
  }
}


## ---- summary counts (FDR<=0.05)
if (length(gsea_results_list) > 0) {
  
  gsea_all <- dplyr::bind_rows(gsea_results_list, .id = "FactorID") %>%
    mutate(p.adjust = suppressWarnings(as.numeric(p.adjust))) %>%
    filter(!is.na(p.adjust), p.adjust <= 0.05)
  
  gsea_counts <- gsea_all %>%
    group_by(Factor = FactorID) %>%
    summarise(sig_terms = n(), .groups = "drop") %>%
    arrange(desc(sig_terms))
  
  out_counts <- file.path(TABLES_DIR, "MOFA_RNA_GSEA_sig_counts.csv")
  readr::write_csv(gsea_counts, out_counts)
  message("Saved summary counts: ", out_counts)
  
  # Plot (MOFA-style)
  gsea_counts <- gsea_counts %>%
    mutate(Factor = factor(Factor, levels = rev(Factor)))
  
  p <- ggplot(gsea_counts, aes(x = sig_terms, y = Factor, fill = sig_terms)) +
    geom_col(width = 0.8) +
    scale_fill_gradient(low = "#CFE1F2", high = "#0B3C5D", name = "GO:BP term count") +
    labs(
      title = "MOFA RNA factors (GSEA, FDR ≤ 0.05): significant GO:BP terms",
      x = "Number of enriched GO:BP terms",
      y = NULL
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )
  
  print(p)
  
  out_png <- file.path(GSEA_PLOTS_DIR, "MOFA_RNA_GSEA_sig_terms_horizontal.png")
  ggsave(out_png, p, width = 10, height = 5, dpi = 300, bg = "white")
  message("Saved: ", out_png)
  
} else {
  message("No MOFA RNA GSEA results were produced.")
}

message("STEP 16 DONE")


## ---------------------------
## 17) MOFA METH — PREPARE GENE LISTS FOR ORA (per factor)
## ---------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

# Load W (METH)
W_all  <- readRDS(file.path(mofa_dir, "rds", "W_all.rds"))
W_METH <- as.matrix(W_all[["METH"]])
storage.mode(W_METH) <- "double"

stopifnot(!is.null(colnames(W_METH)), !is.null(rownames(W_METH)))

# Output dir
ora_gene_dir <- file.path(mofa_dir, "ORA", "tables")
dir.create(ora_gene_dir, recursive = TRUE, showWarnings = FALSE)

# Parameters
TOP_N <- 300  

for (j in seq_len(ncol(W_METH))) {
  
  w <- W_METH[, j]
  names(w) <- rownames(W_METH)
  
  # Rank by absolute loading
  ord <- order(abs(w), decreasing = TRUE)
  top_genes <- names(w)[ord][seq_len(min(TOP_N, length(w)))]
  
  out_df <- data.frame(
    gene = top_genes,
    stringsAsFactors = FALSE
  )
  
  out_file <- file.path(ora_gene_dir, sprintf("MOFA_METH_ORA_genes_F%d.csv", j))
  write_csv(out_df, out_file)
}

message("MOFA METH ORA gene lists prepared.")

## =========================================================
## 18_MOFA_ORA_METH_GO_BP.R — MOFA METH ORA (GO:BP) across ALL active factors
## Purpose:
##   - Run ORA on METH factor weights using top-|weights| genes
##   - No manual factor selection: uses "active" definition
## Input:
##   - results/MOFA/rds/W_all.rds  (expects W_all[["METH"]] genes x K; rownames are gene symbols)
## Output:
##   - results/MOFA/ORA/tables/MOFA_METH_ORA_<FACTOR>.csv
##   - results/MOFA/tables/MOFA_METH_ORA_sig_counts.csv
##   - results/MOFA/plots/MOFA_METH_ORA_sig_terms_horizontal.png
## =========================================================

source("scripts/00_setup_paths.R")
stopifnot(exists("mofa_dir"), dir.exists(mofa_dir))

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(forcats)
})

## ---------------------------
## Output dirs (match your structure)
## ---------------------------
## ---------------------------
## Output dirs (MOFA ORA)
## ---------------------------
ORA_METH_DIR   <- file.path(mofa_dir, "ORA", "tables")
ORA_PLOTS_DIR  <- file.path(mofa_dir, "ORA", "plots")

dir.create(ORA_METH_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(ORA_PLOTS_DIR, recursive = TRUE, showWarnings = FALSE)

## ---------------------------
## Load METH weights (genes x K)
## ---------------------------
W_all  <- readRDS(file.path(mofa_dir, "rds", "W_all.rds"))
W_METH <- as.matrix(W_all[["METH"]])
storage.mode(W_METH) <- "double"

stopifnot(!is.null(colnames(W_METH)))
stopifnot(!is.null(rownames(W_METH)))  # should be gene symbols (may contain suffix)

# Clean gene symbols if needed
rownames(W_METH) <- toupper(rownames(W_METH))
rownames(W_METH) <- gsub("_METH$", "", rownames(W_METH))
rownames(W_METH) <- gsub("_RNA$",  "", rownames(W_METH))
rownames(W_METH) <- gsub("\\|.*$", "", rownames(W_METH))

K <- ncol(W_METH)
message("Loaded MOFA W_METH: ", paste(dim(W_METH), collapse = " x "), " (K=", K, ")")

## ---------------------------
## Active factor detection
## ---------------------------
eps <- 1e-8
meth_sd  <- apply(W_METH, 2, sd)
meth_max <- apply(abs(W_METH), 2, max, na.rm = TRUE)

active_idx     <- which(meth_sd > 0 & meth_max > eps)
active_factors <- colnames(W_METH)[active_idx]

message("Active METH factors: ", length(active_factors), " / ", K)
print(active_factors)

## ---------------------------
## Helpers: SYMBOL -> ENTREZ, build ORA gene set
## ---------------------------
sym2ent <- function(sym) {
  sym <- unique(na.omit(as.character(sym)))
  if (length(sym) == 0L) return(character(0))
  
  sym_up   <- toupper(sym)
  all_keys <- AnnotationDbi::keys(org.Hs.eg.db, keytype = "SYMBOL")
  sym_ok   <- intersect(sym_up, all_keys)
  
  if (length(sym_ok) == 0L) return(character(0))
  
  res <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys      = sym_ok,
    keytype   = "SYMBOL",
    column    = "ENTREZID",
    multiVals = "first"
  )
  unique(na.omit(as.character(res)))
}

prep_ora_genes <- function(w_named, topN = 200) {
  stopifnot(!is.null(names(w_named)))
  ord <- order(abs(w_named), decreasing = TRUE)
  top_sym <- names(w_named)[ord][seq_len(min(topN, length(ord)))]
  sym2ent(top_sym)
}

## ---------------------------
## Run ORA per active factor
## ---------------------------
TOPN_ORA <- 300
ora_results_list <- list()

for (f in active_factors) {
  
  w <- W_METH[, f]
  names(w) <- rownames(W_METH)
  
  genes_entrez <- prep_ora_genes(w, topN = TOPN_ORA)
  
  message("=== ORA METH ", f, " | topN=", TOPN_ORA, " | ENTREZ genes: ", length(genes_entrez), " ===")
  
  if (length(genes_entrez) < 20) {
    message("Skip ", f, ": <20 ENTREZ genes after mapping.")
    next
  }
  
  enr <- tryCatch(
    enrichGO(
      gene          = genes_entrez,
      OrgDb         = org.Hs.eg.db,
      keyType       = "ENTREZID",
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.05,
      readable      = TRUE
    ),
    error = function(e) {
      message("ERROR in ", f, ": ", e$message)
      NULL
    }
  )
  
  if (!is.null(enr) && nrow(as.data.frame(enr)) > 0) {
    df <- as.data.frame(enr)
    df$Factor <- f
    
    out_csv <- file.path(ORA_METH_DIR, paste0("MOFA_METH_ORA_", f, ".csv"))
    readr::write_csv(df, out_csv)
    
    ora_results_list[[f]] <- df
    message("Saved: ", out_csv)
  } else {
    message("No significant GO:BP terms for ", f, " at FDR <= 0.05.")
  }
}

## ---------------------------
## Summary: significant term counts per factor (FDR <= 0.05)
## ---------------------------
if (length(ora_results_list) > 0) {
  
  ora_all <- dplyr::bind_rows(ora_results_list, .id = "FactorID") %>%
    dplyr::mutate(p.adjust = suppressWarnings(as.numeric(p.adjust))) %>%
    dplyr::filter(!is.na(p.adjust), p.adjust <= 0.05)
  
  ora_counts <- ora_all %>%
    dplyr::group_by(Factor = FactorID) %>%
    dplyr::summarise(sig_terms = dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(sig_terms))
  
  out_counts <- file.path(TABLES_DIR, "MOFA_METH_ORA_sig_counts.csv")
  readr::write_csv(ora_counts, out_counts)
  message("Saved summary counts: ", out_counts)
  
  # Plot (MOFA-style)
  ora_counts <- ora_counts %>%
    mutate(Factor = factor(Factor, levels = rev(Factor)))
  
  p <- ggplot(ora_counts, aes(x = sig_terms, y = Factor, fill = sig_terms)) +
    geom_col(width = 0.8) +
    scale_fill_gradient(low = "#cfe1f2", high = "#0b3c5d", name = "GO:BP term count") +
    labs(
      title = "MOFA methylation factors (ORA, FDR ≤ 0.05): significant GO:BP terms",
      x = "Number of enriched GO:BP terms",
      y = NULL
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )
  
  print(p)
  
  out_png <- file.path(ORA_PLOTS_DIR, "MOFA_METH_ORA_sig_terms_horizontal.png")
  ggsave(out_png, p, width = 10, height = 5, dpi = 300, bg = "white")
  message("Saved ORA plot: ", out_png)
  
} else {
  message("No METH ORA results were produced (no active factors or no significant terms).")
}

message("MOFA METH ORA (GO:BP) DONE.")
