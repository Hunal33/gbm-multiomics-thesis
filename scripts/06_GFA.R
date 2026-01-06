## ============================================================
## 06_GFA.R — Prepare GFA inputs (GitHub-safe)
## Purpose:
##   - Load paired RNA + METH inputs (already MOFA-compatible)
##   - Feature-wise scaling per view (row-wise z-score)
##   - Save GFA input object (NO model fitting yet)
## Outputs:
##   - results/GFA/rds/gfa_input_list.rds
##   - results/GFA/rds/gfa_input_scaled_{RNA,METH}.rds  (optional, helpful)
## ============================================================

stopifnot(file.exists("scripts/00_setup_paths.R"))
source("scripts/00_setup_paths.R")

suppressPackageStartupMessages({
  library(Matrix)
})

## ---------------------------
## 1) LOAD INPUTS
## ---------------------------
rna_mat     <- readRDS(RNA_INPUT_FILE)
meth_mat    <- readRDS(METH_INPUT_FILE)
sample_meta <- readRDS(SAMPLE_META_FILE)

rna_mat  <- as.matrix(rna_mat)
meth_mat <- as.matrix(meth_mat)

storage.mode(rna_mat)  <- "double"
storage.mode(meth_mat) <- "double"

## ---------------------------
## 2) SAFETY CHECKS (paired, aligned)
## ---------------------------
stopifnot(
  !is.null(colnames(rna_mat)),
  !is.null(colnames(meth_mat)),
  identical(colnames(rna_mat), colnames(meth_mat)),
  ncol(rna_mat) == nrow(sample_meta)
)

message("Loaded RNA:  ", nrow(rna_mat),  " features x ", ncol(rna_mat),  " samples")
message("Loaded METH: ", nrow(meth_mat), " features x ", ncol(meth_mat), " samples")
message("Samples (meta): ", nrow(sample_meta))

## ---------------------------
## 3) BUILD GFA MATRICES (CORRECT)
## ---------------------------
# Transpose ONCE: samples x features
X_rna  <- t(rna_mat)
X_meth <- t(meth_mat)

## Scale per feature (GFA-style)
X_rna  <- scale(X_rna,  center = TRUE, scale = TRUE)
X_meth <- scale(X_meth, center = TRUE, scale = TRUE)

# Replace NA from zero-variance features
X_rna[is.na(X_rna)]   <- 0
X_meth[is.na(X_meth)] <- 0

stopifnot(nrow(X_rna) == nrow(X_meth))

message("GFA input matrices:")
message("  RNA :  ", paste(dim(X_rna),  collapse = " x "))
message("  METH:  ", paste(dim(X_meth), collapse = " x "))

## ---------------------------
## 4) BUILD INPUT LIST
## ---------------------------
gfa_input <- list(
  RNA  = X_rna,
  METH = X_meth
)

## ---------------------------
## 5) ENSURE OUTPUT FOLDERS
## ---------------------------
gfa_subdirs <- c(
  file.path(gfa_dir, "rds"),
  file.path(gfa_dir, "tables"),
  file.path(gfa_dir, "plots"),
  file.path(gfa_dir, "GSEA", "tables"),
  file.path(gfa_dir, "GSEA", "plots"),
  file.path(gfa_dir, "ORA",  "tables"),
  file.path(gfa_dir, "ORA",  "plots")
)
for (d in gfa_subdirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

## ---------------------------
## 6) SAVE (CANONICAL)
## ---------------------------
out_rds_dir <- file.path(gfa_dir, "rds")

saveRDS(gfa_input, file.path(out_rds_dir, "gfa_input_list.rds"))
saveRDS(X_rna,     file.path(out_rds_dir, "gfa_input_scaled_RNA.rds"))
saveRDS(X_meth,    file.path(out_rds_dir, "gfa_input_scaled_METH.rds"))

message("Saved:")
message(" - gfa_input_list.rds")
message(" - gfa_input_scaled_RNA.rds")
message(" - gfa_input_scaled_METH.rds")

X_list <- readRDS(file.path(gfa_dir,"rds","gfa_input_list.rds"))
dim(X_list$RNA)

## ============================================================
## 02_GFA_fit_model.R — Fit GFA model (opts required)
## Purpose:
##   - Load prepared GFA inputs (gfa_input_list.rds)
##   - Fit GFA with fixed K (match MOFA: 15)
##   - Save fitted model object for reproducibility
## Outputs:
##   results/GFA/rds/
##     - gfa_fit.rds
##     - gfa_runtime_seconds.txt
##     - gfa_fit_sessionInfo.txt
## ============================================================

suppressPackageStartupMessages({
  library(GFA)
})

## ---------------------------
## 0) PARAMETERS
## ---------------------------
K    <- 15
SEED <- 42

## ---------------------------
## 1) LOAD INPUT (samples x features per view)
## ---------------------------
in_file <- file.path(gfa_dir, "rds", "gfa_input_list.rds")
stopifnot(file.exists(in_file))

X_list <- readRDS(in_file)
stopifnot(is.list(X_list), all(c("RNA","METH") %in% names(X_list)))

X_list$RNA  <- as.matrix(X_list$RNA)
X_list$METH <- as.matrix(X_list$METH)
storage.mode(X_list$RNA)  <- "double"
storage.mode(X_list$METH) <- "double"

stopifnot(nrow(X_list$RNA) == nrow(X_list$METH))

N  <- nrow(X_list$RNA)
Dr <- ncol(X_list$RNA)
Dm <- ncol(X_list$METH)

message("Loaded X_list (samples x features):")
message("  RNA :  ", N, " x ", Dr)
message("  METH:  ", N, " x ", Dm)

## ---------------------------
## 2) OPTIONS (required by this GFA)
## ---------------------------
set.seed(SEED)

opts <- GFA::getDefaultOpts()

# iterations (you can tune later)
opts$iter.max    <- 2000
opts$iter.burnin <- 800
opts$verbose     <- 2

# set K (some versions use k, some K; we set both)
opts$k <- K
opts$K <- K

# noise prior (recommended)
opts <- GFA::informativeNoisePrior(X_list, opts)

message("Fitting GFA with K = ", K, " ...")

## ---------------------------
## 3) FIT + RUNTIME
## ---------------------------
t0 <- Sys.time()

# IMPORTANT: your GFA requires opts; keep this exact call style
gfa_fit <- GFA::gfa(X_list, opts, K = K)

t1 <- Sys.time()
runtime_sec <- as.numeric(difftime(t1, t0, units = "secs"))

message("GFA fit done. Runtime (sec): ", round(runtime_sec, 1))

## ---------------------------
## 4) SAVE OUTPUTS
## ---------------------------
out_rds_dir <- file.path(gfa_dir, "rds")
dir.create(out_rds_dir, recursive = TRUE, showWarnings = FALSE)

saveRDS(gfa_fit, file.path(out_rds_dir, "gfa_fit.rds"))

writeLines(
  as.character(runtime_sec),
  con = file.path(out_rds_dir, "gfa_runtime_seconds.txt")
)

sink(file.path(out_rds_dir, "gfa_fit_sessionInfo.txt"))
cat("K=", K, "\n")
cat("SEED=", SEED, "\n")
cat("runtime_sec=", runtime_sec, "\n\n")
print(sessionInfo())
sink()

message("Saved: ", file.path(out_rds_dir, "gfa_fit.rds"))

## ============================================================
## 03_GFA_extract_ZW.R — Extract Z (scores) and W (loadings)
## Purpose:
##   - Load gfa_fit.rds
##   - Extract posterior mean Z (samples x K) from posterior$X
##   - Extract W (features_total x K) from posterior$w (2D or 3D)
##   - Split W into RNA / METH using D
##   - Restore feature names from saved scaled inputs
##   - Save:
##       results/GFA/rds/Z_all.rds
##       results/GFA/rds/W_all.rds  (list: RNA, METH)
## ============================================================

## ============================================================
## 03_GFA_extract_ZW.R — Extract Z and W (robust; reconstruct W if missing)
## Outputs:
##   results/GFA/rds/Z_all.rds   (samples x K)
##   results/GFA/rds/W_all.rds   (list: RNA, METH; features x K)
## ============================================================

stopifnot(file.exists("scripts/00_setup_paths.R"))
source("scripts/00_setup_paths.R")

suppressPackageStartupMessages({
  library(GFA)
})

## ---------------------------
## 1) LOAD FIT + INPUT LIST (needed if W must be reconstructed)
## ---------------------------
fit_file <- file.path(gfa_dir, "rds", "gfa_fit.rds")
stopifnot(file.exists(fit_file))
gfa_fit <- readRDS(fit_file)

in_file <- file.path(gfa_dir, "rds", "gfa_input_list.rds")
stopifnot(file.exists(in_file))
X_list <- readRDS(in_file)
stopifnot(is.list(X_list), all(c("RNA","METH") %in% names(X_list)))

X_list$RNA  <- as.matrix(X_list$RNA)   # samples x features
X_list$METH <- as.matrix(X_list$METH)
storage.mode(X_list$RNA)  <- "double"
storage.mode(X_list$METH) <- "double"

N  <- nrow(X_list$RNA)
Dr <- ncol(X_list$RNA)
Dm <- ncol(X_list$METH)
stopifnot(nrow(X_list$METH) == N)

## ---------------------------
## 2) Z = posterior mean of X  (iter x N x K)  =>  N x K
## ---------------------------
stopifnot(!is.null(gfa_fit$posterior), !is.null(gfa_fit$posterior$X))

X_obj <- gfa_fit$posterior$X
if (!is.null(dim(X_obj)) && length(dim(X_obj)) == 3) {
  Z <- apply(X_obj, c(2, 3), mean)
} else {
  Z <- as.matrix(X_obj)
}
Z <- as.matrix(Z)
storage.mode(Z) <- "double"

message("Z dim: ", paste(dim(Z), collapse = " x "))

## ---------------------------
## 3) W: try to read from model; if missing -> reconstruct from Y and Z
##    Goal: W_all is (D_total x K)
## ---------------------------
W_obj <- NULL

# best-case: posterior$w
if (!is.null(gfa_fit$posterior$w)) {
  W_obj <- gfa_fit$posterior$w
  message("W source: posterior$w")
}

# second: top-level w
if (is.null(W_obj) && !is.null(gfa_fit$w)) {
  W_obj <- gfa_fit$w
  message("W source: gfa_fit$w")
}

if (!is.null(W_obj)) {
  # posterior$w may be 3D (iter x D x K) or 2D (D x K)
  if (!is.null(dim(W_obj)) && length(dim(W_obj)) == 3) {
    W_all <- apply(W_obj, c(2, 3), mean)
  } else {
    W_all <- as.matrix(W_obj)
  }
  W_all <- as.matrix(W_all)
  storage.mode(W_all) <- "double"
  
} else {
  ## ---- reconstruct W from data (Y) and Z (ridge) ----
  message("W not stored in gfa_fit. Reconstructing W from Y and Z (ridge regression).")
  
  # Y: samples x D_total
  Y <- cbind(X_list$RNA, X_list$METH)
  stopifnot(nrow(Y) == nrow(Z))
  
  K <- ncol(Z)
  lambda <- 1e-6
  
  # Compute B = (Z'Z + λI)^{-1} Z'Y   =>  K x D_total
  A <- crossprod(Z) + diag(lambda, K)        # K x K
  B <- solve(A, crossprod(Z, Y))             # K x D_total
  
  W_all <- t(B)                              # D_total x K
  W_all <- as.matrix(W_all)
  storage.mode(W_all) <- "double"
}

message("W_all dim: ", paste(dim(W_all), collapse = " x "))

## sanity: must match D_total
stopifnot(nrow(W_all) == (Dr + Dm), nrow(Z) == N)

## ---------------------------
## 4) Split W_all into RNA / METH (features x K)
## ---------------------------
W_RNA  <- W_all[seq_len(Dr), , drop = FALSE]
W_METH <- W_all[(Dr + 1):(Dr + Dm), , drop = FALSE]

## restore names from X_list (samples x features => feature names in colnames)
if (!is.null(colnames(X_list$RNA)))  rownames(W_RNA)  <- colnames(X_list$RNA)
if (!is.null(colnames(X_list$METH))) rownames(W_METH) <- colnames(X_list$METH)

# Factor names
colnames(Z)      <- paste0("F", seq_len(ncol(Z)))
colnames(W_RNA)  <- paste0("F", seq_len(ncol(W_RNA)))
colnames(W_METH) <- paste0("F", seq_len(ncol(W_METH)))

# Sample names
if (!is.null(rownames(X_list$RNA))) rownames(Z) <- rownames(X_list$RNA)

message("W_RNA dim:  ", paste(dim(W_RNA),  collapse = " x "))
message("W_METH dim: ", paste(dim(W_METH), collapse = " x "))

## ---------------------------
## 5) Save outputs
## ---------------------------
out_dir <- file.path(gfa_dir, "rds")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

saveRDS(Z, file.path(out_dir, "Z_all.rds"))
saveRDS(list(RNA = W_RNA, METH = W_METH), file.path(out_dir, "W_all.rds"))

message("Saved:")
message(" - ", file.path(out_dir, "Z_all.rds"))
message(" - ", file.path(out_dir, "W_all.rds"))

## ============================================================
## 04_GFA_variance_explained_and_class.R
## Purpose:
##   - Load GFA Z/W (small-letter files: z_all.rds, w_all.rds)
##   - Compute variance explained (VE) per factor and per view
##   - Classify factors (RNA-specific / METH-specific / Shared / Weak)
##   - Save ranked VE table + class table
##   - Save top15 Z/W (renamed F1..F15)
##   - Plot VE stacked per factor + VE by class
## Outputs (all under results/GFA/):
##   rds/VE_per_factor_ranked.rds
##   tables/GFA_factor_classes.csv
##   rds/z_top15.rds, rds/w_top15.rds
##   plots/GFA_VE_stacked_per_factor.png
##   plots/GFA_VE_by_factor_class.png
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(forcats)
})

## ---------------------------
## 0) Ensure standard folders (GFA root)
## ---------------------------
dir.create(file.path(gfa_dir, "rds"),    recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(gfa_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(gfa_dir, "plots"),  recursive = TRUE, showWarnings = FALSE)

RDS_DIR    <- file.path(gfa_dir, "rds")
TABLES_DIR <- file.path(gfa_dir, "tables")
PLOTS_DIR  <- file.path(gfa_dir, "plots")

## ---------------------------
## 1) Load Z / W (small-letter!)
## ---------------------------
z_file <- file.path(RDS_DIR, "z_all.rds")
w_file <- file.path(RDS_DIR, "w_all.rds")

stopifnot(file.exists(z_file), file.exists(w_file))

Z <- readRDS(z_file)                 # N x K
W <- readRDS(w_file)                 # list(RNA, METH) each: features x K

stopifnot(is.matrix(Z) || is.data.frame(Z))
Z <- as.matrix(Z); storage.mode(Z) <- "double"

stopifnot(is.list(W), all(c("RNA","METH") %in% names(W)))
W_RNA  <- as.matrix(W[["RNA"]]);  storage.mode(W_RNA)  <- "double"
W_METH <- as.matrix(W[["METH"]]); storage.mode(W_METH) <- "double"

stopifnot(ncol(Z) == ncol(W_RNA), ncol(Z) == ncol(W_METH))
K <- ncol(Z)

message("Loaded Z:      ", paste(dim(Z), collapse=" x "))
message("Loaded W_RNA:  ", paste(dim(W_RNA), collapse=" x "))
message("Loaded W_METH: ", paste(dim(W_METH), collapse=" x "))

## ---------------------------
## 2) Load X_list to compute VE on centered data
##    (We use scaled inputs saved in rds)
## ---------------------------
X_rna  <- readRDS(file.path(RDS_DIR, "gfa_input_scaled_RNA.rds"))   # N x Dr
X_meth <- readRDS(file.path(RDS_DIR, "gfa_input_scaled_METH.rds"))  # N x Dm

X_rna  <- as.matrix(X_rna);  storage.mode(X_rna)  <- "double"
X_meth <- as.matrix(X_meth); storage.mode(X_meth) <- "double"

stopifnot(nrow(X_rna) == nrow(Z), nrow(X_meth) == nrow(Z))
stopifnot(ncol(X_rna) == nrow(W_RNA), ncol(X_meth) == nrow(W_METH))

## Center only (MOFA-like VE; no re-scaling here)
Xr_c <- scale(X_rna,  center = TRUE, scale = FALSE)
Xm_c <- scale(X_meth, center = TRUE, scale = FALSE)

tot_var_rna  <- sum(Xr_c^2)
tot_var_meth <- sum(Xm_c^2)

## ---------------------------
## 3) VE per factor and per view
##    For each factor k:
##      Xhat_rna_k  = Z[,k] %*% t(W_RNA[,k])
##      Xhat_meth_k = Z[,k] %*% t(W_METH[,k])
##      VE_view_k = sum(Xhat^2) / sum(X_centered^2)
## ---------------------------
ve_rna  <- numeric(K)
ve_meth <- numeric(K)

for (k in seq_len(K)) {
  zk <- Z[, k, drop = FALSE]  # N x 1
  
  Xr_hat <- zk %*% t(W_RNA[,  k, drop = FALSE])   # N x Dr
  Xm_hat <- zk %*% t(W_METH[, k, drop = FALSE])   # N x Dm
  
  ve_rna[k]  <- sum(Xr_hat^2) / tot_var_rna
  ve_meth[k] <- sum(Xm_hat^2) / tot_var_meth
}

ve_df <- tibble(
  Factor   = paste0("F", seq_len(K)),
  VE_RNA   = ve_rna  * 100,
  VE_METH  = ve_meth * 100
) %>%
  mutate(
    VE_TOTAL = VE_RNA + VE_METH
  ) %>%
  arrange(desc(VE_TOTAL)) %>%
  mutate(Rank = row_number())

saveRDS(ve_df, file.path(RDS_DIR, "VE_per_factor_ranked.rds"))
readr::write_csv(ve_df, file.path(TABLES_DIR, "GFA_VE_per_factor_ranked.csv"))
message("Saved VE table.")

## ---------------------------
## 4) Factor class (MOFA-parallel)
##    - Weak: total VE < 2%
##    - RNA-specific: RNA_share >= 0.80
##    - METH-specific: RNA_share <= 0.20
##    - else Shared
## ---------------------------
weak_thr <- 2.0
rna_thr  <- 0.80
meth_thr <- 0.20

ve_df2 <- ve_df %>%
  mutate(
    RNA_share = VE_RNA / pmax(VE_TOTAL, 1e-9),
    Class = case_when(
      VE_TOTAL < weak_thr   ~ "Weak",
      RNA_share >= rna_thr  ~ "RNA-specific",
      RNA_share <= meth_thr ~ "METH-specific",
      TRUE                  ~ "Shared"
    )
  )

saveRDS(ve_df2, file.path(RDS_DIR, "GFA_factor_classes.rds"))
readr::write_csv(ve_df2, file.path(TABLES_DIR, "GFA_factor_classes.csv"))
message("Saved factor classes.")

## ---------------------------
## 5) Save top-15 Z/W using VE rank order
## ---------------------------
TOPK <- min(15L, K)

top_factors <- ve_df2 %>%
  arrange(Rank) %>%
  slice_head(n = TOPK)

# Which original columns to keep (by factor label)
keep <- top_factors$Factor

Z_top <- Z[, keep, drop = FALSE]
W_top <- list(
  RNA  = W_RNA[,  keep, drop = FALSE],
  METH = W_METH[, keep, drop = FALSE]
)

# Rename columns F1..F_TOPK consistently (like MOFA)
new_names <- paste0("F", seq_len(ncol(Z_top)))
colnames(Z_top) <- new_names
W_top <- lapply(W_top, function(m) { colnames(m) <- new_names; m })

saveRDS(Z_top, file.path(RDS_DIR, "z_top15.rds"))
saveRDS(W_top, file.path(RDS_DIR, "w_top15.rds"))

message("Saved top15: z_top15.rds, w_top15.rds")

## ---------------------------
## 6) Plot: VE stacked per factor (ranked)
## ---------------------------
ve_long <- ve_df2 %>%
  dplyr::arrange(dplyr::desc(VE_TOTAL)) %>%
  dplyr::mutate(Factor = factor(Factor, levels = Factor)) %>%
  dplyr::select(Factor, VE_RNA, VE_METH) %>%
  tidyr::pivot_longer(
    cols = c(VE_RNA, VE_METH),
    names_to = "View",
    values_to = "Percent"
  ) %>%
  dplyr::mutate(
    View = dplyr::recode(View, VE_RNA = "RNA", VE_METH = "METH"),
    View = factor(View, levels = c("RNA", "METH"))
  )

p1 <- ggplot2::ggplot(ve_long, ggplot2::aes(x = Factor, y = Percent, fill = View)) +
  ggplot2::geom_col(position = ggplot2::position_stack(reverse = TRUE), width = 0.75) +
  ggplot2::scale_fill_manual(
    values = c("RNA" = "#0B3C5D", "METH" = "#CFE1F2"),
    name = "View"
  ) +
  ggplot2::labs(
    title = "GFA: % variance explained per factor (by view)",
    y = "% variance explained",
    x = NULL
  ) +
  ggplot2::theme_minimal(base_size = 13) +
  ggplot2::theme(
    plot.title   = ggplot2::element_text(face = "bold", hjust = 0.5),
    axis.text.x  = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )

print(p1)


ggsave(file.path(PLOTS_DIR, "GFA_VE_stacked_per_factor.png"),
       p1, width = 8, height = 4, dpi = 300, bg = "white")

message("Saved plot: GFA_VE_stacked_per_factor.png")

## ---------------------------
## 7) Plot: VE by factor class (stacked by view)
## ---------------------------
weak_thr <- 2.0    # <2% total VE
rna_thr  <- 0.8    # RNA share ≥ 80%
meth_thr <- 0.2    # RNA share ≤ 20%


ve_df2 <- ve_df2 %>%
  mutate(
    RNA_share = VE_RNA / pmax(VE_TOTAL, 1e-9),
    Class = case_when(
      VE_TOTAL < weak_thr        ~ "weak",
      RNA_share >= rna_thr       ~ "RNA-specific",
      RNA_share <= meth_thr      ~ "METH-specific",
      TRUE                       ~ "shared"
    )
  )

counts_txt <- ve_df2 %>%
  count(Class) %>%
  mutate(Class = factor(Class, levels = c("METH-specific","RNA-specific","Shared","Weak"))) %>%
  arrange(Class)

subtitle_txt <- paste0(
  paste0(as.character(counts_txt$Class), "=", counts_txt$n),
  collapse = " | "
)

p2 <- ggplot(ve_class_sum, aes(x = Class, y = Percent, fill = View)) +
  geom_col(position = position_stack(reverse = TRUE), width = 0.75) +
  scale_fill_manual(
    values = c(
      "RNA"  = "#0B3C5D",
      "METH" = "#CFE1F2"
    ),
    name = "View"
  ) +
  labs(
    title    = "GFA: Explained variance by factor class (by view)",
    subtitle = subtitle_txt,
    y = "Total % variance explained",
    x = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "right"
  )

print(p2)
ggsave(
  file.path(PLOTS_DIR, "GFA_VE_by_factor_class.png"),
  p2,
  width = 7,
  height = 4,
  dpi = 300,
  bg = "white"
)

message("Saved plot: GFA_VE_by_factor_class.png")
message("STEP 04 DONE")

## ============================================================
## 01_gfa_rna_gsea.R — GFA RNA GSEA (GO:BP)
## Inputs:
##   results/GFA/rds/w_top15.rds   (list: RNA, METH; features x K)
## Outputs:
##   results/GFA/GSEA_RNA/tables/gfa_rna_gsea_f*.csv
##   results/GFA/GSEA_RNA/gfa_rna_gsea_sig_counts.csv
## ============================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(dplyr)
  library(readr)
})

## ---------------------------
## paths
## ---------------------------
root_dir <- file.path("results", "GFA")
rds_dir  <- file.path(root_dir, "rds")

in_w <- file.path(rds_dir, "w_top15.rds")
stopifnot(file.exists(in_w))

out_dir   <- file.path(root_dir, "GSEA")
tables_dir <- file.path(out_dir, "tables")

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

## ---------------------------
## load W (top15) and take RNA
## ---------------------------
W_top <- readRDS(in_w)
stopifnot(is.list(W_top), all(c("RNA","METH") %in% names(W_top)))

W_RNA <- as.matrix(W_top[["RNA"]])
storage.mode(W_RNA) <- "double"
stopifnot(!is.null(rownames(W_RNA)), !is.null(colnames(W_RNA)))

message("W_RNA dim: ", paste(dim(W_RNA), collapse = " x "))

## ---------------------------
## helpers: SYMBOL -> ENTREZ and rank prep
## ---------------------------
sym2ent <- function(sym) {
  sym <- unique(na.omit(as.character(sym)))
  if (length(sym) == 0L) return(character(0))
  
  sym_up   <- toupper(sym)
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
  ent  <- sym2ent(names(w_named))
  keep <- !is.na(ent) & !duplicated(ent)
  if (!any(keep)) return(numeric(0))
  r <- setNames(as.numeric(w_named[keep]), ent[keep])
  r[order(r, decreasing = TRUE)]
}

## ---------------------------
## factor filter: sd > 0
## ---------------------------
factor_sd <- apply(W_RNA, 2, sd, na.rm = TRUE)
good_factors <- which(factor_sd > 0)

message("RNA factors with sd>0: ", paste(colnames(W_RNA)[good_factors], collapse = ", "))

## ---------------------------
## GSEA loop
## ---------------------------
gsea_list <- list()

for (k in good_factors) {
  
  f <- colnames(W_RNA)[k]
  
  w <- W_RNA[, k]
  names(w) <- rownames(W_RNA)
  
  if (max(abs(w), na.rm = TRUE) < 1e-6) {
    message("Skip ", f, ": near-zero weights")
    next
  }
  
  ranks <- prep_gsea_rank(w)
  message("=== GFA RNA ", f, " | ranked genes: ", length(ranks), " ===")
  
  if (length(ranks) < 50) {
    message("Skip ", f, ": <50 genes after mapping")
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
      message("ERROR ", f, ": ", e$message)
      NULL
    }
  )
  
  if (!is.null(gse) && nrow(as.data.frame(gse)) > 0) {
    df <- as.data.frame(gse)
    df$factor <- f
    
    out_csv <- file.path(tables_dir, paste0("gfa_rna_gsea_", tolower(f), ".csv"))
    readr::write_csv(df, out_csv)
    
    gsea_list[[f]] <- df
    message("Saved: ", out_csv)
  } else {
    message("No significant GO:BP terms for ", f)
  }
}

## ---------------------------
## summary counts (FDR<=0.05)
## ---------------------------
if (length(gsea_list) > 0) {
  gsea_all <- dplyr::bind_rows(gsea_list, .id = "FactorID") %>%
    mutate(p.adjust = suppressWarnings(as.numeric(p.adjust))) %>%
    filter(!is.na(p.adjust), p.adjust <= 0.05)
  
  gsea_counts <- gsea_all %>%
    group_by(Factor = FactorID) %>%
    summarise(sig_terms = n(), .groups = "drop") %>%
    arrange(desc(sig_terms))
  
  out_counts <- file.path(out_dir, "gfa_rna_gsea_sig_counts.csv")
  readr::write_csv(gsea_counts, out_counts)
  message("Saved summary: ", out_counts)
} else {
  message("No RNA GSEA outputs produced.")
}

message("01_gfa_rna_gsea DONE")

## ============================================================
## 01_gfa_rna_gsea_sig_counts_plot.R
## Purpose:
##   - Plot: number of significant GO:BP terms per GFA RNA factor
## Input :
##   - results/GFA/GSEA/gfa_rna_gsea_sig_counts.csv
## Output:
##   - results/GFA/GSEA/plots/GFA_RNA_GSEA_sig_counts.png
## ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(forcats)
})

IN_FILE  <- file.path("results", "GFA", "GSEA", "gfa_rna_gsea_sig_counts.csv")
OUT_DIR  <- file.path("results", "GFA", "GSEA", "plots")
OUT_PNG  <- file.path(OUT_DIR, "GFA_RNA_GSEA_sig_counts.png")

stopifnot(file.exists(IN_FILE))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

df <- readr::read_csv(IN_FILE, show_col_types = FALSE)

## Beklenen kolon adları: Factor, sig_terms
## (Eğer sende farklıysa: names(df) ile bak, aşağıyı ona göre map’leriz.)
stopifnot(all(c("Factor", "sig_terms") %in% names(df)))

df_plot <- df %>%
  mutate(
    sig_terms = as.numeric(sig_terms),
    Factor = as.character(Factor)
  ) %>%
  filter(!is.na(sig_terms), sig_terms > 0) %>%
  arrange(desc(sig_terms)) %>%
  mutate(Factor = forcats::fct_inorder(Factor))

p <- ggplot(df_plot, aes(x = sig_terms, y = forcats::fct_rev(Factor), fill = sig_terms)) +
  geom_col(width = 0.75) +
  scale_fill_gradient(
    low  = "#CFE1F2",
    high = "#0B3C5D",
    name = "GO:BP term count"
  ) +
  labs(
    title = "GFA RNA factors (GSEA, FDR \u2264 0.05): significant GO:BP terms",
    x = "Number of enriched GO:BP terms",
    y = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  )

print(p)

ggsave(OUT_PNG, p, width = 9, height = 5, dpi = 300, bg = "white")
message("Saved: ", OUT_PNG)


## =========================================================
## 02_gfa_ora_meth_go_bp.R — GFA METH ORA (GO:BP) across active factors
## Purpose:
##   - Run ORA on GFA METH factor weights using top-|weights| genes
## Input:
##   - results/GFA/rds/W_all.rds   (expects list: $RNA, $METH; each features x K)
## Output:
##   - results/GFA/ORA/tables/gfa_meth_ora_f*.csv
##   - results/GFA/ORA/tables/gfa_meth_ora_sig_counts.csv
##   - results/GFA/ORA/plots/gfa_meth_ora_sig_terms_horizontal.png
## =========================================================

source("scripts/00_setup_paths.R")
stopifnot(exists("gfa_dir"), dir.exists(gfa_dir))

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
## Output dirs (GFA ORA) — keep lowercase filenames
## ---------------------------
ora_meth_dir  <- file.path(gfa_dir, "ORA", "tables")
ora_plots_dir <- file.path(gfa_dir, "ORA", "plots")

dir.create(ora_meth_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(ora_plots_dir, recursive = TRUE, showWarnings = FALSE)

## ---------------------------
## Load METH weights (features x K)
## ---------------------------
w_all_file <- file.path(gfa_dir, "rds", "W_all.rds")
stopifnot(file.exists(w_all_file))

W_all  <- readRDS(w_all_file)
stopifnot(is.list(W_all), all(c("RNA","METH") %in% names(W_all)))

W_METH <- as.matrix(W_all[["METH"]])
storage.mode(W_METH) <- "double"

stopifnot(!is.null(rownames(W_METH)), !is.null(colnames(W_METH)))

# Clean gene symbols (safe)
rownames(W_METH) <- toupper(rownames(W_METH))
rownames(W_METH) <- gsub("_METH$", "", rownames(W_METH))
rownames(W_METH) <- gsub("_RNA$",  "", rownames(W_METH))
rownames(W_METH) <- gsub("\\|.*$", "", rownames(W_METH))

K <- ncol(W_METH)
message("Loaded GFA W_METH: ", paste(dim(W_METH), collapse = " x "), " (K=", K, ")")

## ---------------------------
## Active factor detection
## ---------------------------
eps <- 1e-8
meth_sd  <- apply(W_METH, 2, sd, na.rm = TRUE)
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
  
  sym_up <- toupper(sym)
  
  # mapIds is safer than select() for 1:many situations
  res <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys      = sym_up,
    keytype   = "SYMBOL",
    column    = "ENTREZID",
    multiVals = "first"
  )
  
  unique(na.omit(as.character(res)))
}

prep_ora_genes <- function(w_named, topN = 300) {
  stopifnot(!is.null(names(w_named)))
  ord <- order(abs(w_named), decreasing = TRUE)
  top_sym <- names(w_named)[ord][seq_len(min(topN, length(ord)))]
  sym2ent(top_sym)
}

## Universe (recommended for ORA): all genes present in this view
universe_entrez <- sym2ent(rownames(W_METH))
message("Universe ENTREZ (METH): ", length(universe_entrez))

## ---------------------------
## Run ORA per active factor
## ---------------------------
TOPN_ORA <- 300
ora_results_list <- list()

for (f in active_factors) {
  
  w <- W_METH[, f]
  names(w) <- rownames(W_METH)
  
  genes_entrez <- prep_ora_genes(w, topN = TOPN_ORA)
  
  message("=== ORA GFA METH ", f, " | topN=", TOPN_ORA, " | ENTREZ genes: ", length(genes_entrez), " ===")
  
  if (length(genes_entrez) < 20) {
    message("Skip ", f, ": <20 ENTREZ genes after mapping.")
    next
  }
  
  enr <- tryCatch(
    enrichGO(
      gene          = genes_entrez,
      universe      = universe_entrez,
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
    
    out_csv <- file.path(ora_meth_dir, paste0("gfa_meth_ora_", tolower(f), ".csv"))
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
  
  out_counts <- file.path(ora_meth_dir, "gfa_meth_ora_sig_counts.csv")
  readr::write_csv(ora_counts, out_counts)
  message("Saved summary counts: ", out_counts)
  
  ## Plot (same palette as MOFA)
  ora_counts <- ora_counts %>%
    mutate(Factor = factor(Factor, levels = rev(Factor)))
  
  p <- ggplot(ora_counts, aes(x = sig_terms, y = Factor, fill = sig_terms)) +
    geom_col(width = 0.8) +
    scale_fill_gradient(low = "#cfe1f2", high = "#0b3c5d", name = "GO:BP term count") +
    labs(
      title = "GFA methylation factors (ORA, FDR \u2264 0.05): significant GO:BP terms",
      x = "Number of enriched GO:BP terms",
      y = NULL
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )
  
  print(p)
  
  out_png <- file.path(ora_plots_dir, "gfa_meth_ora_sig_terms_horizontal.png")
  ggsave(out_png, p, width = 10, height = 5, dpi = 300, bg = "white")
  message("Saved ORA plot: ", out_png)
  
} else {
  message("No METH ORA results were produced (no active factors or no significant terms).")
}

message("GFA METH ORA (GO:BP) DONE.")
