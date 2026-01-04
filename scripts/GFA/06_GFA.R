## =========================================================
## GFA.R — STEP 0: Paths + libraries (GitHub/repo-safe)
## Purpose:
##   - Define repo-root-safe paths (no setwd to Desktop)
##   - Create fixed output folders for GFA run
##   - Load required libraries
##   - Sanity-check that MOFA inputs exist
## =========================================================

## ---------------------------
## 0A) PATHS (project-safe; ONLY from 00_setup_paths.R)
## ---------------------------
source("scripts/00_setup_paths.R")

# Canonicalize (use setup aliases if needed)
if (!exists("ROOT") && exists("project_dir")) ROOT <- project_dir
if (!exists("RESULTS_DIR") && exists("results_dir")) RESULTS_DIR <- results_dir
if (!exists("RESULTS_DIR")) stop("RESULTS_DIR/results_dir not found. Check scripts/00_setup_paths.R")

RUN_ID <- "run_paired01"

## ---- Inputs: MOFA input objects must live in setup-defined paired_inputs_dir ----
if (!exists("paired_inputs_dir")) stop("paired_inputs_dir not found. Check scripts/00_setup_paths.R")
IN_DIR <- paired_inputs_dir
stopifnot(dir.exists(IN_DIR))

## ---- Outputs: use setup-defined gfa_dir (results/GFA) ----
if (!exists("gfa_dir")) stop("gfa_dir not found. Check scripts/00_setup_paths.R")

GFA_RUN_DIR <- file.path(gfa_dir, RUN_ID)
PLOTS_DIR   <- file.path(GFA_RUN_DIR, "plots")
TABLES_DIR  <- file.path(GFA_RUN_DIR, "tables")
RDS_DIR     <- file.path(GFA_RUN_DIR, "rds")

dir.create(GFA_RUN_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOTS_DIR,   recursive = TRUE, showWarnings = FALSE)
dir.create(TABLES_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(RDS_DIR,     recursive = TRUE, showWarnings = FALSE)

## ---- MOFA fixed input file paths (stored under paired_inputs_dir) ----
RNA_IN  <- file.path(IN_DIR, "rna_mofa_input.rds")
METH_IN <- file.path(IN_DIR, "meth_mofa_input.rds")
META_IN <- file.path(IN_DIR, "sample_meta_mofa.rds")  # optional

message("ROOT:         ", ROOT)
message("RESULTS_DIR:  ", RESULTS_DIR)
message("IN_DIR:       ", IN_DIR)
message("GFA_RUN_DIR:  ", GFA_RUN_DIR)
message("PLOTS_DIR:    ", PLOTS_DIR)
message("TABLES_DIR:   ", TABLES_DIR)
message("RDS_DIR:      ", RDS_DIR)
message("RNA_IN:       ", RNA_IN)
message("METH_IN:      ", METH_IN)
message("META_IN:      ", META_IN)

stopifnot(file.exists(RNA_IN))
stopifnot(file.exists(METH_IN))
# META_IN is optional (checked later)

## ---------------------------
## 0C) LIBRARIES
## ---------------------------

suppressPackageStartupMessages({
  library(GFA)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(forcats)
  library(pheatmap)
  library(grid)
})

message("STEP 0 DONE: paths + libraries loaded + inputs found.")

## =========================================================
## GFA.R — STEP 1: Load MOFA inputs -> build GFA X_list -> scale -> save
## Purpose:
##   - Read paired MOFA inputs (features x samples)
##   - Enforce identical sample order across views
##   - Transpose to GFA format (samples x features)
##   - Scale each feature (mean=0, sd=1), replace NA (zero-variance) with 0
##   - Save prepared objects for reproducibility
## Outputs:
##   results/GFA/run_paired01/rds/
##     - X_list_scaled.rds
##     - rna_X_scaled.rds
##     - meth_X_scaled.rds
##     - sample_meta_gfa.rds
## =========================================================

## ---------------------------
## 1A) Helpers
## ---------------------------

to_num_matrix <- function(x) {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x)) stop("Input is not a matrix/data.frame.")
  storage.mode(x) <- "double"
  x
}

scale_and_clean <- function(m) {
  m <- scale(m, center = TRUE, scale = TRUE)
  m[is.na(m)] <- 0
  m
}

## ---------------------------
## 1B) Load MOFA inputs (paired)
## ---------------------------

rna_mat  <- to_num_matrix(readRDS(RNA_IN))   # features x samples
meth_mat <- to_num_matrix(readRDS(METH_IN))  # features x samples

message("Loaded RNA  (features x samples):  ", paste(dim(rna_mat),  collapse = " x "))
message("Loaded METH (features x samples):  ", paste(dim(meth_mat), collapse = " x "))

stopifnot(!is.null(colnames(rna_mat)), !is.null(colnames(meth_mat)))

## ---------------------------
## 1C) Enforce identical sample set + order
## ---------------------------

common_samples <- intersect(colnames(rna_mat), colnames(meth_mat))
if (length(common_samples) == 0L) {
  stop("No overlapping samples between RNA and METH (check colnames in MOFA inputs).")
}

# Use RNA order as the reference (stable + reproducible)
common_samples <- colnames(rna_mat)[colnames(rna_mat) %in% common_samples]

rna_mat  <- rna_mat[,  common_samples, drop = FALSE]
meth_mat <- meth_mat[, common_samples, drop = FALSE]
meth_mat <- meth_mat[, colnames(rna_mat), drop = FALSE]

stopifnot(identical(colnames(rna_mat), colnames(meth_mat)))
message("Paired samples aligned (n): ", ncol(rna_mat))

## ---------------------------
## 1D) Build GFA inputs (samples x features) + scale
## ---------------------------

X_rna  <- t(rna_mat)   # samples x RNA features
X_meth <- t(meth_mat)  # samples x METH features

X_rna_scaled  <- scale_and_clean(X_rna)
X_meth_scaled <- scale_and_clean(X_meth)

stopifnot(nrow(X_rna_scaled) == nrow(X_meth_scaled))
message("GFA-ready matrices (samples x features):")
message("  RNA :  ", paste(dim(X_rna_scaled),  collapse = " x "))
message("  METH:  ", paste(dim(X_meth_scaled), collapse = " x "))

X_list <- list(
  RNA  = X_rna_scaled,
  METH = X_meth_scaled
)

## ---------------------------
## 1E) Sample metadata (optional but saved)
## ---------------------------

sample_meta <- if (file.exists(META_IN)) {
  readRDS(META_IN)
} else {
  data.frame(sample_id = rownames(X_rna_scaled), stringsAsFactors = FALSE)
}

# Ensure rownames exist for samples (nice for downstream)
if (is.null(rownames(X_rna_scaled))) rownames(X_rna_scaled) <- common_samples
if (is.null(rownames(X_meth_scaled))) rownames(X_meth_scaled) <- common_samples
rownames(X_rna_scaled)  <- common_samples
rownames(X_meth_scaled) <- common_samples

## ---------------------------
## 1F) Save prepared inputs (fixed names)
## ---------------------------

saveRDS(X_list,         file.path(RDS_DIR, "X_list_scaled.rds"))
saveRDS(X_rna_scaled,   file.path(RDS_DIR, "rna_X_scaled.rds"))
saveRDS(X_meth_scaled,  file.path(RDS_DIR, "meth_X_scaled.rds"))
saveRDS(sample_meta,    file.path(RDS_DIR, "sample_meta_gfa.rds"))

message("Saved prepared inputs to: ", RDS_DIR)
message("STEP 1 DONE")

## =========================================================
## GFA.R — STEP 2: Fit GFA model (from scratch) + save gfa_fit
## Purpose:
##   - Load prepared GFA inputs (X_list_scaled.rds)
##   - Fit GFA with a fixed K (start with 15 to match MOFA)
##   - Save fitted model object for reproducibility
## Outputs:
##   results/GFA/run_paired01/rds/
##     - gfa_fit.rds
##     - gfa_runtime_seconds.txt
## =========================================================

## ---------------------------
## 2A) Load prepared inputs
## ---------------------------

X_list <- readRDS(file.path(RDS_DIR, "X_list_scaled.rds"))
stopifnot(is.list(X_list), all(c("RNA","METH") %in% names(X_list)))

# Safety: dimensions must match in N (samples)
stopifnot(nrow(X_list$RNA) == nrow(X_list$METH))

N  <- nrow(X_list$RNA)
Dr <- ncol(X_list$RNA)
Dm <- ncol(X_list$METH)

message("Loaded X_list (samples x features):")
message("  RNA :  ", N, " x ", Dr)
message("  METH:  ", N, " x ", Dm)

## ---------------------------
## 2B) GFA options (reproducible)
## ---------------------------

set.seed(42)

opts <- GFA::getDefaultOpts()

# Iterations (moderate; adjust later if needed)
opts$iter.max    <- 2000
opts$iter.burnin <- 800
opts$verbose     <- 2

# Number of factors (match MOFA first)
K <- 15
opts$k <- K
opts$K <- K

# Noise prior (recommended for real-valued matrices)
opts <- GFA::informativeNoisePrior(X_list, opts)

message("Fitting GFA with K = ", K, " factors ...")

## ---------------------------
## 2C) Fit model + runtime
## ---------------------------

t0 <- Sys.time()

gfa_fit <- GFA::gfa(X_list, opts, K = K)

t1 <- Sys.time()
runtime_sec <- as.numeric(difftime(t1, t0, units = "secs"))

message("GFA fit done. Runtime (sec): ", round(runtime_sec, 1))

## ---------------------------
## 2D) Save fitted model
## ---------------------------

saveRDS(gfa_fit, file.path(RDS_DIR, "gfa_fit.rds"))
writeLines(as.character(runtime_sec), file.path(RDS_DIR, "gfa_runtime_seconds.txt"))

message("Saved: ", file.path(RDS_DIR, "gfa_fit.rds"))
message("STEP 2 DONE")


# how many active factors  defined: 
gfa_fit <- readRDS(file.path(RDS_DIR, "gfa_fit.rds"))

## --- helper: extract W (global loadings) safely ---
get_W_matrix <- function(fit) {
  if (!is.null(fit$W))        return(as.matrix(fit$W))
  if (!is.null(fit$w))        return(as.matrix(fit$w))
  if (!is.null(fit$loadings)) return(as.matrix(fit$loadings))
  stop("Could not find W in gfa_fit (tried $W, $w, $loadings).")
}

W_all <- get_W_matrix(gfa_fit)  # (D_total x K)
K <- ncol(W_all)
cat("W_all dim (D_total x K):", nrow(W_all), "x", K, "\n")

## --- get group indices to split W into RNA / METH blocks ---
if (is.null(gfa_fit$groups)) stop("gfa_fit$groups not found; cannot split W into views.")
grp <- gfa_fit$groups

## Ensure names match our views
if (is.null(names(grp))) names(grp) <- names(X_list)
stopifnot(all(c("RNA","METH") %in% names(grp)))

W_RNA  <- W_all[ grp[["RNA"]],  , drop = FALSE ]
W_METH <- W_all[ grp[["METH"]], , drop = FALSE ]

## --- define activity by threshold on loadings ---
eps <- 1e-6  # activity threshold
rna_active  <- apply(W_RNA,  2, function(w) any(abs(w) > eps, na.rm = TRUE))
meth_active <- apply(W_METH, 2, function(w) any(abs(w) > eps, na.rm = TRUE))

class <- ifelse(rna_active & meth_active, "Shared",
                ifelse(rna_active & !meth_active, "RNA-specific",
                       ifelse(!rna_active & meth_active, "METH-specific", "Inactive")))

## --- summary ---
tab <- table(class)
print(tab)

cat("\nActive factors (non-Inactive):", sum(class != "Inactive"), "out of", K, "\n")
cat("Inactive factors:", sum(class == "Inactive"), "out of", K, "\n")

## Optional: list which factors are inactive
inactive_idx <- which(class == "Inactive")
if (length(inactive_idx) > 0) {
  cat("Inactive factor indices:", paste(inactive_idx, collapse = ", "), "\n")
}

## =========================================================
## GFA.R — STEP 3A: Extract W (loadings) from fitted GFA model
## Purpose:
##   - Extract global W from gfa_fit (D_total x K)
##   - Split W into W_RNA and W_METH using gfa_fit$groups
##   - Save W objects with fixed filenames
## Outputs (RDS_DIR):
##   - W_all.rds
##   - W_RNA_all.rds
##   - W_METH_all.rds
## =========================================================

## ---------------------------
## 3A.1 Load fitted model + inputs (for safety)
## ---------------------------

gfa_fit <- readRDS(file.path(RDS_DIR, "gfa_fit.rds"))
X_list  <- readRDS(file.path(RDS_DIR, "X_list_scaled.rds"))

stopifnot(is.list(X_list), all(c("RNA","METH") %in% names(X_list)))

## ---------------------------
## 3A.2 Extract global W
## ---------------------------

get_W_from_gfa <- function(fit) {
  if (!is.null(fit$W))        return(as.matrix(fit$W))
  if (!is.null(fit$w))        return(as.matrix(fit$w))
  if (!is.null(fit$loadings)) return(as.matrix(fit$loadings))
  stop("Could not find W in gfa_fit (tried $W, $w, $loadings).")
}

W_all <- get_W_from_gfa(gfa_fit)

## Safety: numeric matrix
W_all <- as.matrix(W_all)
storage.mode(W_all) <- "double"

K <- ncol(W_all)
D <- nrow(W_all)

message("Extracted W_all (D_total x K): ", D, " x ", K)

## ---------------------------
## 3A.3 Split W into views using group indices
## ---------------------------

if (is.null(gfa_fit$groups)) stop("gfa_fit$groups not found. Cannot split W by view.")
grp <- gfa_fit$groups

# Ensure group names match the X_list view names
if (is.null(names(grp))) names(grp) <- names(X_list)

stopifnot(all(c("RNA","METH") %in% names(grp)))

W_RNA  <- W_all[ grp[["RNA"]],  , drop = FALSE ]
W_METH <- W_all[ grp[["METH"]], , drop = FALSE ]

message("W_RNA  (D_RNA x K):  ", nrow(W_RNA),  " x ", ncol(W_RNA))
message("W_METH (D_METH x K): ", nrow(W_METH), " x ", ncol(W_METH))

## Safety: feature counts should match X_list columns
stopifnot(nrow(W_RNA)  == ncol(X_list$RNA))
stopifnot(nrow(W_METH) == ncol(X_list$METH))

## ---------------------------
## 3A.4 Save with fixed filenames
## ---------------------------

saveRDS(W_all,  file.path(RDS_DIR, "W_all.rds"))
saveRDS(W_RNA,  file.path(RDS_DIR, "W_RNA_all.rds"))
saveRDS(W_METH, file.path(RDS_DIR, "W_METH_all.rds"))

message("Saved:")
message("  ", file.path(RDS_DIR, "W_all.rds"))
message("  ", file.path(RDS_DIR, "W_RNA_all.rds"))
message("  ", file.path(RDS_DIR, "W_METH_all.rds"))

message("STEP 3A DONE")

## =========================================================
## GFA.R — STEP 3B: Check and extract Z (factor scores)
## Purpose:
##   - Check whether the fitted GFA model provides Z directly
##   - Extract Z ONLY if it exists in the model object
##   - DO NOT compute or project Z post-hoc
## Outputs (RDS_DIR, if available):
##   - Z_all.rds   (samples x K)
## =========================================================

## ---------------------------
## 3B.1 Load fitted model
## ---------------------------

gfa_fit <- readRDS(file.path(RDS_DIR, "gfa_fit.rds"))
X_list  <- readRDS(file.path(RDS_DIR, "X_list_scaled.rds"))

N_expected <- nrow(X_list$RNA)
K_expected <- ncol(readRDS(file.path(RDS_DIR, "W_all.rds")))

message("Expected N (samples): ", N_expected)
message("Expected K (factors): ", K_expected)

## --- helper: safely get dims for candidate objects ---
dim_safe <- function(x) {
  if (is.null(x)) return(c(NA, NA))
  if (is.vector(x)) return(c(length(x), NA))
  d <- dim(x)
  if (is.null(d)) c(NA, NA) else d
}

## --- list candidate slots commonly used across versions ---
candidates <- list(
  Z = gfa_fit$Z,
  z = gfa_fit$z,
  X = gfa_fit$X,
  x = gfa_fit$x,
  scores = gfa_fit$scores,
  latent = gfa_fit$latent,
  E = gfa_fit$E
)

cand_dims <- data.frame(
  slot = names(candidates),
  nrow = sapply(candidates, function(a) dim_safe(a)[1]),
  ncol = sapply(candidates, function(a) dim_safe(a)[2])
)

print(cand_dims)

## --- pick candidate that matches expected sample x factor shape ---
ok <- which(cand_dims$nrow == N_expected & cand_dims$ncol == K_expected)

if (length(ok) == 0) {
  message("No Z found with dim N x K = ", N_expected, " x ", K_expected)
  message("We will proceed WITHOUT model-provided Z.")
} else {
  slot_name <- cand_dims$slot[ok[1]]
  Z_all <- as.matrix(candidates[[slot_name]])
  storage.mode(Z_all) <- "double"
  
  message("Using slot: ", slot_name, "  dim = ", nrow(Z_all), " x ", ncol(Z_all))
  
  saveRDS(Z_all, file.path(RDS_DIR, "Z_all.rds"))
  message("Saved corrected Z_all: ", file.path(RDS_DIR, "Z_all.rds"))
}

message("STEP 3B FIX DONE")

## =========================================================
## GFA.R — STEP 4A: Variance explained (MOFA-like) + rank factors
## GFA.R — STEP 4B: Factor class assignment (sparsity-based)
## Purpose:
##   - Compute VE per factor per view using Z and W
##   - Rank factors by total VE and label as F1..FK
##   - Create factor classes (RNA-specific / METH-specific / Shared / Inactive)
## Outputs (RDS_DIR):
##   - VE_per_factor_ranked.rds
##   - VE_per_factor_ranked.csv
##   - Z_ranked_F1toFK.rds
##   - W_RNA_ranked_F1toFK.rds
##   - W_METH_ranked_F1toFK.rds
##   - GFA_factor_classes_sparsity.csv
## =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

## ---------------------------
## STEP 4A.1 Load inputs
## ---------------------------

X_list <- readRDS(file.path(RDS_DIR, "X_list_scaled.rds"))
Z_all  <- readRDS(file.path(RDS_DIR, "Z_all.rds"))
W_RNA  <- readRDS(file.path(RDS_DIR, "W_RNA_all.rds"))
W_METH <- readRDS(file.path(RDS_DIR, "W_METH_all.rds"))

to_num <- function(m) { m <- as.matrix(m); storage.mode(m) <- "double"; m }

X_rna  <- to_num(X_list$RNA)   # N x D_r
X_meth <- to_num(X_list$METH)  # N x D_m
Z_all  <- to_num(Z_all)        # N x K
W_RNA  <- to_num(W_RNA)        # D_r x K
W_METH <- to_num(W_METH)       # D_m x K

## ---------------------------
## STEP 4A.2 Safety checks
## ---------------------------

stopifnot(
  nrow(X_rna) == nrow(X_meth),
  nrow(Z_all) == nrow(X_rna),
  ncol(Z_all) == ncol(W_RNA),
  ncol(Z_all) == ncol(W_METH),
  nrow(W_RNA) == ncol(X_rna),
  nrow(W_METH) == ncol(X_meth)
)

N <- nrow(Z_all)
K <- ncol(Z_all)

message("N (samples): ", N)
message("K (factors): ", K)
message("RNA  X dim (N x D_r):  ", paste(dim(X_rna),  collapse = " x "))
message("METH X dim (N x D_m):  ", paste(dim(X_meth), collapse = " x "))

## ---------------------------
## STEP 4A.3 Center data (for VE denominator)
## ---------------------------

Xr_c <- scale(X_rna,  center = TRUE, scale = FALSE)
Xm_c <- scale(X_meth, center = TRUE, scale = FALSE)

tot_var_rna  <- sum(Xr_c^2)
tot_var_meth <- sum(Xm_c^2)

## ---------------------------
## STEP 4A.4 Compute VE per factor per view
## ---------------------------

ve_rna  <- numeric(K)
ve_meth <- numeric(K)

for (k in seq_len(K)) {
  zk <- Z_all[, k, drop = FALSE]  # N x 1
  
  Xhat_r <- zk %*% t(W_RNA[,  k, drop = FALSE])   # N x D_r
  Xhat_m <- zk %*% t(W_METH[, k, drop = FALSE])   # N x D_m
  
  ve_rna[k]  <- sum(Xhat_r^2) / tot_var_rna
  ve_meth[k] <- sum(Xhat_m^2) / tot_var_meth
}

ve_df <- data.frame(
  k_index  = seq_len(K),
  VE_RNA   = ve_rna  * 100,
  VE_METH  = ve_meth * 100,
  stringsAsFactors = FALSE
)
ve_df$VE_TOTAL <- ve_df$VE_RNA + ve_df$VE_METH

## ---------------------------
## STEP 4A.5 Rank by VE_TOTAL and relabel F1..FK
## ---------------------------

ve_df <- ve_df[order(ve_df$VE_TOTAL, decreasing = TRUE), , drop = FALSE]
ve_df$Rank   <- seq_len(nrow(ve_df))
ve_df$Factor <- paste0("F", ve_df$Rank)

rank_idx <- ve_df$k_index  # original factor indices in ranked order

Z_ranked      <- Z_all[,  rank_idx, drop = FALSE]
W_RNA_ranked  <- W_RNA[,  rank_idx, drop = FALSE]
W_METH_ranked <- W_METH[, rank_idx, drop = FALSE]

colnames(Z_ranked)      <- ve_df$Factor
colnames(W_RNA_ranked)  <- ve_df$Factor
colnames(W_METH_ranked) <- ve_df$Factor

## ---------------------------
## STEP 4A.6 Save VE + ranked Z/W
## ---------------------------

saveRDS(ve_df, file.path(RDS_DIR, "VE_per_factor_ranked.rds"))
write.csv(ve_df, file.path(RDS_DIR, "VE_per_factor_ranked.csv"), row.names = FALSE)

saveRDS(Z_ranked,      file.path(RDS_DIR, "Z_ranked_F1toFK.rds"))
saveRDS(W_RNA_ranked,  file.path(RDS_DIR, "W_RNA_ranked_F1toFK.rds"))
saveRDS(W_METH_ranked, file.path(RDS_DIR, "W_METH_ranked_F1toFK.rds"))

message("STEP 4A DONE: saved VE + ranked Z/W into: ", RDS_DIR)

## =========================================================
## STEP 4B: Factor class assignment (sparsity-based)
## IMPORTANT: This produces the CSV that STEP 5 needs.
## =========================================================

eps <- 1e-6

rna_active  <- apply(W_RNA,  2, function(w) any(abs(w) > eps, na.rm = TRUE))
meth_active <- apply(W_METH, 2, function(w) any(abs(w) > eps, na.rm = TRUE))

factor_class <- ifelse(rna_active & meth_active, "Shared",
                       ifelse(rna_active,              "RNA-specific",
                              ifelse(meth_active,             "METH-specific",
                                     "Inactive")))

class_df <- data.frame(
  Factor = paste0("F", seq_len(K)),  # original factor order
  Class  = factor_class,
  stringsAsFactors = FALSE
)

print(table(class_df$Class))

write.csv(
  class_df,
  file.path(RDS_DIR, "GFA_factor_classes_sparsity.csv"),
  row.names = FALSE
)

message("STEP 4B DONE: saved factor classes to: ",
        file.path(RDS_DIR, "GFA_factor_classes_sparsity.csv"))

## =========================================================
## GFA.R — STEP 5: PLOTS (MOFA-like)
## Purpose:
##   1) VE stacked barplot (Top-15, by view)
##   2) VE by factor class (by-view)
## Outputs (PLOTS_DIR):
##   - GFA_VE_stacked_top15.png
##   - GFA_VE_by_class_top15.png
## =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(forcats)
  library(readr)
})

## ---------------------------
## 5.1 Load ranked VE + factor classes
## ---------------------------

ve_df <- readRDS(file.path(RDS_DIR, "VE_per_factor_ranked.rds"))

class_df <- read.csv(file.path(RDS_DIR, "GFA_factor_classes_sparsity.csv"),
                     stringsAsFactors = FALSE)

# Keep Top-15 (you have K=15 anyway, but keep robust)
TOPN <- min(15L, nrow(ve_df))
ve_top <- ve_df %>% slice_head(n = TOPN)

# Join class info (Factor names F1..FK are in ve_top)
ve_top <- ve_top %>%
  left_join(class_df, by = "Factor")

stopifnot(all(!is.na(ve_top$Class)))

## ---------------------------
## 5.2 Plot 1 — VE stacked barplot (Top-15, by view)
## ---------------------------

ve_long <- ve_top %>%
  transmute(
    Factor,
    `RNA expression`  = VE_RNA,
    `DNA methylation` = VE_METH
  ) %>%
  pivot_longer(-Factor, names_to = "View", values_to = "Percent") %>%
  mutate(Factor = fct_relevel(Factor, ve_top$Factor))

p_stacked <- ggplot(ve_long, aes(x = Factor, y = Percent, fill = View)) +
  geom_col() +
  scale_fill_manual(
    values = c(`RNA expression` = "#1B4F72",
               `DNA methylation` = "#cfe1f2"),
    name = "View"
  ) +
  labs(
    title = "GFA: % variance explained per factor (Top-15, by view)",
    y = "% variance explained",
    x = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title  = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )

print(p_stacked)

ggsave(
  filename = file.path(PLOTS_DIR, "GFA_VE_stacked_top15.png"),
  plot     = p_stacked,
  width    = 8,
  height   = 4,
  dpi      = 300,
  bg       = "white"
)

## ---------------------------
## 5.3 Plot 2 — VE by factor class (by-view)
## ---------------------------

# Order classes like MOFA
ve_top$Class <- factor(
  ve_top$Class,
  levels = c("METH-specific", "RNA-specific", "Shared", "Inactive")
)

agg_class <- ve_top %>%
  group_by(Class) %>%
  summarise(
    RNA  = sum(VE_RNA,  na.rm = TRUE),
    METH = sum(VE_METH, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(RNA, METH), names_to = "View", values_to = "Percent") %>%
  mutate(View = factor(View, levels = c("RNA", "METH")))

p_class <- ggplot(agg_class, aes(x = Class, y = Percent, fill = View)) +
  geom_col(position = "stack") +
  scale_fill_manual(
    values = c(RNA = "#0b3c5d", METH = "#cfe1f2"),
    name = "View"
  ) +
  labs(
    title = "GFA: Explained variance by factor class (Top-15, by-view)",
    y = "Total % variance explained",
    x = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title      = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom"
  )

print(p_class)

ggsave(
  filename = file.path(PLOTS_DIR, "GFA_VE_by_class_top15.png"),
  plot     = p_class,
  width    = 7.2,
  height   = 4.2,
  dpi      = 300,
  bg       = "white"
)

message("STEP 5 DONE — plots saved to: ", PLOTS_DIR)


## =========================================================
## STEP 6) GFA — Factor score heatmap + dominance (Top-K)
## Fix: annotation_colors must cover ALL annotation columns
## =========================================================

stopifnot(
  exists("Z_all"), exists("W_RNA"), exists("W_METH"), exists("PLOTS_DIR")
)

## ---- 1) Numeric safety ----
Z_all  <- as.matrix(Z_all);  storage.mode(Z_all)  <- "double"
W_RNA  <- as.matrix(W_RNA);  storage.mode(W_RNA)  <- "double"
W_METH <- as.matrix(W_METH); storage.mode(W_METH) <- "double"

stopifnot(ncol(Z_all) == ncol(W_RNA), ncol(Z_all) == ncol(W_METH))

K <- ncol(Z_all)

## ---- 2) Ensure factor names exist (F1..FK) ----
if (is.null(colnames(Z_all))) {
  colnames(Z_all) <- paste0("F", seq_len(K))
}
if (is.null(colnames(W_RNA)))  colnames(W_RNA)  <- colnames(Z_all)
if (is.null(colnames(W_METH))) colnames(W_METH) <- colnames(Z_all)

## ---- 3) Dominance + contributions (computed from W, aligned to Z) ----
common_factors <- Reduce(intersect, list(colnames(Z_all), colnames(W_RNA), colnames(W_METH)))
stopifnot(length(common_factors) == K)

W_RNA_use  <- W_RNA[,  common_factors, drop = FALSE]
W_METH_use <- W_METH[, common_factors, drop = FALSE]
Z_use      <- Z_all[,  common_factors, drop = FALSE]

rna_act  <- colSums(W_RNA_use^2,  na.rm = TRUE)
meth_act <- colSums(W_METH_use^2, na.rm = TRUE)
tot_act  <- rna_act + meth_act
ratio    <- rna_act / pmax(tot_act, 1e-12)

dom_class <- ifelse(
  ratio >= 0.65, "RNA-dominant",
  ifelse(ratio <= 0.35, "METH-dominant", "Balanced")
)

ann_col <- data.frame(
  Dominance    = factor(dom_class, levels = c("RNA-dominant", "Balanced", "METH-dominant")),
  RNA_contrib  = round(rna_act  / pmax(tot_act, 1e-12), 2),
  METH_contrib = round(meth_act / pmax(tot_act, 1e-12), 2),
  row.names = colnames(Z_use)
)

## ---- 4) Z-score + clip (MOFA style) ----
Zz <- scale(Z_use)

clip <- function(m, lim = 2) {
  m[m >  lim] <-  lim
  m[m < -lim] <- -lim
  m
}
Zz <- clip(Zz, lim = 2)

## ---- 5) Annotation colors (prevents "subscript out of bounds") ----

ann_col <- data.frame(
  Dominance    = dom_class,
  RNA_contrib  = round(rna_act  / pmax(tot_act, 1e-12), 2),
  METH_contrib = round(meth_act / pmax(tot_act, 1e-12), 2)
)
rownames(ann_col) <- colnames(Zz)

## ---- 6) Draw heatmap ----
plots_dir <- PLOTS_DIR
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
hp <- pheatmap(
  Zz,
  color = colorRampPalette(c("#2c7fb8", "#ffffbf", "#d95f0e"))(100),
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  clustering_method = "ward.D2",
  show_rownames = FALSE,
  show_colnames = TRUE,
  angle_col = 90,
  fontsize_col = 10,
  main = paste0("GFA factor scores (Top-15, z-scored)"),
  annotation_col = ann_col,
  annotation_colors = ann_colors,
  border_color = NA,
  silent = TRUE
)

grid::grid.newpage()
grid::grid.draw(hp$gtable)

# save the file:
out_png <- file.path(plots_dir, "GFA_factor_scores_top15_zscored.png")

png(out_png, width = 3200, height = 1800, res = 300, bg = "white")
grid::grid.newpage()
grid::grid.draw(hp$gtable)
dev.off()

stopifnot(file.exists(out_png))
message("Saved heatmap: ", out_png)

## =========================================================
## STEP 7) CROSS-OMICS CORRELATION HEATMAP
## Purpose:
##   - Correlate RNA-specific vs METH-specific factors (Z scores)
##   - Plot clustered correlation matrix (pheatmap)
## Requires (already created earlier):
##   - RDS_DIR, plots_dir
##   - Z_ranked (or Z_all) with colnames = F1..FK
##   - class_df with columns: Factor, Class  (RNA-specific / METH-specific / Shared / Inactive)
## Output:
##   - plots_dir/GFA_crossomics_corr_RNA_vs_METH_specific.png
## =========================================================

suppressPackageStartupMessages({
  library(pheatmap)
  library(grid)
})

## ---- 7.1 Load inputs (use your CURRENT objects if already in memory) ----
if (!exists("Z_ranked")) {
  # If you saved ranked Z from Step 4:
  if (file.exists(file.path(RDS_DIR, "Z_ranked_F1toFK.rds"))) {
    Z_ranked <- readRDS(file.path(RDS_DIR, "Z_ranked_F1toFK.rds"))
  } else {
    # Fallback: use Z_all (then you must ensure factor names match class_df$Factor)
    Z_ranked <- readRDS(file.path(RDS_DIR, "Z_all.rds"))
  }
}
if (!exists("class_df")) {
  # If you saved classes earlier, load it (adjust filename if different in your repo)
  class_path <- file.path(RDS_DIR, "GFA_factor_classes_sparsity.csv")
  class_df <- read.csv(class_path, stringsAsFactors = FALSE)
}

Z_ranked <- as.matrix(Z_ranked); storage.mode(Z_ranked) <- "double"

## ---- 7.2 Ensure factor name alignment (critical) ----
# Expect class_df$Factor like "F1", "F2", ... and Z_ranked colnames same
if (is.null(colnames(Z_ranked))) stop("Z_ranked has no colnames (expected F1..FK).")

common_f <- intersect(colnames(Z_ranked), class_df$Factor)
if (length(common_f) < 2) {
  stop("Too few overlapping factors between Z and class_df. Check factor naming.")
}

# Restrict to common factors only
Z_use <- Z_ranked[, common_f, drop = FALSE]
class_df_use <- class_df[match(common_f, class_df$Factor), , drop = FALSE]

## ---- 7.3 Pick RNA-specific vs METH-specific factors ----
rna_f  <- class_df_use$Factor[class_df_use$Class %in% c("RNA-specific", "RNA_specific")]
meth_f <- class_df_use$Factor[class_df_use$Class %in% c("METH-specific", "METH_specific")]

if (length(rna_f) == 0 || length(meth_f) == 0) {
  stop("Need at least 1 RNA-specific and 1 METH-specific factor for cross-omics correlation.")
}

Z_rna  <- Z_use[, rna_f,  drop = FALSE]
Z_meth <- Z_use[, meth_f, drop = FALSE]

## ---- 7.4 Correlation matrix (RNA-specific vs METH-specific) ----
cor_mat <- cor(Z_rna, Z_meth, use = "pairwise.complete.obs", method = "pearson")

# Label axes to make the plot self-explanatory
rownames(cor_mat) <- paste0("RNA_",  rownames(cor_mat))
colnames(cor_mat) <- paste0("METH_", colnames(cor_mat))

## ---- 7.5 Plot (MOFA-like diverging palette) ----
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

hm <- pheatmap(
  cor_mat,
  color = colorRampPalette(c("#b2182b", "#f7f7f7", "#2166ac"))(100),
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "ward.D2",
  border_color = NA,
  fontsize = 10,
  main = "Cross-omics correlation: RNA-specific vs METH-specific GFA factors",
  silent = TRUE
)

grid::grid.newpage()
grid::grid.draw(hm$gtable)

## ---- 7.6 Save PNG (fixed path) ----
out_png <- file.path(plots_dir, "GFA_crossomics_corr_RNA_vs_METH_specific.png")
png(out_png, width = 2200, height = 1600, res = 300, bg = "white")
grid::grid.newpage()
grid::grid.draw(hm$gtable)
dev.off()

message("Saved cross-omics correlation heatmap: ", out_png)
message("RNA-specific factors used: ", paste(rna_f, collapse = ", "))
message("METH-specific factors used: ", paste(meth_f, collapse = ", "))

## =========================================================
## GFA.R — STEP 8A: RNA GSEA (GO:BP) across ALL active factors
## Purpose:
##   - Run GSEA on RNA factor weights (ranked, signed weights)
##   - No manual factor selection: uses "active" definition (non-zero signal)
## Inputs:
##   - results/GFA/run_paired01/rds/W_RNA_ranked_F1toFK.rds
## Outputs:
##   - results/GFA/run_paired01/GSEA_RNA/GFA_RNA_GSEA_<FACTOR>.csv
##   - results/GFA/run_paired01/tables/GFA_RNA_GSEA_sig_counts.csv
## =========================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(dplyr)
  library(readr)
})

## ---------------------------
## 8.0 Safety: required PATH vars must exist (from STEP 0A)
## ---------------------------
stopifnot(exists("GFA_RUN_DIR"), exists("RDS_DIR"), exists("TABLES_DIR"))
stopifnot(dir.exists(GFA_RUN_DIR), dir.exists(RDS_DIR), dir.exists(TABLES_DIR))

## ---------------------------
## 8.1 Output dir 
## ---------------------------
GSEA_RNA_DIR <- file.path(GFA_RUN_DIR, "GSEA_RNA")
dir.create(GSEA_RNA_DIR, recursive = TRUE, showWarnings = FALSE)

# Write-permission test (catches OneDrive/permission issues early)
test_file <- file.path(GSEA_RNA_DIR, "TEST_WRITE.txt")
writeLines("test", con = test_file)
stopifnot(file.exists(test_file))

message("GSEA_RNA_DIR ready: ", GSEA_RNA_DIR)

## ---------------------------
## 8.2 Load RNA weights (features x K)
## ---------------------------
W_RNA <- readRDS(file.path(RDS_DIR, "W_RNA_ranked_F1toFK.rds"))
W_RNA <- as.matrix(W_RNA)
storage.mode(W_RNA) <- "double"

stopifnot(!is.null(colnames(W_RNA)))
stopifnot(!is.null(rownames(W_RNA)))  # should be gene symbols

K <- ncol(W_RNA)
message("Loaded W_RNA: ", paste(dim(W_RNA), collapse = " x "), " (K=", K, ")")

## ---------------------------
## 8.3 Active factor detection (no manual selection)
## ---------------------------
eps <- 1e-8
rna_sd  <- apply(W_RNA, 2, sd)
rna_max <- apply(abs(W_RNA), 2, max, na.rm = TRUE)

active_idx     <- which(rna_sd > 0 & rna_max > eps)
active_factors <- colnames(W_RNA)[active_idx]

message("Active RNA factors: ", length(active_factors), " / ", K)
print(active_factors)

## ---------------------------
## 8.4 Helpers: SYMBOL -> ENTREZ + GSEA rank vector
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
  stopifnot(!is.null(names(w_named)))
  ent <- sym2ent(names(w_named))
  
  keep <- !is.na(ent) & !duplicated(ent)
  if (!any(keep)) return(numeric(0))
  
  r <- setNames(as.numeric(w_named[keep]), ent[keep])
  r[order(r, decreasing = TRUE)]
}

## ---------------------------
## 8.5 Run GSEA for each active factor
## ---------------------------
gsea_results_list <- list()

for (f in active_factors) {
  
  w <- W_RNA[, f]
  names(w) <- rownames(W_RNA)
  
  ranks <- prep_gsea_rank(w)
  
  message("=== GSEA RNA ", f, " | ranked genes: ", length(ranks), " ===")
  
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
    
    out_csv <- file.path(GSEA_RNA_DIR, paste0("GFA_RNA_GSEA_", f, ".csv"))
    readr::write_csv(df, out_csv)
    
    gsea_results_list[[f]] <- df
    message("Saved: ", out_csv)
    
  } else {
    message("No significant GO:BP terms for ", f, " at FDR <= 0.05.")
  }
}

## ---------------------------
## 8.6 Summary: significant term counts per factor (FDR <= 0.05)
## ---------------------------
if (length(gsea_results_list) > 0) {
  
  gsea_all <- dplyr::bind_rows(gsea_results_list, .id = "FactorID") |>
    dplyr::mutate(p.adjust = suppressWarnings(as.numeric(p.adjust))) |>
    dplyr::filter(!is.na(p.adjust), p.adjust <= 0.05)
  
  gsea_counts <- gsea_all |>
    dplyr::group_by(Factor = FactorID) |>
    dplyr::summarise(sig_terms = dplyr::n(), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(sig_terms))
  
  out_counts <- file.path(TABLES_DIR, "GFA_RNA_GSEA_sig_counts.csv")
  readr::write_csv(gsea_counts, out_counts)
  
  message("Saved summary counts: ", out_counts)
  
} else {
  message("No RNA GSEA results were produced (no active factors or no significant terms).")
}

message("STEP 8 DONE")

## =========================================================
## GFA.R — STEP 8B: Summarise RNA GSEA outputs (GO:BP)
## Purpose:
##   - Read all per-factor GSEA CSVs from GSEA_RNA_DIR
##   - Summarise number of significant terms (FDR <= 0.05) per factor
## Outputs:
##   - results/GFA/run_paired01/tables/GFA_RNA_GSEA_sig_counts.csv
##   - (optional) barplot saved into plots/
## =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
})

# Load summary table
gfa_counts <- readr::read_csv(
  file.path(TABLES_DIR, "GFA_RNA_GSEA_sig_counts.csv"),
  show_col_types = FALSE
)

# Safety
stopifnot(all(c("Factor", "sig_terms") %in% colnames(gfa_counts)))

# Order factors: top -> bottom (like MOFA)
gfa_counts <- gfa_counts %>%
  arrange(desc(sig_terms)) %>%
  mutate(Factor = factor(Factor, levels = rev(Factor)))

# Plot (MOFA-style)
p <- ggplot(gfa_counts, aes(x = sig_terms, y = Factor, fill = sig_terms)) +
  geom_col(width = 0.8) +
  scale_fill_gradient(
    low  = "#cfe1f2",
    high = "#0b3c5d",
    name = "GO:BP term count"
  ) +
  labs(
    title = "GFA RNA factors (GSEA, FDR ≤ 0.05): significant GO:BP terms",
    x = "Number of enriched GO:BP terms",
    y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

print(p)

# Save
out_png <- file.path(PLOTS_DIR, "GFA_RNA_GSEA_sig_terms_horizontal.png")
ggsave(out_png, p, width = 10, height = 5, dpi = 300, bg = "white")

message("Saved: ", out_png)

## =========================================================
## GFA.R — STEP 9: METH ORA (GO:BP) across ALL active factors
## Purpose:
##   - Run ORA on METH factor weights using top-|weights| genes
##   - No manual factor selection: uses "active" definition
## Input (RDS_DIR):
##   - W_METH_ranked_F1toFK.rds  (genes x K)
## Output (GFA_RUN_DIR):
##   - ORA_METH/GFA_METH_ORA_<FACTOR>.csv
##   - tables/GFA_METH_ORA_sig_counts.csv
## =========================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(dplyr)
  library(readr)
})

## ---------------------------
## 9.1 Fixed output paths (DO NOT change your PATH block)
## ---------------------------

ORA_METH_DIR <- file.path(GFA_RUN_DIR, "ORA_METH")
dir.create(ORA_METH_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLES_DIR,  recursive = TRUE, showWarnings = FALSE)

## ---------------------------
## 9.2 Load METH weights (genes x K)
## ---------------------------

W_METH <- readRDS(file.path(RDS_DIR, "W_METH_ranked_F1toFK.rds"))
W_METH <- as.matrix(W_METH)
storage.mode(W_METH) <- "double"

stopifnot(!is.null(colnames(W_METH)))
stopifnot(!is.null(rownames(W_METH)))  # should be gene symbols

K <- ncol(W_METH)
message("Loaded W_METH: ", paste(dim(W_METH), collapse = " x "), " (K=", K, ")")

## ---------------------------
## 9.3 Active factor detection (no manual selection)
## ---------------------------

eps <- 1e-8
meth_sd  <- apply(W_METH, 2, sd)
meth_max <- apply(abs(W_METH), 2, max, na.rm = TRUE)

active_idx <- which(meth_sd > 0 & meth_max > eps)
active_factors <- colnames(W_METH)[active_idx]

message("Active METH factors: ", length(active_factors), " / ", K)
print(active_factors)

## ---------------------------
## 9.4 Helpers: SYMBOL -> ENTREZ, build ORA gene set
## ---------------------------

sym2ent <- function(sym) {
  sym <- unique(na.omit(as.character(sym)))
  if (length(sym) == 0L) return(character(0))
  
  sym_up <- toupper(sym)
  all_keys <- AnnotationDbi::keys(org.Hs.eg.db, keytype = "SYMBOL")
  sym_ok <- intersect(sym_up, all_keys)
  
  if (length(sym_ok) == 0L) return(rep(NA_character_, length(sym)))
  
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

prep_ora_genes <- function(w_named, topN = 200) {
  stopifnot(!is.null(names(w_named)))
  
  # Select top-|weights| genes
  ord <- order(abs(w_named), decreasing = TRUE)
  top_sym <- names(w_named)[ord][seq_len(min(topN, length(ord)))]
  
  ent <- sym2ent(top_sym)
  ent <- unique(na.omit(as.character(ent)))
  
  ent
}

## ---------------------------
## 9.5 Run ORA per active factor
## ---------------------------

TOPN_ORA <- 200
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
    
    out_csv <- file.path(ORA_METH_DIR, paste0("GFA_METH_ORA_", f, ".csv"))
    readr::write_csv(df, out_csv)
    
    ora_results_list[[f]] <- df
    message("Saved: ", out_csv)
  } else {
    message("No significant GO:BP terms for ", f, " at FDR <= 0.05.")
  }
}

## ---------------------------
## 9.6 Summary: significant term counts per factor (FDR <= 0.05)
## ---------------------------

if (length(ora_results_list) > 0) {
  
  ora_all <- dplyr::bind_rows(ora_results_list, .id = "FactorID") |>
    dplyr::mutate(p.adjust = suppressWarnings(as.numeric(p.adjust))) |>
    dplyr::filter(!is.na(p.adjust), p.adjust <= 0.05)
  
  ora_counts <- ora_all |>
    dplyr::group_by(Factor = FactorID) |>
    dplyr::summarise(sig_terms = dplyr::n(), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(sig_terms))
  
  out_counts <- file.path(TABLES_DIR, "GFA_METH_ORA_sig_counts.csv")
  readr::write_csv(ora_counts, out_counts)
  
  message("Saved summary counts: ", out_counts)
  
} else {
  message("No METH ORA results were produced (no active factors or no significant terms).")
}

message("STEP 9 DONE")

## =========================================================
## GFA METH factors (ORA, FDR ≤ 0.05): significant GO:BP terms
## MOFA-style horizontal barplot
## =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(forcats)
})

# Load ORA summary table
meth_ora_counts <- readr::read_csv(
  file.path(TABLES_DIR, "GFA_METH_ORA_sig_counts.csv"),
  show_col_types = FALSE
)

# Safety check
stopifnot(all(c("Factor", "sig_terms") %in% colnames(meth_ora_counts)))

# Order factors by decreasing number of enriched terms
meth_ora_counts <- meth_ora_counts %>%
  arrange(desc(sig_terms)) %>%
  mutate(Factor = factor(Factor, levels = rev(Factor)))

# Plot (MOFA-style)
p_meth_ora <- ggplot(
  meth_ora_counts,
  aes(x = sig_terms, y = Factor, fill = sig_terms)
) +
  geom_col(width = 0.8) +
  scale_fill_gradient(
    low  = "#cfe1f2",
    high = "#0b3c5d",
    name = "GO:BP term count"
  ) +
  labs(
    title = "GFA methylation factors (ORA, FDR ≤ 0.05): significant GO:BP terms",
    x = "Number of enriched GO:BP terms",
    y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

print(p_meth_ora)

# Save plot
out_png <- file.path(PLOTS_DIR, "GFA_METH_ORA_sig_terms_horizontal.png")
ggsave(
  filename = out_png,
  plot     = p_meth_ora,
  width    = 10,
  height   = 5,
  dpi      = 300,
  bg       = "white"
)

message("Saved ORA plot: ", out_png)
