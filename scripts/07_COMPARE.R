
## =============================================================
## 01_compare_GOBP_GSEA.R — COMPARE: GO:BP GSEA (MOFA vs GFA)
## Purpose:
##   - Read MOFA RNA GSEA (GO:BP) + GFA RNA GSEA (GO:BP)
##   - Make pathway plots for comparison
## Outputs (ONLY):
##   - results/COMPARE/plots/*
##   - results/COMPARE/tables/*
##   - results/COMPARE/rds/*
## =============================================================

source("scripts/00_setup_paths.R")


## ---------------------------
## 0) COMPARE output folders (already exist in your screenshot)
## ---------------------------
compare_plots_dir  <- file.path(compare_dir, "plots")
compare_tables_dir <- file.path(compare_dir, "tables")
compare_rds_dir    <- file.path(compare_dir, "rds")

dir.create(compare_plots_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(compare_tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(compare_rds_dir,    recursive = TRUE, showWarnings = FALSE)

## ---------------------------
## 1) Helpers
## ---------------------------
std_gsea_cols <- function(df) {
  # Normalize column names used by clusterProfiler outputs
  # Accepts common variants: p.adjust / p.adjust., etc.
  cn <- colnames(df)
  
  # Try to keep original, but ensure these exist: ID, Description, NES, p.adjust, Factor
  # clusterProfiler typical: ID, Description, setSize, enrichmentScore, NES, pvalue, p.adjust, qvalues, rank, leading_edge, core_enrichment
  
  if (!"p.adjust" %in% cn) {
    pa <- cn[grepl("^p\\.adjust$", cn)]
    if (length(pa) == 1) colnames(df)[match(pa, colnames(df))] <- "p.adjust"
  }
  if (!"Description" %in% cn && "description" %in% cn) {
    colnames(df)[match("description", colnames(df))] <- "Description"
  }
  if (!"NES" %in% cn && "nes" %in% cn) {
    colnames(df)[match("nes", colnames(df))] <- "NES"
  }
  df
}

extract_factor_from_filename <- function(x) {
  f <- basename(x)
  f <- gsub("\\.csv$", "", f)
  
  # MOFA_RNA_GSEA_Factor12.csv -> F12
  # MOFA_RNA_GSEA_F12.csv     -> F12
  f <- sub(".*Factor", "F", f)
  f <- sub(".*_F", "F", f)
  
  f
}


read_gsea_dir <- function(dir_path, method_label) {
  
  all_file <- file.path(dir_path, paste0(method_label, "_GSEA_RNA_all.csv"))
  
  if (file.exists(all_file)) {
    df <- readr::read_csv(all_file, show_col_types = FALSE) %>% std_gsea_cols()
    if (!"Factor" %in% colnames(df)) {
      stop("Found ", basename(all_file), " but it has no 'Factor' column.")
    }
    df <- df %>% mutate(Method = method_label)
    return(df)
  }
  
  # fallback: per-factor files
  patt <- if (method_label == "MOFA") {
    "^MOFA_RNA_GSEA_(Factor|F)\\d+\\.csv$"
  } else {
    "^gfa_rna_gsea_f\\d+\\.csv$"
  }
  
  files <- list.files(dir_path, pattern = patt, full.names = TRUE)
  if (length(files) == 0) {
    stop("No GSEA files found in: ", dir_path, " (expected ", patt, " or *_GSEA_RNA_all.csv)")
  }
  
  df <- purrr::map_dfr(files, function(f) {
    readr::read_csv(f, show_col_types = FALSE) %>%
      std_gsea_cols() %>%
      mutate(Factor = extract_factor_from_filename(f))
  }) %>%
    mutate(Method = method_label)
  
  df
}

top_terms_per_factor <- function(df, padj_cutoff = 0.05, top_n = 8) {
  df %>%
    filter(!is.na(NES), !is.na(p.adjust)) %>%
    filter(p.adjust <= padj_cutoff) %>%
    group_by(Method, Factor) %>%
    arrange(desc(abs(NES)), .by_group = TRUE) %>%
    slice_head(n = top_n) %>%
    ungroup()
}

jaccard <- function(a, b) {
  a <- unique(a); b <- unique(b)
  if (length(a) == 0 && length(b) == 0) return(NA_real_)
  inter <- length(intersect(a, b))
  uni   <- length(union(a, b))
  if (uni == 0) return(NA_real_)
  inter / uni
}

## ---------------------------
## 2) Load MOFA + GFA GO:BP GSEA
## ---------------------------
mofa_gsea_dir <- file.path(mofa_dir, "GSEA", "tables")
gfa_gsea_dir  <- file.path(gfa_dir,  "GSEA", "tables")

mofa_gsea <- read_gsea_dir(mofa_gsea_dir, "MOFA")
gfa_gsea  <- read_gsea_dir(gfa_gsea_dir,  "GFA")

gsea_all <- bind_rows(mofa_gsea, gfa_gsea)

# Minimal sanity checks
stopifnot(all(c("Method", "Factor", "Description", "NES", "p.adjust") %in% colnames(gsea_all)))

## Save combined table (raw-ish)
readr::write_csv(gsea_all, file.path(compare_tables_dir, "COMPARE_GSEA_GO_BP_all.csv"))

## =========================================================
## Plot A (NEW) — Global Top GO:BP terms (NO factor split)
##   - Collapse across factors
##   - Pick best (lowest p.adjust) record per term
##   - Plot top 20 terms per method
## =========================================================

padj_cutoff <- 0.05
top_n_terms <- 20

# --- split data
gsea_mofa <- gsea_global_top %>% filter(Method == "MOFA")
gsea_gfa  <- gsea_global_top %>% filter(Method == "GFA")

# --- MOFA plot
p_mofa <- ggplot(gsea_mofa, aes(x = NES, y = reorder(Term, NES), fill = NES)) +
  geom_col(width = 0.85) +
  scale_fill_gradient(low = "#c6dbef", high = "#0B3C5D") +
  labs(
    title = "MOFA2 RNA — GSEA (GO:BP): Top 20 enriched terms",
    subtitle = paste0("Collapsed across factors; p.adjust ≤ ", padj_cutoff),
    x = "Normalized Enrichment Score (NES)",
    y = NULL,
    fill = "NES"
  ) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 9))
p_mofa
ggsave(
  filename = file.path(compare_plots_dir, "MOFA_GSEA_GO_BP_global_top20.png"),
  plot = p_mofa,
  width = 10, height = 6.5, dpi = 300
)

# --- GFA plot
p_gfa <- ggplot(gsea_gfa, aes(x = NES, y = reorder(Term, NES), fill = NES)) +
  geom_col(width = 0.85) +
  scale_fill_gradient(low = "#c6dbef", high = "#0B3C5D") +
  labs(
    title = "GFA RNA — GSEA (GO:BP): Top 20 enriched terms",
    subtitle = paste0("Collapsed across factors; p.adjust ≤ ", padj_cutoff),
    x = "Normalized Enrichment Score (NES)",
    y = NULL,
    fill = "NES"
  ) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 9))
p_gfa
ggsave(
  filename = file.path(compare_plots_dir, "GFA_GSEA_GO_BP_global_top20.png"),
  plot = p_gfa,
  width = 10, height = 6.5, dpi = 300
)


## =============================================================
## 02_compare_ORA_GOBP_METH.R — COMPARE: METH ORA (GO:BP) MOFA vs GFA
## Purpose:
##   - Read MOFA METH ORA tables + GFA METH ORA tables
##   - Collapse across factors (no factor plots)
##   - Save top 20 terms per method + barplots (separate PNGs)
## Outputs (ONLY):
##   - results/COMPARE/tables/*
##   - results/COMPARE/plots/*
## =============================================================


## ---------------------------
## 0) COMPARE output folders
## ---------------------------
compare_plots_dir  <- file.path(compare_dir, "plots")
compare_tables_dir <- file.path(compare_dir, "tables")
dir.create(compare_plots_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(compare_tables_dir, recursive = TRUE, showWarnings = FALSE)

## ---------------------------
## 1) Helpers
## ---------------------------
extract_factor_from_filename <- function(x) {
  f <- basename(x)
  f <- gsub("\\.csv$", "", f)
  f <- sub(".*Factor", "F", f)  # MOFA...Factor12 -> F12
  f <- sub(".*_f", "F", f)      # gfa..._f12      -> F12
  f <- sub(".*_F", "F", f)      # ..._F12         -> F12 (fallback)
  f
}

std_ora_cols <- function(df) {
  # Make column names consistent across MOFA/GFA outputs
  cn <- colnames(df)
  
  # Common ORA columns you might have: ID/Description/p.adjust/pvalue/Count/GeneRatio/etc.
  # Enrichr-like: Term, Adjusted.P.value, Overlap, P.value, Combined.Score...
  # We'll standardize minimal needed: Description + p.adjust + Count (if exists)
  
  # Description / Term
  if (!"Description" %in% cn) {
    if ("Term" %in% cn) colnames(df)[match("Term", cn)] <- "Description"
    if ("term" %in% cn) colnames(df)[match("term", cn)] <- "Description"
  }
  
  # p.adjust
  if (!"p.adjust" %in% cn) {
    cand <- cn[cn %in% c("Adjusted.P.value", "adj_p", "FDR", "padj", "p_adj")]
    if (length(cand) >= 1) colnames(df)[match(cand[1], cn)] <- "p.adjust"
  }
  
  # Count (optional)
  if (!"Count" %in% cn) {
    cand2 <- cn[cn %in% c("Overlap", "overlap", "Genes", "genes")]
    # do nothing; Count is optional for plotting
  }
  
  df
}

read_ora_dir <- function(dir_path, method_label) {
  
  # We read per-factor tables only (as in your folders)
  patt <- if (method_label == "MOFA") {
    "^MOFA_METH_ORA_Factor\\d+\\.csv$"
  } else {
    "^gfa_meth_ora_f\\d+\\.csv$"
  }
  
  files <- list.files(dir_path, pattern = patt, full.names = TRUE)
  if (length(files) == 0) {
    stop("No ORA files found in: ", dir_path, " (expected ", patt, ")")
  }
  
  df <- purrr::map_dfr(files, function(f) {
    readr::read_csv(f, show_col_types = FALSE) %>%
      std_ora_cols() %>%
      mutate(Factor = extract_factor_from_filename(f))
  }) %>%
    mutate(Method = method_label)
  
  df
}

## ---------------------------
## 2) Load MOFA + GFA ORA tables
## ---------------------------
mofa_ora_dir <- file.path(mofa_dir, "ORA", "tables")
gfa_ora_dir  <- file.path(gfa_dir,  "ORA", "tables")

mofa_ora <- read_ora_dir(mofa_ora_dir, "MOFA")
gfa_ora  <- read_ora_dir(gfa_ora_dir,  "GFA")

ora_all <- bind_rows(mofa_ora, gfa_ora)

# Minimal sanity checks
stopifnot(all(c("Method", "Factor", "Description") %in% colnames(ora_all)))
if (!"p.adjust" %in% colnames(ora_all)) {
  stop("ORA tables do not contain an adjusted p-value column. Expected 'p.adjust' (or mappable variant).")
}

readr::write_csv(ora_all, file.path(compare_tables_dir, "COMPARE_ORA_GO_BP_all.csv"))

## ---------------------------
## 3) Global collapse across factors (NO factor split)
## ---------------------------
padj_cutoff <- 0.05
top_n_terms <- 20

ora_global <- ora_all %>%
  filter(!is.na(Description), !is.na(p.adjust)) %>%
  filter(p.adjust <= padj_cutoff) %>%
  dplyr::select(where(~ !is.list(.x))) %>%
  group_by(Method, Description) %>%
  arrange(p.adjust, .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup()

readr::write_csv(
  ora_global,
  file.path(compare_tables_dir, "COMPARE_ORA_GO_BP_global_collapsed.csv")
)

ora_global_top <- ora_global %>%
  group_by(Method) %>%
  arrange(p.adjust, .by_group = TRUE) %>%
  slice_head(n = top_n_terms) %>%
  ungroup() %>%
  mutate(Term = stringr::str_wrap(Description, width = 55),
         neglog10_padj = -log10(p.adjust))

readr::write_csv(
  ora_global_top,
  file.path(compare_tables_dir, "COMPARE_ORA_GO_BP_global_top20.csv")
)

## ---------------------------
## 4) Plots (separate PNGs)
## ---------------------------
ora_mofa <- ora_global_top %>% filter(Method == "MOFA")
ora_gfa  <- ora_global_top %>% filter(Method == "GFA")

p_mofa <- ggplot(ora_mofa, aes(x = neglog10_padj, y = reorder(Term, neglog10_padj))) +
  geom_col(fill = "#0B3C5D", width = 0.85) +
  labs(
    title = "MOFA2 METH — ORA (GO:BP): Top 20 enriched terms",
    subtitle = paste0("Collapsed across factors; p.adjust ≤ ", padj_cutoff),
    x = "-log10(p.adjust)",
    y = NULL
  ) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 9))
p_mofa
ggsave(
  filename = file.path(compare_plots_dir, "MOFA_ORA_GO_BP_global_top20.png"),
  plot = p_mofa,
  width = 10, height = 6.5, dpi = 300
)

p_gfa <- ggplot(ora_gfa, aes(x = neglog10_padj, y = reorder(Term, neglog10_padj))) +
  geom_col(fill = "#0B3C5D", width = 0.85) +
  labs(
    title = "GFA METH — ORA (GO:BP): Top 20 enriched terms",
    subtitle = paste0("Collapsed across factors; p.adjust ≤ ", padj_cutoff),
    x = "-log10(p.adjust)",
    y = NULL
  ) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 9))
p_gfa
ggsave(
  filename = file.path(compare_plots_dir, "GFA_ORA_GO_BP_global_top20.png"),
  plot = p_gfa,
  width = 10, height = 6.5, dpi = 300
)

message("DONE: ORA compare tables+plots written to results/COMPARE/")



## =============================================================
## UpSet (MOFA vs GFA) — GO:BP overlap (RNA GSEA + METH ORA)
## Outputs ONLY:
##   results/COMPARE/plots/UPSET_GSEA_GO_BP_MOFA_vs_GFA.png
##   results/COMPARE/plots/UPSET_ORA_GO_BP_MOFA_vs_GFA.png
## =============================================================

compare_plots_dir  <- file.path(compare_dir, "plots")
compare_tables_dir <- file.path(compare_dir, "tables")
dir.create(compare_plots_dir, recursive = TRUE, showWarnings = FALSE)

padj_cutoff <- 0.05

## ---------------------------
## 1) RNA — GSEA UpSet
## ---------------------------
gsea_file <- file.path(compare_tables_dir, "COMPARE_GSEA_GO_BP_global_collapsed.csv")
gsea_df <- readr::read_csv(gsea_file, show_col_types = FALSE) %>%
  dplyr::filter(!is.na(Description), !is.na(p.adjust), p.adjust <= padj_cutoff)

mofa_terms_rna <- unique(gsea_df %>% dplyr::filter(Method == "MOFA") %>% dplyr::pull(Description))
gfa_terms_rna  <- unique(gsea_df %>% dplyr::filter(Method == "GFA")  %>% dplyr::pull(Description))

upset_input_rna <- UpSetR::fromList(list(
  GFA_RNA  = gfa_terms_rna,
  MOFA_RNA = mofa_terms_rna
))

out_png_rna <- file.path(compare_plots_dir, "UPSET_GSEA_GO_BP_MOFA_vs_GFA.png")

png(filename = out_png_rna, width = 1400, height = 900, res = 150)

UpSetR::upset(
  upset_input_rna,
  sets = c("GFA_RNA", "MOFA_RNA"),
  order.by = "freq",
  decreasing = FALSE,
  keep.order = TRUE,
  mainbar.y.label = "Intersection size",
  sets.x.label = "Set size",
  main.bar.color = "#0B3C5D",
  sets.bar.color = "#0B3C5D",
  matrix.color   = "#0B3C5D",
  matrix.dot.alpha = 0.9,
  text.scale = c(1.6, 1.6, 1.2, 1.2, 1.4, 1.2),
  mainbar.y.max = NA
)

# Title should be under the main bar plot (like the grey example)
grid::grid.text(
  "GSEA – Overlap of GO:BP terms between MOFA2 and GFA (RNA)",
  x = 0.6, y = 0.29,
  gp = grid::gpar(fontsize = 11, fontface = "bold")
)


# Bottom label (like the grey example)
grid::grid.text(
  "Number of enriched GO:BP terms",
  x = 0.62, y = 0.06,
  gp = grid::gpar(fontsize = 11)
)

dev.off()

## ---------------------------
## METH — ORA UpSet (FULL, with labels + title inside plot area)
## ---------------------------

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(UpSetR)
  library(grid)
})

padj_cutoff <- 0.05

ora_file <- file.path(compare_tables_dir, "COMPARE_ORA_GO_BP_global_collapsed.csv")
stopifnot(file.exists(ora_file))

ora_df <- readr::read_csv(ora_file, show_col_types = FALSE) %>%
  dplyr::filter(!is.na(Description), !is.na(p.adjust), p.adjust <= padj_cutoff)

mofa_terms_meth <- unique(ora_df %>% dplyr::filter(Method == "MOFA") %>% dplyr::pull(Description))
gfa_terms_meth  <- unique(ora_df %>% dplyr::filter(Method == "GFA")  %>% dplyr::pull(Description))

upset_input_meth <- UpSetR::fromList(list(
  MOFA_METH = mofa_terms_meth,
  GFA_METH  = gfa_terms_meth
))

out_png_meth <- file.path(compare_plots_dir, "UPSET_ORA_GO_BP_MOFA_vs_GFA.png")

# counts to decide decreasing so SHARED ends up last as often as possible
dfu <- upset_input_meth
n_shared    <- sum(dfu$MOFA_METH == 1 & dfu$GFA_METH == 1)
n_mofa_only <- sum(dfu$MOFA_METH == 1 & dfu$GFA_METH == 0)
n_gfa_only  <- sum(dfu$MOFA_METH == 0 & dfu$GFA_METH == 1)

decreasing_flag <- if (n_shared == min(c(n_shared, n_mofa_only, n_gfa_only))) TRUE else FALSE
if (n_shared == max(c(n_shared, n_mofa_only, n_gfa_only))) decreasing_flag <- FALSE

png(filename = out_png_meth, width = 1400, height = 900, res = 150)

UpSetR::upset(
  upset_input_meth,
  sets = c("MOFA_METH", "GFA_METH"),
  order.by = "freq",
  decreasing = decreasing_flag,
  keep.order = TRUE,
  
  mainbar.y.label = "Intersection size",
  sets.x.label    = "Set size",
  
  main.bar.color = "#0B3C5D",
  sets.bar.color = "#0B3C5D",
  matrix.color   = "#0B3C5D",
  matrix.dot.alpha = 0.9,
  
  text.scale = c(1.6, 1.6, 1.2, 1.2, 1.4, 1.2),
  mainbar.y.max = NA
)

# title inside plot area (same placement style as your RNA)
grid::grid.text(
  "ORA – Overlap of GO:BP terms between MOFA2 and GFA (Methylation)",
  x = 0.6, y = 0.29,
  gp = grid::gpar(fontsize = 11, fontface = "bold")
)

# bottom label (like the grey example)
grid::grid.text(
  "Number of enriched GO:BP terms",
  x = 0.5, y = 0.06,
  gp = grid::gpar(fontsize = 12)
)

dev.off()
message("Saved: ", out_png_meth)
