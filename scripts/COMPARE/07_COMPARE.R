## =============================================================
## 07_COMPARE.R — GO:BP compare assembly (MOFA2 vs GFA)
## Purpose:
##   - Read MOFA2 and GFA enrichment outputs (GO:BP)
##   - Standardize into compare-ready tables (RNA: GSEA, METH: ORA)
##   - Produce global Top-20 plots and UpSet overlaps
## Output (ONLY under results/COMPARE/run_paired01):
##   - tables/*.csv
##   - rds/*.rds
##   - plots/*.png
## =============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(forcats)
})

## =============================================================
## PATHS (project-safe; ONLY from 00_setup_paths.R)
## =============================================================
source("scripts/00_setup_paths.R")

# Canonicalize (use setup aliases if needed)
if (!exists("ROOT") && exists("project_dir")) ROOT <- project_dir
if (!exists("RESULTS_DIR") && exists("results_dir")) RESULTS_DIR <- results_dir
if (!exists("RESULTS_DIR")) stop("RESULTS_DIR/results_dir not found. Check scripts/00_setup_paths.R")

RUN_ID <- "run_paired01"

## ---- Compare output root (setup-defined) ----
if (!exists("compare_dir")) stop("compare_dir not found. Check scripts/00_setup_paths.R")
COMPARE_DIR <- file.path(compare_dir, RUN_ID)

TABLES_DIR  <- file.path(COMPARE_DIR, "tables")
RDS_DIR     <- file.path(COMPARE_DIR, "rds")
PLOTS_DIR   <- file.path(COMPARE_DIR, "plots")

dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RDS_DIR,    recursive = TRUE, showWarnings = FALSE)
dir.create(PLOTS_DIR,  recursive = TRUE, showWarnings = FALSE)

## ---- MOFA run paths (setup-defined) ----
if (!exists("mofa_dir")) stop("mofa_dir not found. Check scripts/00_setup_paths.R")
MOFA_RUN_DIR <- file.path(mofa_dir, RUN_ID)
MOFA_GSEA_RNA_DIR  <- file.path(MOFA_RUN_DIR, "GSEA_RNA")
MOFA_ORA_METH_DIR  <- file.path(MOFA_RUN_DIR, "ORA_METH")

## Keep the variable name used later in your script
IN_FILE <- file.path(MOFA_GSEA_RNA_DIR, "MOFA_RNA_GSEA_ALL.csv")

## ---- GFA run paths (setup-defined) ----
if (!exists("gfa_dir")) stop("gfa_dir not found. Check scripts/00_setup_paths.R")
GFA_RUN_DIR <- file.path(gfa_dir, RUN_ID)
GFA_GSEA_RNA_DIR <- file.path(GFA_RUN_DIR, "GSEA_RNA")
GFA_ORA_METH_DIR <- file.path(GFA_RUN_DIR, "ORA_METH")

message("ROOT:            ", ROOT)
message("RESULTS_DIR:     ", RESULTS_DIR)
message("COMPARE_DIR:     ", COMPARE_DIR)
message("TABLES_DIR:      ", TABLES_DIR)
message("RDS_DIR:         ", RDS_DIR)
message("PLOTS_DIR:       ", PLOTS_DIR)

message("MOFA_RUN_DIR:    ", MOFA_RUN_DIR)
message("MOFA_GSEA_RNA:   ", MOFA_GSEA_RNA_DIR)
message("MOFA_ORA_METH:   ", MOFA_ORA_METH_DIR)
message("MOFA IN_FILE:    ", IN_FILE)

message("GFA_RUN_DIR:     ", GFA_RUN_DIR)
message("GFA_GSEA_RNA:    ", GFA_GSEA_RNA_DIR)
message("GFA_ORA_METH:    ", GFA_ORA_METH_DIR)

stopifnot(dir.exists(RESULTS_DIR))
stopifnot(dir.exists(MOFA_RUN_DIR))
stopifnot(dir.exists(MOFA_GSEA_RNA_DIR))
stopifnot(dir.exists(GFA_RUN_DIR))
stopifnot(dir.exists(GFA_GSEA_RNA_DIR))
stopifnot(file.exists(IN_FILE))

## =============================================================
## (REST OF YOUR SCRIPT CONTINUES BELOW — UNCHANGED)
## =============================================================
