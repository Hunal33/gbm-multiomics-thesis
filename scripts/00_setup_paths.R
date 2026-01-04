# ==============================================================
# 00_setup_paths.R
# Purpose:
#   Define ALL project paths (data + scripts + results)
#   Create output folders safely
#   Provide mode-switch paths for transcriptomics / methylation runs
# Author: Havva Ünal
# ==============================================================

## -----------------------------
## 0) PROJECT ROOT (DEFINE ONCE)
## -----------------------------
ROOT <- "C:/Users/arcal/OneDrive/Masaüstü/gbm-multiomics-thesis"

## -----------------------------
## 1) TOP-LEVEL FOLDERS
## -----------------------------
DATA_DIR    <- file.path(ROOT, "data")
SCRIPTS_DIR <- file.path(ROOT, "scripts")
RESULTS_DIR <- file.path(ROOT, "results")
FIGURES_DIR <- file.path(ROOT, "figures")

## -----------------------------
## 2) DATA FOLDERS (INPUTS)
## -----------------------------
# Transcriptomics
tx_dir      <- file.path(DATA_DIR, "transcriptomics")
tx_unpaired <- file.path(tx_dir, "unpaired")
tx_paired   <- file.path(tx_dir, "paired")

# Methylation
meth_dir      <- file.path(DATA_DIR, "methylation")
meth_unpaired <- file.path(meth_dir, "unpaired")
meth_paired   <- file.path(meth_dir, "paired")

## -----------------------------
## 3) RESULTS OUTPUT FOLDERS
## -----------------------------
tx_run_unpaired   <- file.path(RESULTS_DIR, "transcriptomics_unpaired")
tx_run_paired     <- file.path(RESULTS_DIR, "transcriptomics_paired")

meth_run_unpaired <- file.path(RESULTS_DIR, "methylation_unpaired")
meth_run_paired   <- file.path(RESULTS_DIR, "methylation_paired")

paired_inputs_dir <- file.path(RESULTS_DIR, "paired_inputs")

mofa_dir    <- file.path(RESULTS_DIR, "MOFA")
gfa_dir     <- file.path(RESULTS_DIR, "GFA")
compare_dir <- file.path(RESULTS_DIR, "COMPARE")  # matches your current structure

## -----------------------------
## 4) MODE SWITCHES (ONE LINE CONTROL)
## -----------------------------
# Choose which transcriptomics inputs/outputs your DE script should use
TX_MODE   <- "unpaired"   # "unpaired" or "paired"
METH_MODE <- "paired"     # "unpaired" or "paired" (set when needed)

# Resolve mode-dependent paths
tx_countsdir <- if (TX_MODE == "paired") tx_paired else tx_unpaired
tx_run_dir   <- if (TX_MODE == "paired") tx_run_paired else tx_run_unpaired

meth_input_dir <- if (METH_MODE == "paired") meth_paired else meth_unpaired
meth_run_dir   <- if (METH_MODE == "paired") meth_run_paired else meth_run_unpaired

## -----------------------------
## 5) CREATE OUTPUT FOLDERS (SAFE)
## -----------------------------
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGURES_DIR, showWarnings = FALSE, recursive = TRUE)

dir.create(tx_run_unpaired, showWarnings = FALSE, recursive = TRUE)
dir.create(tx_run_paired,   showWarnings = FALSE, recursive = TRUE)

dir.create(meth_run_unpaired, showWarnings = FALSE, recursive = TRUE)
dir.create(meth_run_paired,   showWarnings = FALSE, recursive = TRUE)

dir.create(paired_inputs_dir, showWarnings = FALSE, recursive = TRUE)

dir.create(mofa_dir,    showWarnings = FALSE, recursive = TRUE)
dir.create(gfa_dir,     showWarnings = FALSE, recursive = TRUE)
dir.create(compare_dir, showWarnings = FALSE, recursive = TRUE)

## -----------------------------
## 6) BACKWARD-COMPATIBILITY ALIASES
## -----------------------------
project_dir <- ROOT
data_dir    <- DATA_DIR
scripts_dir <- SCRIPTS_DIR
results_dir <- RESULTS_DIR
figures_dir <- FIGURES_DIR

## -----------------------------
## 7) DEBUG PRINT (QUICK CHECK)
## -----------------------------
message("ROOT           : ", ROOT)
message("DATA_DIR       : ", DATA_DIR)
message("SCRIPTS_DIR    : ", SCRIPTS_DIR)
message("RESULTS_DIR    : ", RESULTS_DIR)
message("FIGURES_DIR    : ", FIGURES_DIR)

message("TX_MODE        : ", TX_MODE)
message("tx_countsdir   : ", tx_countsdir)
message("tx_run_dir     : ", tx_run_dir)

message("METH_MODE      : ", METH_MODE)
message("meth_input_dir : ", meth_input_dir)
message("meth_run_dir   : ", meth_run_dir)

message("paired_inputs  : ", paired_inputs_dir)
message("MOFA dir       : ", mofa_dir)
message("GFA dir        : ", gfa_dir)
message("COMPARE dir    : ", compare_dir)
