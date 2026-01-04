# gbm-multiomics-thesis
Systematic benchmarking of MOFA2 and GFA using paired RNA-seq and promoter-level DNA methylation data from TCGA-GBM.
## Repository structure

scripts/
├── 00_setup_paths.R        # Central path configuration
├── paired/                 # Paired RNA–methylation preprocessing
├── unpaired/               # Unpaired RNA-seq preprocessing
├── methylation/            # DNA methylation processing (CpG → promoter)
├── MOFA/                   # MOFA2 model fitting and downstream analyses
├── GFA/                    # GFA model fitting and downstream analyses
├── COMPARE/                # Cross-model comparison (MOFA2 vs GFA)
├── intersection/           # Functional overlap and intersection analyses

rds/
├── Intermediate R objects saved during analysis

tables/
├── Result tables (GSEA / ORA / comparison outputs)
