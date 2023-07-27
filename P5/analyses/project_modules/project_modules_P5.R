.libPaths("/home/shared/R/x86_64-pc-linux-gnu-library/4.0")

setwd("/mnt/bdata/rebecca/collabs/ghashghaei/SC_expr/P5/analyses/project_modules")

source("/mnt/bdata/rebecca/collabs/ghashghaei/SC_expr/project_modules_fxns.R")

expr <- readRDS("../../cluster/expr_SC_P5.RDS")
fm_dir <- "/mnt/bdata/rebecca/collabs/ghashghaei/bulk_expr/P5/analyses/FM/featureCounts_P5_Modules"
projectname <- "P5_SC_RNAseq"
expr_type <- "normalized_counts"
pval_cut <- 1e-20
n_genes <- 100

expr$Cluster <- expr$customclassif

project_modules_rel_group(projectname, expr, fm_dir, expr_type, pval_cut, n_genes)
