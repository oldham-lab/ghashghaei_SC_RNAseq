setwd("/mnt/bdata/rebecca/collabs/ghashghaei/SC_expr/P5/analyses/project_modules")

source("/mnt/bdata/rebecca/collabs/ghashghaei/SC_expr/project_modules_fxn.R")

expr <- readRDS("../../cluster/expr_SC_P5.RDS")
fm_dir <- "/mnt/bdata/rebecca/collabs/ghashghaei/bulk_expr/P5/analyses/FM/featureCounts_P5_Modules"
projectname <- "P5_SC_RNAseq"
expr_type <- "normalized_counts"
n_genes <- 20
pval_cut <- 1e-10

expr$Cluster <- expr$customclassif

project_modules(projectname, expr, fm_dir, expr_type, pval_cut, n_genes)