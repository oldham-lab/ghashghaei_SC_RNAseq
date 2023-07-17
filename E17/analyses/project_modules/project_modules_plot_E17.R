setwd("/mnt/bdata/rebecca/collabs/ghashghaei/SC_expr/E17/analyses/project_modules")

source("/mnt/bdata/rebecca/collabs/ghashghaei/SC_expr/project_modules_plot_fxn.R")

df <- fread("data/E17_SC_RNAseq_normalized_counts_module_projections_pval_cut_1e-10_top_20_module_genes.csv", data.table=F)

expr <- readRDS("../../cluster/expr_SC_E17.RDS")

cellinfo <- expr[[]]
cellinfo$Background <- "Homo"
cellinfo$Background[grep("F+", cellinfo$Genotype, fixed=T)] <- "Het"
cellinfo$EGFR_Status <- cellinfo$Condition

projectname <- "E17_SC_RNAseq"
expr_type <- "normalized_counts"
pval_cut <- 1e-20
n_genes <- 20
setsource <- "all"

plot_module_projections(projectname, df, cellinfo, pval_cut, setsource)