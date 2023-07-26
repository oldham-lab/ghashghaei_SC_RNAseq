setwd("/mnt/bdata/rebecca/collabs/ghashghaei/SC_expr/E17/analyses/project_modules")

source("/mnt/bdata/rebecca/collabs/ghashghaei/SC_expr/project_modules_plot_fxns.R")

expr <- readRDS("../../cluster/expr_SC_E17.RDS")
projectname <- "E17_SC_RNAseq"
expr_type <- "normalized_counts"
pval_cut <- 1e-20
setsource <- "MO"

cellinfo <- expr[[]]
cellinfo$Background <- "Homo"
cellinfo$Background[grep("F+", cellinfo$Genotype, fixed=T)] <- "Het"
cellinfo$EGFR_Status <- cellinfo$Condition

n_genes <- 100
df <- fread("data/E17_SC_RNAseq_normalized_counts_module_projections_pval_cut_1e-20_top_100_module_genes.csv", data.table=F)

plot_module_projections(projectname, df, cellinfo, expr, expr_type, pval_cut, n_genes, setsource)

n_genes <- 100
df <- fread("data/E17_SC_RNAseq_normalized_counts_rel_group_module_projections_pval_cut_1e-20_top_100_module_genes.csv", data.table=F)

plot_rel_group_module_projections(projectname, df, cellinfo, expr, expr_type, pval_cut, n_genes, setsource)
