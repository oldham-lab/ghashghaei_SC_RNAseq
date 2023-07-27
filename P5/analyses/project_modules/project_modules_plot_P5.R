setwd("/mnt/bdata/rebecca/collabs/ghashghaei/SC_expr/P5/analyses/project_modules")

source("/mnt/bdata/rebecca/collabs/ghashghaei/SC_expr/project_modules_plot_fxns.R")

expr <- readRDS("../../cluster/expr_SC_P5.RDS")
projectname <- "P5_SC_RNAseq"
expr_type <- "normalized_counts"
pval_cut <- 1e-20
setsource <- "MO"

cellinfo <- expr[[]]
cellinfo$Group <- gsub("het", "F+", cellinfo$Group)
cellinfo$Group <- gsub("hom", "FF", cellinfo$Group)
cellinfo$Group <- gsub("wt", "WT", cellinfo$Group)
cellinfo$Group <- gsub("TDT", "Td", cellinfo$Group)
cellinfo$Background <- "Homo"
cellinfo$Background[grep("F+", cellinfo$Group, fixed=T)] <- "Het"

n_genes <- 100
df <- fread("data/P5_SC_RNAseq_normalized_counts_module_projections_pval_cut_1e-20_top_100_module_genes.csv", data.table=F)
df$Group <- gsub("het", "F+", df$Group)
df$Group <- gsub("hom", "FF", df$Group)
df$Group <- gsub("wt", "WT", df$Group)
df$Group <- gsub("TDT", "Td", df$Group)

plot_module_projections(projectname, df, cellinfo, expr, expr_type, pval_cut, n_genes, setsource)

n_genes <- 100
df <- fread("data/P5_SC_RNAseq_normalized_counts_rel_group_module_projections_pval_cut_1e-20_top_100_module_genes.csv", data.table=F)
df$Group <- gsub("het", "F+", df$Group)
df$Group <- gsub("hom", "FF", df$Group)
df$Group <- gsub("wt", "WT", df$Group)
df$Group <- gsub("TDT", "Td", df$Group)

plot_rel_group_module_projections(projectname, df, cellinfo, expr, expr_type, pval_cut, n_genes, setsource)