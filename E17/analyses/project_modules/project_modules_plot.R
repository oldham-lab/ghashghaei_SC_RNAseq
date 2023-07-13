setwd("/mnt/bdata/rebecca/collabs/ghashghaei/SC_expr/E17/analyses/project_modules")

source("project_modules_plot_fxn.R")

df <- fread("data/E17_SC_RNAseq_counts_module_projections.csv", data.table=F)

expr <- readRDS("../../cluster/expr_SC_E17.RDS")

cellinfo <- expr[[]]
cellinfo <- cellinfo[!is.element(cellinfo$Cluster, "EPEN"),]
cellinfo$Cluster[cellinfo$Cluster=="NA"] <- "SPhase"
cellinfo$Background <- "Homo"
cellinfo$Background[grep("F+", cellinfo$Genotype, fixed=T)] <- "Het"

pvalcut <- 1e-20
setsource <- "all"

plot_module_projections(df, cellinfo, pvalcut, setsource)