setwd("/mnt/bdata/rebecca/collabs/ghashghaei/SC_expr/E17/analyses/FM/data/E17_SC_RNAseq_Modules")

source("/home/rebecca/code/GSEA/GSEA_functions.R")

n_threads <- 7
pvalcut1 <- NULL
exclude <- "none"
min_size <- 5
max_size <- 500
map_kme <- T
kme_unique_id <- "SYMBOL"
kme_unique_id_col <- 1 
tables_path <- "/home/rebecca/omicon/mapping_tables/2022-11-02"

##############################################################################################
############################################# MO #############################################
##############################################################################################

load("/home/rebecca/gene_sets/MO/MO_sets_mapped.RData")

set_source <- "MO"
legend <- MO_legend
gene_sets <- MO_sets_mapped

GSHGloop(set_source, gene_sets, legend, kmecut1="topmodposbc", topx=NULL, exclude, pvalcut1, min_size, max_size, map_sets=F, map_kme, kme_unique_id, kme_unique_id_col, platform=NULL, tables_path, n_threads)

cat("\nMyGSHGloop topmodposbc complete\n") 

##############################################################################################
############################################# Broad ##########################################
##############################################################################################

load("/home/rebecca/gene_sets/broad/broad_7.4_sets_mapped.RData")

set_source <- "BROAD"
legend <- broad_legend
gene_sets <- broad_sets_mapped

GSHGloop(set_source, gene_sets, legend, kmecut1="topmodposbc", topx=NULL, exclude, pvalcut1, min_size, max_size, map_sets=F, map_kme, kme_unique_id, kme_unique_id_col, platform=NULL, tables_path, n_threads)

cat("\nMyGSHGloop topmodposbc complete\n") 

