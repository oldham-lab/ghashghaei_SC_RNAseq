library(dplyr)
library(Seurat)
library(data.table)
library(future.apply)

options(future.globals.maxSize=Inf)
plan(multicore, workers=10)

source("/home/rebecca/code/misc/normalize_fxns.R")

project_modules <- function(projectname, expr, fm_dir, expr_type, pval_cut, n_genes){
  
  ## Access expression data:
  mat <- expr@assays$RNA@counts
  if(expr_type == "normalized_counts"){
    mat <- log2_sparse(normalize_sparse(as(mat,"TsparseMatrix")))
  } 
  ## We will be stratifying data by cell type and experimental group:
  group_index <- tapply(1:ncol(expr), expr$Group, "[")
  clust_index <- tapply(1:ncol(expr), expr$Cluster, "[")
  networks <- list.files(path=fm_dir, pattern="signum", full.names=T)
  networks <- networks[unlist(lapply(networks, function(x) length(list.files(x))))>0]
  cat("\n")
  print("Finding maximally enriched modules...")
  network_list <- future_lapply(1:length(networks), FUN=function(i){
    print(i)
    ## Get maximally enriched module per gene set:
    enrich_paths <- list.files(path=networks[i], pattern="GSHyperG", full.names=T)
    enrich_list <- lapply(1:length(enrich_paths), function(j){
      enrich1 <- fread(enrich_paths[j], data.table=F)
      enrich <- suppressMessages({reshape2::melt(enrich1[,c(1, 8:ncol(enrich1))])})
      colnames(enrich)[2:3] <- c("Module", "Pval")
      enrich <- enrich[enrich$Pval < pval_cut,]
      if(nrow(enrich) > 0){
        ## Select maximially enriched module per set:
        enrich <- enrich %>% 
          dplyr::group_by(SetID) %>%
          dplyr::slice_min(Pval, with_ties=F) %>%
          dplyr::arrange(Pval)
        enrich <- merge(enrich, enrich1[,1:7], by="SetID")
        kme_paths <- list.files(path=networks[i], pattern="kME", full.names=T)
        return(data.frame(enrich, Network=networks[i], kME=kme_paths))
      } 
    })
    return(do.call(rbind, enrich_list))
  }, future.seed=T) 
  ## Select maximally enriched module per set across networks:
  signif_mods <- do.call(rbind, network_list)
  signif_mods <- signif_mods %>%
    dplyr::group_by(SetID) %>%
    dplyr::slice_min(Pval, with_ties=F)
  ## For each unique network...
  proj_list <- lapply(unique(signif_mods$Network), FUN=function(network){
    print(paste("Projecting modules from", network, "onto cells..."))
    signif_mods1 <- signif_mods[is.element(signif_mods$Network, network),]
    kme <- fread(signif_mods1$kME[1], data.table=F)
    ## Get list of top N genes by kME for maximally enriched modules:
    mod_genes <- lapply(unique(signif_mods1$Module), function(mod){
      index <- order(kme[,colnames(kme) == paste0("kME", mod)], decreasing=T)
      return(kme$Gene[index[1:n_genes]])
    })
    names(mod_genes) <- unique(signif_mods1$Module)
    ## For each experimental group:
    group_proj <- future_lapply(1:length(group_index), FUN=function(j){
      clust_index1 <- lapply(clust_index, function(x) intersect(x, group_index[[j]]))
      ## Stratify expression data by cell type:
      mat_clust <- lapply(clust_index1, function(index) as.matrix(mat[,index]))
      all_mean <- unlist(lapply(mat_clust, mean))
      ## For each module:
      mod_proj_list <- lapply(mod_genes, FUN=function(genes){
        ## Project top module genes onto cells:
        mod_mean <- unlist(lapply(mat_clust, function(mat){
          return(mean(mat[rownames(expr) %in% genes,]))
        }))
        std_err <- unlist(lapply(mat_clust, std_err_fxn))
        ## Normalize relative to all genes:
        proj_index <- mod_mean/all_mean
        std_err <- std_err/all_mean
        ## Normalize relative to all cell types:
        if(max(proj_index) > 0){
          proj_index <- proj_index/max(proj_index)
          std_err <- std_err/max(proj_index)
        }
        return(Projection_Index=proj_index, Std_Err=std_err)
      })
      temp <- data.frame(Module=names(mod_genes), do.call(rbind, mod_proj_list))
      mod_proj <- merge(signif_mods1, temp, by="Module")
      return(data.frame(Group=names(group_index)[j], mod_proj))
    })
    return(do.call(rbind, group_proj))
  })
  df <- do.call(rbind, proj_list)
  fwrite(df, file=paste0("data/", projectname, "_", expr_type, "_module_projections_pval_cut_", pval_cut, "_top_", n_genes, "_module_genes.csv"))
  
}

std_err_fxn <- function(x){sd(x)/sqrt(length(x))}
