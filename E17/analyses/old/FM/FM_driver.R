setwd("/mnt/bdata/rebecca/collabs/ghashghaei/SC_expr/E17/analyses/FM")

library(dplyr)
library(Seurat)

source("/home/rebecca/code/FindModules/FindModules.R")

projectname <- "E17_SC_RNAseq"

expr <- readRDS("../../cluster/expr_SC_E17.RDS")
expr <- subset(expr, subset=Cluster!="EPEN")
expr$Cluster[expr$Cluster=="NA"] <- "SPhase"

cellinfo <-  expr[[]] %>%
  dplyr::mutate(idx=1:n()) %>%
  dplyr::sample_n(size=10e3) %>%
  dplyr::arrange(Cluster)

mat <- expr@assays$SCT@data
mat <- data.frame(Gene=rownames(expr), as.matrix(mat[,cellinfo$idx]))

## Use HVGs:

subset <- is.element(mat$Gene, expr@assays$integrated@var.features)
sum(subset)
# [1] 3000

samplegroups <- as.factor(cellinfo$Cluster)

FindModules(
  projectname=projectname,
  expr=mat,
  geneinfo=c(1),
  sampleindex=c(2:ncol(mat)),
  samplegroups=samplegroups,
  subset=subset,
  simMat=NULL,
  saveSimMat=FALSE,
  simType="Bicor",
  beta=1,
  overlapType="None",
  TOtype="signed",
  TOdenom="min",
  MIestimator="mi.mm",
  MIdisc="equalfreq",
  signumType="rel",
  iterate=TRUE,
  signumvec=c(0.99, 0.98, 0.97, 0.96, 0.95, .94, .93),
  minsizevec=c(8, 10, 12, 15, 20),
  signum=NULL,
  minSize=NULL,
  merge.by="ME",
  merge.param=0.85,
  export.merge.comp=FALSE,
  ZNCcut=2,
  calcSW=FALSE,
  loadTree=FALSE,
  writeKME=TRUE,
  calcBigModStat=FALSE,
  writeModSnap=TRUE
)

