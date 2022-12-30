setwd("/mnt/bdata/rebecca/collabs/ghashghaei/SC_expr/E17/cluster")

library(Seurat)
library(data.table)
library(dplyr)

source("/home/rebecca/code/misc/upper_first.R")

expr <- fread("../counts/QC_counts/expression_QC.csv", data.table=F)
sampleinfo <- read.csv("../sampleinfo_SC_E17.csv")

## Prep metadata:

sample_ids <- sapply(strsplit(colnames(expr[-c(1)]), "_"), function(x) paste(x[-length(x)], collapse="_"))
cellinfo <- data.frame(Cell_ID=colnames(expr)[-c(1)], Label=sample_ids)
sampleinfo <- sampleinfo[!duplicated(sampleinfo$Label),]
cellinfo <- merge(cellinfo, sampleinfo, by="Label")
cellinfo <- cellinfo %>% dplyr::select(-c(Sequencing_Date, Depth, Number_of_Cells_Loaded))

genes <- expr[,c(1)]
expr <- as(expr[,-c(1)], "matrix")
rownames(expr) <- genes
colnames(expr) <- cellinfo$Cell_ID
rownames(cellinfo) <- cellinfo$Cell_ID

############################################# Cluster & integrate #############################################

## Ref: https://satijalab.org/seurat/articles/sctransform_v2_vignette.html

expr <- CreateSeuratObject(counts=expr, meta.data=cellinfo, row.names=genes)
expr[["percent.mt"]] <- PercentageFeatureSet(expr, pattern="^mt-")
expr$Group <- paste(expr$Genotype, expr$Protein)
expr <- subset(expr, subset=percent.mt<10)

## Get genes for cell cycle scoring:

s.genes <- sapply(tolower(cc.genes$s.genes), upper_first)
g2m.genes <- sapply(tolower(cc.genes$g2m.genes), upper_first)

## First look at data w/o batch correction:

expr1 <- expr

set.seed(123)
expr1 <- SCTransform(expr1, vars.to.regress=c('percent.mt'), vst.flavor="v2")
expr1 <- CellCycleScoring(expr1, s.features=s.genes, g2m.features=g2m.genes, assay='SCT', set.ident=T)
expr1 <- RunPCA(expr1, npcs=30)
expr1 <- RunUMAP(expr1, reduction="pca", dims=1:30)

pdf("SCT_UMAP_no_integration.pdf", width=12, height=7)
DimPlot(expr1, reduction="umap", pt.size=.3, label=T)
DimPlot(expr1, reduction="umap", pt.size=.3, split.by="Group", group.by="Phase")
DimPlot(expr1, reduction="umap", pt.size=.3, group.by="Library_Prep_Date")
FeaturePlot(expr1, features="nFeature_RNA") 
FeaturePlot(expr1, features="percent.mt")
dev.off()

## Clearly there are still batch effects w/o integration.

rm(expr1)

## Integrate samples:

## Ref: https://satijalab.org/seurat/articles/integration_introduction.html

split <- SplitObject(expr, split.by="Library_Prep_Date")
split <- lapply(X=split, FUN=function(x){
  set.seed(123)
  x <- SCTransform(x, vars.to.regress=c('percent.mt'), vst.flavor="v2") # 'nFeature_RNA', 'nCount_RNA'
  x <- CellCycleScoring(x, s.features=s.genes, g2m.features=g2m.genes, assay='SCT', set.ident=T)
})
features <- SelectIntegrationFeatures(object.list=split, nfeatures=3000)
split <- PrepSCTIntegration(object.list=split, anchor.features=features)
anchors <- FindIntegrationAnchors(object.list=split, anchor.features=features, normalization.method="SCT")
expr_int <- IntegrateData(anchorset=anchors, normalization.method="SCT")
expr_int <- RunPCA(expr_int, npcs=30)
expr_int <- RunUMAP(expr_int, reduction="pca", dims=1:30)
expr_int <- FindNeighbors(expr_int, reduction="pca", dims=1:30)
res <- 1
expr_int <- FindClusters(expr_int, resolution=res)

pdf(paste0("SCT_UMAP_cluster_res_", res, ".pdf"), width=12, height=7)
DimPlot(expr_int, reduction="umap", pt.size=.3, label=T)
DimPlot(expr_int, reduction="umap", pt.size=.3, split.by="Group", group.by="Phase")
DimPlot(expr_int, reduction="umap", pt.size=.3, group.by="Library_Prep_Date")
FeaturePlot(expr_int, features="nFeature_RNA") 
FeaturePlot(expr_int, features="percent.mt")
dev.off()

## Visualize canonical markers:

DefaultAssay(object=expr_int) <- "SCT"

pdf("SCT_canonical_markers.pdf", width=10, height=7)

## OG

FeaturePlot(
  object=expr_int, 
  features=c("Olig1", "Olig2", "Sox10"),
  label=T,
  repel=F,
  label.size=2,
  slot="data",
  pt.size=.1,
  order=T
)

## OPC

FeaturePlot(
  object=expr_int, 
  features=c("Pdgfra", "Cspg4"), 
  label=T,
  repel=F,
  label.size=2,
  pt.size=.1,
  order=T
)

## ASC

FeaturePlot(
  object=expr_int, 
  features=c("Apoe", "Slc1a3"), 
  label=T,
  repel=F,
  label.size=2,
  slot="data",
  pt.size=.1,
  order=T
)

## EXC

FeaturePlot(
  object=expr_int, 
  features=c("Neurod2", "Neurod6", "Slc17a6", "Slc17a7"),
  label=T,
  repel=F,
  label.size=2,
  pt.size=.1,
  order=T
)

## INH

FeaturePlot(
  object=expr_int, 
  features=c("Gad2", "Gad1", "Dlx6os1"), 
  label=T,
  repel=F,
  label.size=2,
  pt.size=.1,
  order=T
)

## EPEN

FeaturePlot(
  object=expr_int, 
  features=c("Foxj1"), 
  label=T,
  repel=F,
  label.size=2,
  pt.size=.1,
  order=T
)

## ENDO (note: we don't expend any of these cells in this data)

FeaturePlot(
  object=expr_int, 
  features=c("Pecam1", "Cd44"), 
  label=T,
  repel=F,
  label.size=2,
  slot="data",
  pt.size=.1,
  order=T
)

## MIC (note: we don't expend any of these cells in this data)

FeaturePlot(
  object=expr_int,
  features=c("Hpgds", "Mafb", "Ski"),
  label=T,
  repel=F,
  label.size=2,
  slot="data",
  pt.size=.1,
  order=T
)

dev.off()

## Consolidate clusters based on marker expression:

expr_int$Cluster <- "NA"
expr_int$Cluster[is.element(expr_int$seurat_clusters, c(12, 19))] <- "ASC"
expr_int$Cluster[is.element(expr_int$seurat_clusters, 22)] <- "OG"
expr_int$Cluster[is.element(expr_int$seurat_clusters, c(0:4, 6, 8, 11, 16:18, 21, 23))] <- "EXC"
expr_int$Cluster[is.element(expr_int$seurat_clusters, c(5, 7, 9:10, 13:15, 20, 24))] <- "INH"

## Add some more finely resolved annotations based on expression level:

expr_mat <- GetAssayData(object=expr_int, assay="SCT", slot="data")
expr_int$Cluster[expr_mat[rownames(expr_mat)=="Pdgfra"]>=1] <- "OPC"
expr_int$Cluster[expr_mat[rownames(expr_mat)=="Foxj1"]>=1] <- "EPEN"

Idents(expr_int) <- "Cluster"

## Plot annotated clusters:

pdf("SCT_annotated_clusters.pdf", width=10, height=7)
DimPlot(expr_int, reduction="umap", pt.size=.3, shuffle=T)
DimPlot(expr_int, reduction="umap", pt.size=.3, split.by="Group", shuffle=T)
dev.off()

## Make summary tables breaking down # nuclei per cell types and/or cell cycle per experimental group:

nuclei_per_group <- expr_int[[]] %>%
  dplyr::group_by(Group) %>%
  dplyr::summarise(No.Total_Nuclei=n()) %>%
  arrange(No.Total_Nuclei)

nuclei_per_cc <- expr_int[[]] %>%
  dplyr::group_by(Phase) %>%
  dplyr::summarise(No.Total_Nuclei=n()) %>%
  arrange(No.Total_Nuclei)

nuclei_per_ct <- expr_int[[]] %>%
  dplyr::group_by(Cluster) %>%
  dplyr::summarise(No.Total_Nuclei=n()) %>%
  arrange(No.Total_Nuclei)

cts_per_group <- expr_int[[]] %>%
  dplyr::group_by(Group) %>%
  dplyr::mutate(Total=n()) %>%
  dplyr::group_by(Group, Condition, Cluster) %>%
  dplyr::summarise(
    No.CT_Nuclei=n(),
    No.Total_Nuclei=unique(Total),
    Percent_CT_Nuclei=unique(n()/Total)*100,
  ) %>%
  arrange(Cluster, Condition, Group)

cc_per_group <- expr_int[[]] %>%
  dplyr::group_by(Group) %>%
  dplyr::mutate(Total=n()) %>%
  dplyr::group_by(Group, Condition, Phase) %>%
  dplyr::summarise(
    No.Phase_Nuclei=n(),
    No.Total_Nuclei=unique(Total),
    Percent_Phase_Nuclei=unique(n()/Total)*100,
  ) %>%
  arrange(Phase, Group, Condition)

cc_per_ct <- expr_int[[]] %>%
  dplyr::group_by(Cluster) %>%
  dplyr::mutate(Total=n()) %>%
  dplyr::group_by(Cluster, Phase) %>%
  dplyr::summarise(
    No.CT_Phase_Nuclei=n(),
    No.Total_CT_Nuclei=unique(Total),
    Percent_CT_Phase_Nuclei=unique(n()/Total)*100,
  ) %>%
  arrange(Cluster, Phase)

pdf("SCT_CC_phase_per_CT.pdf", width=14, height=7)
DimPlot(expr_int, reduction="umap", pt.size=.3, group.by="Phase", shuffle=T)
DimPlot(expr_int, reduction="umap", pt.size=.3, group.by="Phase", split.by="Cluster", shuffle=T)
dev.off()

fwrite(nuclei_per_group, file="SCT_nuclei_per_experimental_group.csv")
fwrite(nuclei_per_cc, file="SCT_nuclei_per_CC_phase.csv")
fwrite(nuclei_per_ct, file="SCT_nuclei_per_CT.csv")
fwrite(cts_per_group, file="SCT_CTs_per_experimental_group.csv")
fwrite(cc_per_group, file="SCT_CC_phase_per_experimental_group.csv")
fwrite(cc_per_ct, file="SCT_CC_phase_per_CT.csv")

## Find cluster markers:

expr_int <- PrepSCTFindMarkers(expr_int)
markers <- FindAllMarkers(expr_int, assay="SCT")

top_markers <- markers %>% 
  dplyr::group_by(cluster) %>% 
  dplyr::slice_head(n=10) %>%
  as.data.frame()

fwrite(top_markers, file="SCT_top_markers_per_CT.csv")

saveRDS(expr_int, file="expr_SC_E17.RDS")

############################################# DE genes within cell types #############################################

## Ref: https://bookdown.org/ytliu13207/SingleCellMultiOmicsDataAnalysis/pseudobulk-expression.html