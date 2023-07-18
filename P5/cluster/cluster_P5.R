setwd("/mnt/bdata/rebecca/collabs/ghashghaei/SC_expr/P5/cluster")

library(dplyr)
library(Seurat)
library(data.table)

source("/home/rebecca/code/misc/upper_first.R")

expr <- readRDS("../counts/expression_QC.RDS")
sampleinfo <- read.csv("../sampleinfo_SC_P5.csv")

## Prep metadata:

sample_ids <- sapply(strsplit(colnames(expr), "_"), function(x) paste(x[-length(x)], collapse="_"))
cellinfo <- data.frame(Cell_ID=colnames(expr), Sample_ID=sample_ids)
cellinfo <- merge(cellinfo, sampleinfo, by="Sample_ID", sort=F)
rownames(cellinfo) <- cellinfo$Cell_ID

table(cellinfo$Sample_ID)
# gfp_het_1 gfp_het_2 gfp_hom_1 gfp_hom_2  gfp_wt_1  gfp_wt_2 tdt_het_1 tdt_het_2 
# 4847      4603      5834      7759      6030      2306      8924      3751 
# tdt_hom_1 tdt_hom_2  tdt_wt_1  tdt_wt_2 
# 4216      6740      5484      2970 

###############################################################################################################
############################################# Cluster & integrate #############################################
###############################################################################################################

## Ref: https://satijalab.org/seurat/articles/sctransform_v2_vignette.html

expr <- CreateSeuratObject(counts=expr, meta.data=cellinfo)
expr[["percent.mt"]] <- PercentageFeatureSet(expr, pattern="^mt-")
expr$Group <- paste(expr$Genotype, expr$Protein)
expr <- subset(expr, subset=percent.mt < 10)
pdf(paste0("expr_nFeature_RNA.pdf"), width=12, height=7)
VlnPlot(expr, group.by = "Sample_ID", features = "nFeature_RNA", pt.size = 0.1)
dev.off()
pdf(paste0("expr_nCount_RNA.pdf"), width=12, height=7)
VlnPlot(expr, group.by = "Sample_ID", features = "nCount_RNA", pt.size = 0.1)
dev.off()
pdf(paste0("expr_percent_mt.pdf"), width=12, height=7)
VlnPlot(expr, group.by = "Sample_ID", features = "percent.mt", pt.size = 0.1)
dev.off()

median(expr$nCount_RNA)
# [1] 15191
median(expr$nFeature_RNA)
# [1] 4457
dim(expr)
# [1] 21910 61306

s.genes <- sapply(tolower(cc.genes$s.genes), upper_first)
g2m.genes <- sapply(tolower(cc.genes$g2m.genes), upper_first)
expr <- CellCycleScoring(expr, s.features=s.genes, g2m.features=g2m.genes, set.ident=T)

## Integrate samples:

## Ref: https://satijalab.org/seurat/articles/integration_introduction.html

split <- SplitObject(expr, split.by="Sample_ID")
split <- lapply(X=split, FUN=function(x){
  set.seed(123)
  x <- SCTransform(x, vars.to.regress=c('percent.mt'), vst.flavor="v2")
})
features <- SelectIntegrationFeatures(object.list=split, nfeatures=2000)
split <- PrepSCTIntegration(object.list=split, anchor.features=features)
anchors <- FindIntegrationAnchors(object.list=split, anchor.features=features, normalization.method="SCT")
expr_int <- IntegrateData(anchorset=anchors, normalization.method="SCT")
expr_int <- RunPCA(expr_int, npcs=30)
expr_int <- RunUMAP(expr_int, reduction="pca", dims=1:30)
expr_int <- FindNeighbors(expr_int, reduction="pca", dims=1:30)
res <- 1
expr_int <- FindClusters(expr_int, resolution=res)

pdf(paste0("SCT_UMAP_cluster_res_", res, ".pdf"), width=12, height=7)
DimPlot(expr_int, reduction="umap", pt.size=.3, label=T, group.by="seurat_clusters")
DimPlot(expr_int, reduction="umap", pt.size=.1, group.by="Group", shuffle=T)
DimPlot(expr_int, reduction="umap", pt.size=.1, group.by="Sample_ID", shuffle=T)
DimPlot(expr_int, reduction="umap", pt.size=.3, split.by="Group", group.by="Phase")
FeaturePlot(expr_int, features="nFeature_RNA")
FeaturePlot(expr_int, features="nCount_RNA")
FeaturePlot(expr_int, features="percent.mt", split.by="Sample_ID")
dev.off()

## Using sc-type auto annotation for a first pass at annotation:

## Ref: https://github.com/IanevskiAleksandr/sc-type/
library(HGNChelper)
# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
# DB file
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Brain"
# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)
es.max = sctype_score(scRNAseqData = expr_int[["integrated"]]@scale.data, 
                      scaled = T, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(expr_int@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(expr_int@meta.data[expr_int@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(expr_int@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] <- "Unknown"
expr_int@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,];
  expr_int@meta.data$customclassif[expr_int@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}
Idents(expr_int) <- "customclassif"
  
pdf("SCT_sc-type_auto_annotated_clsuters.pdf", width=10, height=7)
DimPlot(expr_int, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif')
DimPlot(expr_int, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'seurat_clusters')
dev.off()

## Visualize canonical markers:

DefaultAssay(object=expr_int) <- "RNA"

pdf("SCT_canonical_markers.pdf", width=10, height=7)
## Oligos
FeaturePlot(object=expr_int, features=c("Olig1", "Olig2", "Sox10", "Mog"), label=T, repel=F, label.size=2, slot="data", pt.size=.1, order=T)
## OPCs
FeaturePlot(object=expr_int, features=c("Pdgfra", "Cspg4", "Pcdh15", "Lhfpl3"), label=T, repel=F, label.size=2, slot="data", pt.size=.1, order=T)
## Asctrocytes
FeaturePlot(object=expr_int, features=c("Apoe", "Gfap", "Slc1a3"), label=T, repel=F, label.size=2, slot="data", pt.size=.1, order=T)
## Glutamatergic
FeaturePlot(object=expr_int, features=c("Neurod2", "Neurod6", "Slc17a6", "Slc17a7"), label=T, repel=F, label.size=2, slot="data", pt.size=.1, order=T)
## GABAergic
FeaturePlot(object=expr_int, features=c("Gad2", "Gad1", "Dlx6os1"), label=T, repel=F, label.size=2, slot="data", pt.size=.1, order=T)
## Immature neurons
FeaturePlot(object=expr_int, features=c("Ncam1", "Neurod1"), label=T, repel=F, label.size=2, slot="data", pt.size=.1, order=T)
## Mature neurons
FeaturePlot(object=expr_int, features=c("Dlg4", "Gap43", "Map2", "Rbfox3"), label=T, repel=F, label.size=2, slot="data", pt.size=.1, order=T)
## Ependymal cells
FeaturePlot(object=expr_int, features=c("Foxj1"), label=T, repel=F, label.size=2, slot="data", pt.size=.1, order=T)
## Fibroblasts
FeaturePlot(object=expr_int, features=c("Col1a1", "Col1a2", "Col5a1", "Fbln1"), label=T, repel=F, label.size=2, slot="data", pt.size=.1, order=T)
## Radial glia
FeaturePlot(object=expr_int, features=c("Vim", "Hes1", "Hes5", "Cdh2"), label=T, repel=F, label.size=2, slot="data", pt.size=.1, order=T)
## Neural progenitor cells
FeaturePlot(object=expr_int, features=c("Ascl1", "Sox1", "Smarca4"), label=T, repel=F, label.size=2, slot="data", pt.size=.1, order=T)
## Neural stem cells
FeaturePlot(object=expr_int, features=c("Sox9", "Prom1"), label=T, repel=F, label.size=2, slot="data", pt.size=.1, order=T)
## Cycling
FeaturePlot(object=expr_int, features=c("Hmgb2", "Top2a", "Cdk6"), label=T, repel=F, label.size=2, slot="data", pt.size=.1, order=T)
dev.off()

# ## Annotate clusters based on expression:
# 
# expr_int$My_Anno <- as.character(expr_int$customclassif)
# expr_int$My_Anno[is.element(expr_int$seurat_clusters, c(36))] <- "Radial glia"
# expr_int$My_Anno[is.element(expr_int$seurat_clusters, c(16, 28, 33))] <- "Neural progenitor cells"
# expr_int$My_Anno[is.element(expr_int$seurat_clusters, c(6, 9, 11, 13, 25, 36))] <- "Astrocyte"
# expr_int$My_Anno[is.element(expr_int$seurat_clusters, c(4, 14))] <- "Glutamatergic neuron"
# expr_int$My_Anno[is.element(expr_int$seurat_clusters, c(0:3, 7, 38))] <- "GABAergic neuron"
# expr_int$My_Anno[is.element(expr_int$seurat_clusters, c(20))] <- "Mature neuron"
# expr_int$My_Anno[is.element(expr_int$seurat_clusters, c(23))] <- "Immature neuron"
# expr_int$My_Anno[is.element(expr_int$seurat_clusters, c(30))] <- "Oligodendrocyte"
# expr_int$My_Anno[is.element(expr_int$seurat_clusters, c(8, 19, 27, 32, 27))] <- "OPC"
# Idents(expr_int) <- "My_Anno"

## Get cluster markers

expr_int <- PrepSCTFindMarkers(expr_int, assay="SCT")
markers <- FindAllMarkers(expr_int, assay="SCT")
top_markers <- markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_head(n=15) %>%
  as.data.frame()
write.csv(top_markers, file="top_markers_rd1.csv", row.names=F)

# ## Add additional annotations:
# expr_int$My_Anno[is.element(expr_int$seurat_clusters, c(15))] <- "Ependymal cell" # Tmem212, Rarres2, Cfap299
# expr_int$My_Anno[is.element(expr_int$seurat_clusters, c(10, 26, 34, 21))] <- "Neural progenitor cells" # Hist1h2ap, Top2a, Hmgb2, Pclaf, Ccnd2, Hmgn2, Neurod1, Ascl1
# expr_int$My_Anno[is.element(expr_int$seurat_clusters, c(12, 17))] <- "GABAergic neuron"
# expr_int$My_Anno[is.element(expr_int$seurat_clusters, c(22, 24, 41))] <- "Mature neuron"
# expr_int$My_Anno[is.element(expr_int$seurat_clusters, c(27))] <- "OPC"
# expr_int$My_Anno[is.element(expr_int$seurat_clusters, c(31, 5, 29))] <- "Astrocyte"
# expr_int$My_Anno[is.element(expr_int$seurat_clusters, c(18, 39, 40))] <- "Radial glia" # Rsph1, Sparc, Tmem212, Fabp7, Gli2, Gli3, Notch2
# expr_int$My_Anno[is.element(expr_int$seurat_clusters, c(37))] <- "Pericytes and Vascular leptomeningeal" # Vtn, Igf2, Apod, Dcn, Col1a2
# 
# ## Plot annotated clusters:
# 
# pdf("SCT_annotated_clusters.pdf", width=10, height=7)
# DimPlot(expr_int, reduction="umap", pt.size=.3, shuffle=T)
# DimPlot(expr_int, reduction="umap", pt.size=.3, split.by="Group", shuffle=T)
# dev.off()

## Make summary tables breaking down # nuclei per cell types and/or cell cycle and experimental group:

nuclei_group <- expr_int[[]] %>%
  dplyr::group_by(Group) %>%
  dplyr::summarise(No.Total_Nuclei=n()) %>%
  dplyr::arrange(No.Total_Nuclei)

nuclei_cycle <- expr_int[[]]  %>%
  dplyr::group_by(Phase) %>%
  dplyr::summarise(No.Total_Nuclei=n()) %>%
  dplyr::arrange(No.Total_Nuclei)

nuclei_cts <- expr_int[[]]  %>%
  dplyr::group_by(customclassif) %>%
  dplyr::summarise(No.Total_Nuclei=n()) %>%
  dplyr::arrange(No.Total_Nuclei)

cts_per_group <- expr_int[[]]  %>%
  dplyr::group_by(Group) %>%
  dplyr::mutate(Total=n()) %>%
  dplyr::group_by(Group, EGFR_Status, customclassif) %>%
  dplyr::summarise(
    No.CT_Nuclei=n(),
    No.Total_Nuclei=unique(Total),
    Percent_CT_Nuclei=unique(n()/Total)*100,
  ) %>%
  dplyr::arrange(customclassif, EGFR_Status, Group)

cycles_per_group <- expr_int[[]]  %>%
  dplyr::group_by(Group) %>%
  dplyr::mutate(Total=n()) %>%
  dplyr::group_by(Group, EGFR_Status, Phase) %>%
  dplyr::summarise(
    No.Phase_Nuclei=n(),
    No.Total_Nuclei=unique(Total),
    Percent_Phase_Nuclei=unique(n()/Total)*100,
  ) %>%
  dplyr::arrange(Phase, Group, EGFR_Status)

cycles_per_ct <- expr_int[[]]  %>%
  dplyr::group_by(customclassif) %>%
  dplyr::mutate(Total=n()) %>%
  dplyr::group_by(customclassif, Phase) %>%
  dplyr::summarise(
    No.CT_Phase_Nuclei=n(),
    No.Total_CT_Nuclei=unique(Total),
    Percent_CT_Phase_Nuclei=unique(n()/Total)*100,
  ) %>%
  dplyr::arrange(customclassif, Phase)

## Plot:

pdf("phase_per_cell_type.pdf", width=14, height=7)
DimPlot(expr_int, reduction="umap", pt.size=.3, group.by="Phase", shuffle=T)
DimPlot(expr_int, reduction="umap", pt.size=.3, group.by="Phase", split.by="customclassif", shuffle=T, ncol=3)
dev.off()

## Save:

fwrite(nuclei_group, file="nuclei_per_experimental_group.csv")
fwrite(nuclei_cycle, file="nuclei_per_phase.csv")
fwrite(nuclei_cts, file="nuclei_per_cell_type.csv")
fwrite(cts_per_group, file="cell_types_per_experimental_group.csv")
fwrite(cycles_per_group, file="phase_per_experimental_group.csv")
fwrite(cycles_per_ct, file="phase_per_cell_type.csv")

## Save annotated data:

saveRDS(expr_int, file="expr_SC_P5.RDS")

