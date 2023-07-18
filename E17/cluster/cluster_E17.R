setwd("/mnt/bdata/rebecca/collabs/ghashghaei/SC_expr/E17/cluster")

library(dplyr)
library(Seurat)
library(data.table)

source("/home/rebecca/code/misc/upper_first.R")

expr <- fread("../counts/expression_QC.csv", data.table=F)
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

###############################################################################################################
############################################# Cluster & integrate #############################################
###############################################################################################################

## Ref: https://satijalab.org/seurat/articles/sctransform_v2_vignette.html

expr <- CreateSeuratObject(counts=expr, meta.data=cellinfo, row.names=genes)
expr[["percent.mt"]] <- PercentageFeatureSet(expr, pattern="^mt-")
expr$Group <- paste(expr$Genotype, expr$Protein)
table(expr$Group)
expr <- subset(expr, subset=percent.mt < 10)
pdf(paste0("expr_nFeature_RNA.pdf"), width=12, height=7)
VlnPlot(expr, group.by = "Sample_ID", features = "nFeature_RNA", pt.size = 0.1)
dev.off()
pdf(paste0("expr_percent_mt.pdf"), width=12, height=7)
VlnPlot(expr, group.by = "Sample_ID", features = "percent.mt", pt.size = 0.1)
dev.off()
expr <- subset(expr, subset=nFeature_RNA < 1e3)

median(expr$nCount_RNA)
# [1] 4659
median(expr$nFeature_RNA)
# [1] 2119
dim(expr)
# [1] 14452 46433

## Add cell cycle scores:
s.genes <- sapply(tolower(cc.genes$s.genes), upper_first)
g2m.genes <- sapply(tolower(cc.genes$g2m.genes), upper_first)
expr <- CellCycleScoring(expr, s.features=s.genes, g2m.features=g2m.genes, set.ident=T)

## First look at data w/o batch correction:

# set.seed(123)
# expr1 <- expr
# expr1 <- SCTransform(expr1, vars.to.regress=c('percent.mt'), vst.flavor="v2")
# expr1 <- CellCycleScoring(expr1, s.features=s.genes, g2m.features=g2m.genes, assay='SCT', set.ident=T)
# expr1 <- RunPCA(expr1, npcs=30)
# expr1 <- RunUMAP(expr1, reduction="pca", dims=1:30)
# 
# pdf("SCT_UMAP_no_integration.pdf", width=12, height=7)
# DimPlot(expr1, reduction="umap", pt.size=.3, label=T)
# DimPlot(expr1, reduction="umap", pt.size=.3, split.by="Group", group.by="Phase")
# DimPlot(expr1, reduction="umap", pt.size=.3, group.by="Library_Prep_Date")
# FeaturePlot(expr1, features="nFeature_RNA") 
# FeaturePlot(expr1, features="percent.mt")
# dev.off()
# 
# ## Clearly there are still batch effects w/o integration.
# rm(expr1)

## Integrate samples:

## Ref: https://satijalab.org/seurat/articles/integration_introduction.html

split <- SplitObject(expr, split.by="Label")
split <- lapply(X=split, FUN=function(x){
  set.seed(123)
  x <- SCTransform(x, vars.to.regress=c('percent.mt'), vst.flavor="v2") # 'nFeature_RNA', 'nCount_RNA'
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
DimPlot(expr_int, reduction="umap", pt.size=.3, label=T)
DimPlot(expr_int, reduction="umap", pt.size=.3, split.by="Group", group.by="Phase")
DimPlot(expr_int, reduction="umap", pt.size=.3, group.by="Library_Prep_Date")
DimPlot(expr_int, reduction="umap", pt.size=.3, group.by="Label")
FeaturePlot(expr_int, features="nFeature_RNA") 
FeaturePlot(expr_int, features="percent.mt")
dev.off()

## Using sc-type auto annotation for annotation:

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
es.max = sctype_score(scRNAseqData = expr_int[["integrated"]]@scale.data, scaled = T, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(expr_int@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(expr_int@meta.data[expr_int@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(expr_int@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
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

## Get cluster markers:

expr_int <- PrepSCTFindMarkers(expr_int, assay="SCT")
markers <- FindAllMarkers(expr_int, assay="SCT")
top_markers <- markers %>% 
  dplyr::group_by(cluster) %>% 
  dplyr::slice_head(n=15) %>%
  as.data.frame()
write.csv(top_markers, file="top_markers_rd1.csv", row.names=F)

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
  dplyr::group_by(Group, Condition, customclassif) %>%
  dplyr::summarise(
    No.CT_Nuclei=n(),
    No.Total_Nuclei=unique(Total),
    Percent_CT_Nuclei=unique(n()/Total)*100,
  ) %>%
  dplyr::arrange(customclassif, Condition, Group)

cycles_per_group <- expr_int[[]]  %>%
  dplyr::group_by(Group) %>%
  dplyr::mutate(Total=n()) %>%
  dplyr::group_by(Group, Condition, Phase) %>%
  dplyr::summarise(
    No.Phase_Nuclei=n(),
    No.Total_Nuclei=unique(Total),
    Percent_Phase_Nuclei=unique(n()/Total)*100,
  ) %>%
  dplyr::arrange(Phase, Group, Condition)

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

saveRDS(expr_int, file="expr_SC_E17.RDS")
