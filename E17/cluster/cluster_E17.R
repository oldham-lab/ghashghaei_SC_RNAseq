# E17

## Annotate cell types
### Ref: https://satijalab.org/seurat/articles/sctransform_v2_vignette.html

setwd("/mnt/bdata/rebecca/collabs/ghashghaei/SC_expr/E17/cluster")

source("/home/rebecca/code/misc/upper_first.R") ## Custom function converts first character of string to upper case, rest lowercase

## Merged QC counts for all samples:
expr <- fread("../counts/expression_QC.csv", data.table=F)
## Sample metadata:
sampleinfo <- read.csv("../sampleinfo_SC_E17.csv")

## Prep metadata:
sample_ids <- sapply(strsplit(colnames(expr[-c(1)]), "_"), function(x) paste(x[-length(x)], collapse="_"))
cellinfo <- data.frame(Cell_ID=colnames(expr)[-c(1)], Label=sample_ids)
sampleinfo <- sampleinfo[!duplicated(sampleinfo$Label),]
cellinfo <- merge(cellinfo, sampleinfo, by="Label")
cellinfo <- cellinfo %>% dplyr::select(-c(Sequencing_Date, Depth, Number_of_Cells_Loaded))
## Prep expression data:
genes <- expr[,c(1)]
expr <- as(expr[,-c(1)], "matrix")
rownames(expr) <- genes
colnames(expr) <- cellinfo$Cell_ID
rownames(cellinfo) <- cellinfo$Cell_ID

## Create Seurat object:
  
expr <- CreateSeuratObject(counts=expr, meta.data=cellinfo, row.names=genes)
expr[["percent.mt"]] <- PercentageFeatureSet(expr, pattern="^mt-")
expr$Group <- paste(expr$Genotype, expr$Protein)
## Subset to cells with < 10% mitochondrial reads:
expr <- subset(expr, subset=percent.mt < 10)

## Visualize reads and features per sample:

pdf(paste0("expr_nFeature_RNA.pdf"), width=12, height=7)
VlnPlot(expr, group.by = "Label", features = "nFeature_RNA", pt.size = 0.1)
dev.off()
pdf(paste0("expr_percent_mt.pdf"), width=12, height=7)
VlnPlot(expr, group.by = "Label", features = "percent.mt", pt.size = 0.1)
dev.off()

## How many cells per experimental group?
table(expr$Group)

## Median number reads per nucleus:
median(expr$nCount_RNA)

## Median number of genes expressed per nucleus:
median(expr$nFeature_RNA)

## No. genes and cells in total:
dim(expr)

## Adding cell cycle scores:
  
s.genes <- sapply(tolower(cc.genes$s.genes), upper_first)
g2m.genes <- sapply(tolower(cc.genes$g2m.genes), upper_first)
expr <- CellCycleScoring(expr, s.features=s.genes, g2m.features=g2m.genes, set.ident=T)

##Integrating samples:
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
res <- .7
expr_int <- FindClusters(expr_int, resolution=res)

## Visualize clusters:
pdf(paste0("SCT_UMAP_cluster_res_", res, ".pdf"), width=12, height=7)
DimPlot(expr_int, reduction="umap", pt.size=.3, label=T)
DimPlot(expr_int, reduction="umap", pt.size=.3, split.by="Group", group.by="Phase")
DimPlot(expr_int, reduction="umap", pt.size=.3, group.by="Library_Prep_Date")
DimPlot(expr_int, reduction="umap", pt.size=.3, group.by="Label")
FeaturePlot(expr_int, features="nFeature_RNA")
FeaturePlot(expr_int, features="percent.mt")
dev.off()

## Using sc-type tool for annotation:
## Ref: https://github.com/IanevskiAleksandr/sc-type/
  
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
## Using custom cell type marker list:
marker_genes <- readRDS("/mnt/bdata/rebecca/collabs/ghashghaei/SC_expr/marker_gene_list.RDS")
marker_genes <- lapply(marker_genes, HGNChelper::findExcelGeneSymbols)
gs_list <- list(gs_positive=marker_genes)
es.max = sctype_score(scRNAseqData = expr_int[["integrated"]]@scale.data, scaled = T, gs = gs_list$gs_positive, gs2=NULL, gene_names_to_uppercase=F)
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

## Visualize annotated clusters:
  
pdf("SCT_sc-type_auto_annotated_clsuters.pdf", width=10, height=7)
DimPlot(expr_int, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif')
dev.off()

## Get cluster markers:
  
# expr_int <- PrepSCTFindMarkers(expr_int, assay="SCT")
# markers <- FindAllMarkers(expr_int, assay="SCT")
# top_markers <- markers %>% 
#   dplyr::group_by(cluster) %>% 
#   dplyr::slice_head(n=15) %>%
#   as.data.frame()
# write.csv(top_markers, file="top_markers_rd1.csv", row.names=F)
# top_markers

## Calculate no./fraction of cells across cell types and experimental conditions:
  
metadata <- expr_int[[]]

n_cells_per_group <- metadata %>%
  dplyr::group_by(Group) %>%
  dplyr::summarise(No.Total_Nuclei=n()) %>%
  dplyr::arrange(No.Total_Nuclei)

n_cells_per_cycle <- metadata  %>%
  dplyr::group_by(Phase) %>%
  dplyr::summarise(No.Total_Nuclei=n()) %>%
  dplyr::arrange(No.Total_Nuclei)

n_cells_per_ct <- metadata  %>%
  dplyr::group_by(customclassif) %>%
  dplyr::summarise(No.Total_Nuclei=n()) %>%
  dplyr::arrange(No.Total_Nuclei)


n_cts_per_group <- metadata  %>%
  dplyr::group_by(Group) %>%
  dplyr::mutate(Total=n()) %>%
  dplyr::group_by(Group, Condition, customclassif) %>%
  dplyr::summarise(
    No.CT_Nuclei=n(),
    No.Total_Nuclei=unique(Total),
    Percent_CT_Nuclei=unique(n()/Total)*100,
  ) %>%
  dplyr::arrange(customclassif, Condition, Group)

n_cycles_per_group <- metadata  %>%
  dplyr::group_by(Group) %>%
  dplyr::mutate(Total=n()) %>%
  dplyr::group_by(Group, Condition, Phase) %>%
  dplyr::summarise(
    No.Phase_Nuclei=n(),
    No.Total_Nuclei=unique(Total),
    Percent_Phase_Nuclei=unique(n()/Total)*100,
  ) %>%
  dplyr::arrange(Phase, Group, Condition)

n_cycles_per_ct <- metadata  %>%
  dplyr::group_by(customclassif) %>%
  dplyr::mutate(Total=n()) %>%
  dplyr::group_by(customclassif, Phase) %>%
  dplyr::summarise(
    No.CT_Phase_Nuclei=n(),
    No.Total_CT_Nuclei=unique(Total),
    Percent_CT_Phase_Nuclei=unique(n()/Total)*100,
  ) %>%
  dplyr::arrange(customclassif, Phase)

## Visualize cell cycle by cell type:

pdf("phase_per_cell_type.pdf", width=14, height=7)
DimPlot(expr_int, reduction="umap", pt.size=.3, group.by="Phase", shuffle=T)
DimPlot(expr_int, reduction="umap", pt.size=.3, group.by="Phase", split.by="customclassif", shuffle=T, ncol=3)
dev.off()

# Save:
  
fwrite(n_cells_per_group, file="cells_per_experimental_group_E17.csv")
fwrite(n_cells_per_cycle, file="cells_per_phase_E17.csv")
fwrite(n_cells_per_ct, file="cells_per_cell_type_E17.csv")
fwrite(n_cts_per_group, file="cell_types_per_experimental_group_E17.csv")
fwrite(n_cycles_per_group, file="cycles_per_experimental_group_E17.csv")
fwrite(n_cycles_per_ct, file="cycles_per_cell_type_E17.csv")
saveRDS(expr_int, file="expr_SC_E17.RDS")
