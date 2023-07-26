setwd("/mnt/bdata/rebecca/collabs/ghashghaei/SC_expr")

library(plyr)
library(Seurat)
library(data.table)

source("/home/rebecca/code/misc/upper_first.R")

load("/home/rebecca/gene_sets/MO/MO_sets_mapped_formatted.RData")

categ <- grep("CellType", MO_legend$Category)

############################################# Neuroepithelial cells #############################################

NEs <- c("Nes", "Notch1", "Sox10", "Cdh1", "Hmga2", "Ccnd1")

sets <- list(`Neuroepithelial cells`=NEs)
############################################# Radial glia #############################################

## Using 2015 Kriegstein supp. material:

oRGs <- "ACSBG1;PAQR8;FBN2;RAB3GAP2;CARHSP1;GABRB1;SEL1L3;STOM;GPR126;BRCA1;RHOJ;SEMA5A;LGALS3;HNMT;MT3;ITM2C;HS6ST2;ARAP2;DIO2;EDNRB;FUT8;HOPX;MOXD1;SLCO1C1;RGS16;TNC;EEPD1;TKTL1;LGALS3BP;GNG5;MLC1;NOG;MGST3;CDK6;IL6ST;DOCK1;LOC100507206;VEPH1;FHOD3;LRRTM3;LIFR;PMP22;ETV1;CDCA7L;FABP7;IDI1;MOB3B;FGFBP3;FAM107A;BMP7;TCF7L2;TTYH1;TFPI;APBB2;GULP1;SOAT1;CSDA;SLITRK2;NRG1;FKBP9;SAT1;GFAP;PTPRZ1;ETV5;TMEM132B;ABAT"
vRGs <- "TSPAN12;SNTG1;KIAA1217;LECT1;CRYAB;TBC1D1;TNFRSF19;NR4A1;NAMPT;PPARGC1A;NAPEPLD;FBXO32;MAFF;SAMD4A;TMEM47;CDON;STOX1;CTGF;PLCH1;DDIT3;NCKAP5;DNAJC1;TLE4;SHROOM3;LRIG3;PARD3B;PDGFD;PALLD;BMPR1B;DACH1;SALL1;SHISA2;LDHA"
RGs <- "VIM;PSAT1;ZFP36L1;PHGDH;ATP1A2;CLU;PTN;C8orf4;HES1;SPARC;SOX2;SLC1A3;GPX3;PON2;LTBP1;TFAP2C;SOX9;GLI3;SFRP1;DOK5;SALL3;CREB5;DBI;AXL;DDAH1;ID4;AKAP12;SLC35F1;COL11A1;PAX6;FOS;LIPG"
RGs <- sort(unique(unlist(strsplit(c(oRGs, vRGs, RGs), ";"))))

sets <- c(sets, `Radial glia`=list(RGs))

############################################# Intermediate progenitor cells (Tbr2/Eomes) #############################################

# index <- intersect(categ, grep("intermediate", MO_legend$SetName, ignore.case=T))
# MO_legend[index,]
# sort(MO_sets_mapped[[index]])
# IPCs <- MO_sets_mapped[[index]]

# sets <- list(`Intermediate progenitors`=c(MO_sets_mapped[[index]], "Tbr2"))

sets <- c(sets, list(`Intermediate progenitors`=c("Eomes", "Tbr2")))

############################################# Neurons #############################################
index <- intersect(categ, grep("BARRES_NEURONS", MO_legend$SetName, ignore.case=T))
MO_legend[index,]
sort(MO_sets_mapped[[index[1]]])

sets <- c(sets, `Mature neurons`=list(MO_sets_mapped[[index[1]]]))

############################################# Excitatory neurons #############################################

index <- intersect(categ, grep("excit", MO_legend$SetName, ignore.case=T))
MO_legend[index,]
sort(MO_sets_mapped[[index[1]]])

sets <- c(sets, `Excitatory neurons`=list(MO_sets_mapped[[index[1]]]))

############################################# Interneurons #############################################

index <- intersect(categ, grep("inhib", MO_legend$SetName, ignore.case=T))
MO_legend[index,]
sort(MO_sets_mapped[[index[1]]])

sets <- c(sets, `Inhibitory neurons`=list(MO_sets_mapped[[index[1]]]))

############################################# Astrocytes #############################################

index <- intersect(categ, grep("astro", MO_legend$SetName, ignore.case=T))
MO_legend[index,]
sort(MO_sets_mapped[[index[7]]])

sets <- c(sets, Astrocytes=list(MO_sets_mapped[[index[7]]]))

############################################# Microglia #############################################

index <- intersect(categ, grep("microglia", MO_legend$SetName, ignore.case=T))
MO_legend[index,]
sort(MO_sets_mapped[[index[11]]])

sets <- c(sets, Microglia=list(MO_sets_mapped[[index[11]]]))

############################################# Oligodendrocytes #############################################

index <- intersect(categ, grep("oligo", MO_legend$SetName, ignore.case=T))
MO_legend[index,]
sort(MO_sets_mapped[[index[9]]])
sets <- c(sets, Oligodendrocytes=list(MO_sets_mapped[[index[5]]]))

############################################# OPCs #############################################

index <- intersect(categ, grep("OPC", MO_legend$SetName, ignore.case=T))
MO_legend[index,]
sort(MO_sets_mapped[[index[1]]])

sets <- c(sets, OPCs=list(MO_sets_mapped[[index[1]]]))

############################################# Endothelial cells #############################################

index <- intersect(categ, grep("Endo", MO_legend$SetName, ignore.case=T))
MO_legend[index,]
sort(MO_sets_mapped[[index[2]]])

sets <- c(sets, Endothelial=list(MO_sets_mapped[[index[2]]]))

############################################# Dividing (mitosis) #############################################

# ## Using GO mitotic cell cycle genes:
# mitos <- fread("GO_term_summary_20230722_120753.csv", data.table=F)
# mitos <- mitos[mitos$Evidence != "",]

## Using Seurat list:
cc <- unlist(cc.genes, use.names=F)

sets <- c(sets, Dividing=list(cc))

############################################# Save #############################################

sets <- lapply(sets, function(x) unname(sapply(x, function(y){
  y <- upper_first(y, force_lower=T)
  gsub("rik", "Rik", y)
})))


saveRDS(sets, file="marker_gene_list.RDS")

df <- data.frame(lapply(sets, "length<-", max(lengths(sets))), check.names=F)

write.csv(df, file="marker_genes.csv", row.names=F)
