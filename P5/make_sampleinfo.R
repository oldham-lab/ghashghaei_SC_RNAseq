setwd("/mnt/bdata/rebecca/collabs/ghashghaei/SC_expr/P5/aligned_reads")

samples <- list.dirs(recursive=F)

sample_ids <- sapply(strsplit(samples, "\\/"), "[", 2)
genotype <- sapply(strsplit(sample_ids, "_"), "[", 2)
protein <- sapply(strsplit(sample_ids, "_"), "[", 1)

sampleinfo <- data.frame(Sample_ID=sample_ids, 
                         Genotype=genotype, 
                         Protein=toupper(protein),
                         EGFR_Status="WT")
sampleinfo$EGFR_Status[with(sampleinfo, Protein == "GFP" & Genotype == "het")] <- "Null"
sampleinfo$EGFR_Status[with(sampleinfo, Genotype == "hom")] <- "Null"

write.csv(sampleinfo, file="../sampleinfo_SC_P5.csv", row.names=F)
