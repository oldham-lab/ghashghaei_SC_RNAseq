setwd("/mnt/bdata/rebecca/collabs/ghashghaei/SC_expr/P5/analyses/cell_proportions")

library(dplyr)
library(scales)
library(Seurat)
library(ggplot2)
library(ggsignif)
library(ggpattern)
library(RColorBrewer)

expr <- readRDS("../../cluster/expr_SC_P5.RDS")

expr$Genotype <- gsub("het", "F+", expr$Genotype)
expr$Genotype <- gsub("hom", "FF", expr$Genotype)
expr$Genotype <- gsub("wt", "WT", expr$Genotype)
expr$Background <- "Homo"
expr$Background[grep("F+", expr$Genotype, fixed=T)] <- "Het"
expr$Protein <- gsub("TDT", "Td", expr$Protein)

df <- expr[[]] %>%
  dplyr::group_by(customclassif) %>%
  dplyr::mutate(Group=paste(Genotype, Protein)) %>%
  dplyr::group_by(customclassif, Group, Background, EGFR_Status) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::group_by(Group) %>%
  dplyr::mutate(Frac=n/sum(n))

df$Group <- factor(df$Group, levels=c("WT Td", "WT GFP", "FF Td", "FF GFP", "F+ Td", "F+ GFP"))
df$EGFR_Status <- factor(df$EGFR_Status, levels=c("WT", "Null"))
df$customclassif <- sapply(df$customclassif, function(x) paste(strwrap(x, width=15), collapse="\n"))

chi.test <- function(a, b) {
  return(chisq.test(cbind(a, b)))
}

sig_fxn <- function(x){
  if(x<.001) 
    "***" 
  else if(x<.01) 
    "**"
  else if(x<.05) 
    "*"
  else NA
}

compare_list <- list(c("F+ GFP", "F+ Td"),
                     c("WT GFP", "WT Td"),
                     c("FF GFP", "FF Td"))

pdf("figures/P5_cells_by_group_and_cell_type_barplot.pdf", width=14, height=7)

print(
  ggplot(df, aes(x=Group, y=Frac, fill=EGFR_Status)) +
    geom_col_pattern(aes(pattern=Background),
                     position=position_dodge(0),
                     color="black",
                     pattern_color="white",
                     pattern_density=.04,
                     pattern_spacing=.04,
                     pattern_key_scale_factor=.3) +
    theme_minimal() +
    theme(plot.subtitle=element_text(margin=margin(b=7)),
          legend.title=element_text(face="bold"),
          legend.position="bottom",
          legend.direction="horizontal",
          axis.title.x=element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x=element_line(color="black"),
          axis.title.y=element_text(margin=margin(r=15)),
          plot.margin=margin(.5, .5, .5, .5, "cm"),
          panel.grid.major.x=element_blank(),
          panel.grid.minor.y=element_line(linewidth=.2, color="lightgrey"),
          panel.grid.major.y=element_line(linewidth=.2, color="lightgrey")) +
    labs(title="P5 single-cell RNA-seq", 
         subtitle=paste(comma(nrow(expr)), "total cells")) +
    ylab("Fraction of cells\n(per experimental group)") +
    scale_fill_manual(values=brewer.pal(3, "Set2")) +
    guides(fill=guide_legend(title="Egfr status"),
           pattern=guide_legend(title="Background")) +
    facet_wrap(. ~ customclassif, nrow=1, strip.position="bottom") +
    theme(strip.text.x=element_text(size=9),
          strip.placement="outside",
          panel.spacing.x=unit(0, "lines")) 
)

# geom_signif(data=subset(df, Cluster=="ASC"),
#             comparisons=compare_list, test="chi.test",
#             map_signif_level=sig_fxn,
#             y_position=max(subset(df, Cluster=="ASC")$n)+20) +
# geom_signif(data=subset(df, Cluster=="OG"),
#             comparisons=compare_list, test="chi.test",
#             map_signif_level=sig_fxn,
#             y_position=max(subset(df, Cluster=="OG")$n)+20) +
# geom_signif(data=subset(df, Cluster=="EXC"),
#             comparisons=compare_list, test="chi.test",
#             map_signif_level=sig_fxn,
#             y_position=4500) +
# geom_signif(data=subset(df, Cluster=="INH"),
#             comparisons=compare_list, test="chi.test",
#             map_signif_level=sig_fxn,
#             y_position=max(subset(df, Cluster=="INH")$n)+20) +
# geom_signif(data=subset(df, Cluster=="OPC"),
#             comparisons=compare_list, test="chi.test",
#             map_signif_level=sig_fxn,
#             y_position=max(subset(df, Cluster=="OPC")$n)+20) +
# geom_signif(data=subset(df, Cluster=="SPhase"),
#             comparisons=compare_list, test="chi.test",
#             map_signif_level=sig_fxn,
#             y_position=max(subset(df, Cluster=="SPhase")$n)+20) 

dev.off()

