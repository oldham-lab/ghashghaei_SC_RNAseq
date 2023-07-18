library(dplyr)
library(stringr)
library(ggplot2)
library(gridExtra)
library(data.table)
library(RColorBrewer)

plot_module_projections <- function(projectname, df, cellinfo, expr, expr_type, pval_cut, 
                                    n_genes, setsource=c("all", "MO", "BROAD")){

  if(setsource == "MO") df <- df[grep("^MO", df$SetID),]
  if(setsource == "BROAD") df <- df[!grepl("^MO", df$SetID),]
  temp <- cellinfo %>%
    dplyr::group_by(Group, Background, EGFR_Status) %>%
    dplyr::summarise()
  df <- merge(df, temp, by="Group")
  
  file_path <- paste0("figures/", projectname, "_", expr_type, "_module_projections_", setsource, "_sets_pval_cut_", pval_cut, "_top_", n_genes, "_module_genes.pdf")
  
  pdf(file=file_path, height=11, width=14)
  
  ## For each network:
  lapply(unique(df$Network), FUN=function(network){
    df1 <- df[df$Network %in% network,]
    kme <- fread(df1$kME[1], data.table=F)
    ## Subset kME table to genes present in SC data:
    kme <- kme[kme$Gene %in% rownames(expr),]
    ## For each module:
    lapply(unique(df1$Module), function(mod){
      df2 <- df1[df1$Module %in% mod,]
      plot_title <- paste(mod, "module\n", sapply(strsplit(network, "/"), function(x) x[length(x)]))
      ## Order genes by kME for working module:
      kme1 <- kme[,c(2, grep(paste0(mod, "$"), colnames(kme)))]
      kme1 <- kme1[order(kme1[,2], decreasing=T),]
      ## Make table of top 10 genes:
      seed_df <- data.frame(`Top genes by kME`=kme1$Gene[1:10], check.names=F)
      p_seed <- tableGrob(d=seed_df, theme=ttheme_minimal(base_size=12, padding=unit(c(2, 2), "mm")), rows=NULL)
      enrich <- df2 %>% 
        dplyr::group_by(SetName, Pval) %>%
        dplyr::summarise() %>%
        dplyr::ungroup() %>%
        dplyr::arrange(Pval) %>%
        dplyr::slice(1:10) 
      enrich$SetName <- factor(enrich$SetName, levels=enrich$SetName)
      ## Enrichment p-values barplot:
      p_enrich <- ggplot(enrich, aes(x=SetName, y=-log10(Pval))) +
        geom_bar(stat="identity", color="black", width=.6) +
        geom_hline(yintercept=-log10(.05), color="red") +
        theme_minimal() +
        theme(axis.text.y=element_text(size=8),
              axis.ticks.y=element_line(),
              axis.title.y=element_text(margin=margin(r=10)),
              axis.title.x=element_blank(),
              axis.text.x=element_text(angle=45, size=7, hjust=1, vjust=1),
              panel.border=element_rect(color="black", fill=NA),
              panel.grid.major.y=element_blank(),
              plot.margin=margin(1.5, .5, 1.5, .5, "cm")) +
        scale_x_discrete(labels=function(x) str_wrap(str_replace_all(x, "_" , " "), width=20)) +
        scale_y_continuous(expand=c(0, 0, 0, 0))
      ## Projection barplots:
      df2 <- df2 %>%
        dplyr::mutate(Cell_Type=gsub(".", " ", Cell_Type, fixed=T),
                      Err_Min=Projection_Index - 2*Std_Err,
                      Err_Max=Projection_Index + 2*Std_Err) %>%
        dplyr::group_by(Group, Cell_Type) %>%
        dplyr::slice(1)
      df2$Cell_Type_Label <- sapply(df2$Cell_Type, function(x) paste(strwrap(x, width=15), collapse="\n"))
      df2$EGFR_Status <- factor(df2$EGFR_Status, levels=c("WT", "Null"))
      df2$Background <- factor(df2$Background, levels=c("Homo", "Het"))
      df2$Group <- factor(df2$Group, levels=c("WT Td", "WT GFP", "FF Td", "FF GFP", "F+ Td", "F+ GFP"))
      ## Facet by experimental group:
      p1 <- ggplot(df2, aes(x=Cell_Type_Label, y=Projection_Index, fill=Cell_Type)) +
        geom_bar(stat="identity", color="black", show.legend=F) +
        geom_errorbar(aes(ymin=Err_Min, ymax=Err_Max), width=.2) + # position=position_dodge(.9)
        theme_minimal() +
        theme(plot.subtitle=element_text(margin=margin(b=7)),
              axis.title.x=element_blank(),
              axis.text.x=element_text(angle=45, hjust=1),
              axis.title.y=element_text(margin=margin(r=15)),
              panel.border=element_rect(color="black", fill=NA),
              plot.margin=margin(.5, .5, .5, .5, "cm"),
              panel.grid.major.x=element_blank(),
              panel.grid.minor.y=element_blank(),
              panel.grid.major.y=element_line(linewidth=.2, color="lightgrey"),
              strip.text.x=element_text(size=10)) +
        ylab("Projection index") +
        scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, by=.25)) +
        scale_fill_manual(values=brewer.pal(n_distinct(df2$Cell_Type), "Paired")) +
        facet_wrap(Group ~ ., scales="free_y", ncol=1)
      ## Facet by cell type:
      p2 <- ggplot(df2, aes(x=Group, y=Index, fill=EGFR_Status)) +
        geom_bar(stat="identity", color="black") +
        geom_errorbar(aes(ymin=Err_Min, ymax=Err_Max), width=.2) + # , position=position_dodge(.9)
        theme_minimal() +
        theme(plot.subtitle=element_text(margin=margin(b=7)),
              legend.title=element_text(face="bold"),
              legend.position="bottom",
              legend.direction="horizontal",
              axis.title.x=element_blank(),
              axis.text.x=element_text(size=9),
              axis.ticks.x=element_line(color="black"),
              axis.title.y=element_text(margin=margin(r=15)),
              panel.border=element_rect(color="black", fill=NA),
              plot.margin=margin(.5, .5, .5, .5, "cm"),
              panel.grid.major.x=element_blank(),
              panel.grid.minor.y=element_blank(),
              panel.grid.major.y=element_line(linewidth=.2, color="lightgrey"),
              strip.text.x=element_text(size=9, face="bold")) +
        scale_fill_manual(values=brewer.pal(3, "Set2")) +
        scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, by=.25)) +
        facet_wrap(. ~ Cell_Type, ncol=1, scales="free_y") +
        guides(fill=guide_legend(title="Egfr status"), 
               pattern=guide_legend(title="Background"))
      ## Arrange multiple plots on the same page:
      grid.arrange(p1, arrangeGrob(p_seed, p_enrich, ncol=1, heights=c(.5, 1)), ncol=2,
                   top=grid::textGrob(label=plot_title, gp=grid::gpar(fontsize=13)))
      grid.arrange(p2, arrangeGrob(p_seed, p_enrich, ncol=1, heights=c(.5, 1)), ncol=2,
                   top=grid::textGrob(label=plot_title, gp=grid::gpar(fontsize=13)))
    })
  })
  
  dev.off()
  
}