library(dplyr)
library(stringr)
library(ggplot2)
library(gridExtra)
library(data.table)
library(RColorBrewer)

plot_module_projections <- function(df, cellinfo, pvalcut=.05, setsource="all"){

  if(setsource=="MO"){
    
    df <- df[grep("MO", df$SetID),]
    
  } else if(setsource=="BROAD"){
    
    df <- df[!grepl("MO", df$SetID),]
    
  }
  
  temp <- cellinfo %>%
    dplyr::group_by(Group, Background, Condition) %>%
    dplyr::summarise()
  
  df <- merge(df, temp, by="Group")

  pdf(file=paste0("figures/projection_plots_", setsource, "_sets_pval<", pvalcut, ".pdf"), height=10, width=8.5)
  
  ## For each unique network...
  
  lapply(unique(df$Network), FUN=function(network){

    df1 <- df[is.element(df$Network, network),]
    kme <- fread(df1$kME[1], data.table=F)
    mods <- unique(df1$Module)
    
    ## For each module...
    
    lapply(1:length(mods), function(i){
      
      ## Plot module projections:

      df2 <- df1[is.element(df1$Module, mods[i]),]
      
      if(!is.na(df2$ASC[1])&min(df2$Pval)<pvalcut){

        plot_title <- paste(df2$Module[1], "module\n", 
                            sapply(strsplit(network, "/"), function(x) x[length(x)]))
        
        kme1 <- kme[,c(2, grep(paste0(df2$Module[1], "$"), colnames(kme)))]
        kme1 <- kme1[order(-kme1[,2]),]
        
        ## Make seed genes table:
        
        p_seed <- tableGrob(d=data.frame(`Seed Genes`=kme1$Gene[1:10], check.names=F), 
                            theme=ttheme_minimal(base_size=12, padding=unit(c(2, 2), "mm")),
                            rows=NULL)
        
        ## Make enrichments barplot:
        
        enrich <- df2 %>% 
          dplyr::group_by(SetName, Pval, SetSize) %>%
          dplyr::summarise() %>%
          dplyr::ungroup() %>%
          dplyr::arrange(Pval) %>%
          dplyr::slice(1:10)
        
        enrich$SetName <- factor(enrich$SetName, levels=enrich$SetName)
        
        p_enrich <- ggplot(enrich, aes(x=SetName, y=-log10(Pval))) +
          geom_bar(stat="identity", color="black", width=.6) +
          geom_hline(yintercept=-log10(.05), color="red") +
          theme_minimal() +
          theme(axis.text.y=element_text(size=8),
                axis.ticks.y=element_line(),
                axis.title.y=element_text(margin=margin(t=10)),
                axis.title.x=element_blank(),
                axis.text.x=element_text(angle=45, size=7, hjust=1, vjust=1),
                panel.border=element_rect(color="black", fill=NA),
                panel.grid.major.y=element_blank(),
                plot.margin=margin(1.5, .5, 1.5, .5, "cm")) +
          scale_x_discrete(labels=function(x) str_wrap(str_replace_all(x, "_" , " "), width=30)) +
          scale_y_continuous(expand=c(.01, .01, .01, .01))
        
        ## Make projection barplots:
        
        dat <- reshape2::melt(df2[,c(1, 13:20)]) 
        colnames(dat)[4:5] <- c("Cluster", "Index")
        dat <- dat %>%
          # dplyr::filter(Cluster!="SPhase") %>%
          dplyr::group_by(Group, Background, Condition, Cluster, Index) %>%
          dplyr::summarise()
      
        dat$Condition <- factor(dat$Condition, levels=c("WT", "Null"))
        dat$Background <- factor(dat$Background, levels=c("Homo", "Het"))
        dat$Group <- factor(dat$Group, levels=c("WT Td", "WT GFP", 
                                                "FF Td", "FF GFP", 
                                                "F+ Td", "F+ GFP"))
        
        ## Facet by experimental group:
        
        p1 <- ggplot(dat, aes(x=Cluster, y=Index, fill=Cluster)) +
          geom_col(color="black", show.legend=F) +
          theme_minimal() +
          theme(plot.subtitle=element_text(margin=margin(b=7)),
                legend.title=element_text(face="bold"),
                legend.position="bottom",
                legend.direction="horizontal",
                axis.title.x=element_blank(),
                axis.text.x=element_text(size=9),
                axis.title.y=element_text(margin=margin(r=15)),
                panel.border=element_rect(color="black", fill=NA),
                plot.margin=margin(.5, .5, .5, .5, "cm"),
                panel.grid.major.x=element_blank(),
                panel.grid.minor.y=element_blank(),
                panel.grid.major.y=element_line(linewidth=.2, color="lightgrey")) +
          scale_fill_manual(values=brewer.pal(6, "Paired")) +
          scale_y_continuous(breaks=seq(0, 1, by=.25), expand=c(.01, .1, .1, .1)) +
          facet_wrap(.~Group, ncol=1, scales="free_x") +
          theme(strip.text.x=element_text(size=10)) 
        
        ## Facet by cell type:
        
        p2 <- ggplot(dat, aes(x=Group, y=Index, fill=Condition)) +
          geom_bar(stat="identity", color="black") +
          theme_minimal() +
          theme(plot.subtitle=element_text(margin=margin(b=7)),
                legend.title=element_text(face="bold"),
                legend.position="bottom",
                legend.direction="horizontal",
                axis.title.x=element_blank(),
                axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=9),
                axis.ticks.x=element_line(color="black"),
                axis.title.y=element_text(margin=margin(r=15)),
                panel.border=element_rect(color="black", fill=NA),
                plot.margin=margin(.5, .5, .5, .5, "cm"),
                panel.grid.major.x=element_blank(),
                panel.grid.minor.y=element_blank(),
                panel.grid.major.y=element_line(linewidth=.2, color="lightgrey")) +
          scale_fill_manual(values=brewer.pal(3, "Set2")) +
          facet_wrap(.~Cluster, ncol=1) +
          theme(strip.text.x=element_text(size=11, face="bold")) +
          guides(fill=guide_legend(title="Egfr status"), 
                 pattern=guide_legend(title="Background"))
        
        ## Plot together:
        
        grid.arrange(p1, 
                     arrangeGrob(p_seed, p_enrich, ncol=1, heights=c(.5, 1)), 
                     top=grid::textGrob(label=plot_title, gp=grid::gpar(fontsize=13)),
                     ncol=2)
        
        grid.arrange(p2, 
                     arrangeGrob(p_seed, p_enrich, ncol=1, heights=c(.5, 1)), 
                     top=grid::textGrob(label=plot_title, gp=grid::gpar(fontsize=13)),
                     ncol=2)
        
      }

    })
    
  })
  
  dev.off()
  
}