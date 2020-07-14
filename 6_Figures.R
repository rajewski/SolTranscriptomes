# To make the figure for the pre-ripening transcriptome data
library("DEGreport")
library("ggplot2")
library("cowplot")
theme_set(theme_cowplot())
library("topGO")
library("wesanderson")
library("tidyr")
library("dplyr")
library("Rgraphviz")
library("scales")

# Cluster Plots -----------------------------------------------------------
Cluster_AllOrtho_Noise <- readRDS("DEGAnalysis/RNA-seq/Cluster_AllOrtho_Noise.rds")
Cluster_AllOrtho_DEGByFruit <- readRDS("DEGAnalysis/RNA-seq/Cluster_AllOrtho_DEGByFruit.rds")
Cluster_AllOrtho_DEBbySpecies <- readRDS("DEGAnalysis/RNA-seq/Cluster_AllOrtho_DEGBySpecies.rds")

pal <- wes_palette("Royal1", 2, type="discrete")

# Make Model 1 Cluster 1 and 2 plots
M1 <- Cluster_AllOrtho_Noise$normalized

M1C1_Plot <- ggplot(subset(M1, cluster %in% 1),
                    aes(x=Stage, y=value, fill=pal[2], col=pal[2])) +
  labs(y="Z-score"
       #,
       #title="Cluster 1",
       #subtitle = "(n=580)"
       ) +
  scale_fill_manual(values=pal[2]) +
  scale_color_manual(values=pal[2]) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust=0.5),
        legend.position = "none") +
  geom_violin(alpha=0.5) +
  stat_summary(fun=mean, geom="line", color=pal[1], aes(group=1))

M1C2_Plot <- ggplot(subset(M1, cluster %in% 2),
       aes(x=Stage, y=value, fill=pal[2], col=pal[2])) +
  labs(y="Z-score"
       #,
       #title="Cluster 2",
       #subtitle = "(n=199)"
       ) +
  scale_fill_manual(values=pal[2]) +
  scale_color_manual(values=pal[2]) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust=0.5),
        legend.position = "none") +
  geom_violin(alpha=0.5) +
  stat_summary(fun=mean, geom="line", color=pal[1], aes(group=1))
  
# Make Model 2 Cluster 2, 14, and 19
M2 <- Cluster_AllOrtho_DEGByFruit$normalized
M2All_plot <- ggplot(M2, aes(x=Stage, y=value, col=Fruit, fill=Fruit)) +
  labs(y="Z-score") +
  scale_fill_manual(values=pal) +
  scale_color_manual(values=pal) +
  facet_wrap(~cluster,
             labeller=label_both,
             nrow = 3) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust=0.5),
        legend.position = "none",
        strip.background = element_rect(fill="#FFFFFF")) +
  geom_violin(position=position_dodge(width=0), alpha=0.5) +
  stat_summary(fun=mean, geom="line", aes(group=Fruit))


M2C1_Plot <- ggplot(subset(M2, cluster %in% 1),
                    aes(x=Stage, y=value, col=Fruit, fill=Fruit)) +
  labs(y="Z-score"
       #,
       #title="Cluster 1",
       #subtitle="(n=777)"
  ) +
  scale_fill_manual(values=pal) +
  scale_color_manual(values=pal) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust=0.5),
        legend.position = "none") +
  geom_violin(position=position_dodge(width=0), alpha=0.5) +
  stat_summary(fun=mean, geom="line", aes(group=Fruit))

M2C2_Plot <- ggplot(subset(M2, cluster %in% 2),
               aes(x=Stage, y=value, col=Fruit, fill=Fruit)) +
  labs(y="Z-score"
       #,
       #title="Cluster 2",
       #subtitle="(n=268)"
       ) +
  scale_fill_manual(values=pal) +
  scale_color_manual(values=pal) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust=0.5),
        legend.position = "none") +
  geom_violin(position=position_dodge(width=0), alpha=0.5) +
  stat_summary(fun=mean, geom="line", aes(group=Fruit))

M2C5_Plot <- ggplot(subset(M2, cluster %in% 5),
                    aes(x=Stage, y=value, col=Fruit, fill=Fruit)) +
  labs(y="Z-score",
       title="Cluster 5",
       subtitle="(n=292)"
  ) +
  scale_fill_manual(values=pal) +
  scale_color_manual(values=pal) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust=0.5),
        legend.position = "none") +
  geom_violin(position=position_dodge(width=0), alpha=0.5) +
  stat_summary(fun=mean, geom="line", aes(group=Fruit))

M2C6_Plot <- ggplot(subset(M2, cluster %in% 6),
                    aes(x=Stage, y=value, col=Fruit, fill=Fruit)) +
  labs(y="Z-score",
       title="Cluster 6",
       subtitle="(n=64)"
  ) +
  scale_fill_manual(values=pal) +
  scale_color_manual(values=pal) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust=0.5),
        legend.position = "none") +
  geom_violin(position=position_dodge(width=0), alpha=0.5) +
  stat_summary(fun=mean, geom="line", aes(group=Fruit))


M2C14_Plot <- ggplot(subset(M2, cluster %in% 14),
                aes(x=Stage, y=value, col=Fruit, fill=Fruit)) +
  labs(y="Z-score"
       ,
       title="Cluster 14",
       subtitle="(n=18)"
       ) +
  scale_fill_manual(values=pal) +
  scale_color_manual(values=pal) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust=0.5),
        legend.position = "none") +
  geom_violin(position=position_dodge(width=0), alpha=0.5) +
  stat_summary(fun=mean, geom="line", aes(group=Fruit))

M2C18_Plot <- ggplot(subset(M2, cluster %in% 18),
                    aes(x=Stage, y=value, col=Fruit, fill=Fruit)) +
  labs(y="Z-score",
       title="Cluster 18",
       subtitle="(n=19)"
  ) +
  scale_fill_manual(values=pal) +
  scale_color_manual(values=pal) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust=0.5)) +
  geom_violin(position=position_dodge(width=0), alpha=0.5) +
  stat_summary(fun=mean, geom="line", aes(group=Fruit))
  
M2C19_Plot <- ggplot(subset(M2, cluster %in% 19),
                aes(x=Stage, y=value, col=Fruit, fill=Fruit)) +
  labs(y="Z-score"
       ,
       title="Cluster 19",
       subtitle="(n=24)"
       ) +
  scale_fill_manual(values=pal) +
  scale_color_manual(values=pal) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust=0.5),
        legend.position = "none") +
  geom_violin(position=position_dodge(width=0), alpha=0.5) +
  stat_summary(fun=mean, geom="line", aes(group=Fruit))

M2_Legend <- get_legend(ggplot(subset(M2, cluster %in% 19),
                               aes(x=Stage, y=value, col=Fruit, fill=Fruit)) +
                          scale_fill_manual(values=pal) +
                          scale_color_manual(values=pal) +
                          geom_violin(position=position_dodge(width=0), alpha=0.5)+
                          theme(legend.position = "bottom",
                             legend.justification="center",
                             legend.title = element_blank()))


# GO Plots ----------------------------------------------------------------
# grab the functions to make these plots below in the final section

# M1 Overall Table
M1GO <- GOPlot(GOEnrich(gene2go = "DEGAnalysis/Pfam/Ortho.gene2go.tsv",
                   GOIs="DEGAnalysis/RNA-seq/Lists/AllOrtho_Noise_AllGenes.txt"),
               Title = "Overall GO Enrichment") +
  theme(legend.position = "none")
M2GO <- GOPlot(GOEnrich(gene2go = "DEGAnalysis/Pfam/Ortho.gene2go.tsv",
                   GOIs="DEGAnalysis/RNA-seq/Lists/AllOrtho_DEGByFruit_AllGenes.txt"),
               Title = "Overall GO Enrichment") +
  theme(legend.position = "none")
M3GO <- GOPlot(GOEnrich(gene2go = "DEGAnalysis/Pfam/Ortho.gene2go.tsv",
                   GOIs="DEGAnalysis/RNA-seq/Lists/AllOrtho_DEGBySpecies_AllGenes.txt"),
               Title = "Overall GO Enrichment")
M1C1GO <- GOPlot(GOEnrich(gene2go = "DEGAnalysis/Pfam/Ortho.gene2go.tsv",
                      GOIs="DEGAnalysis/RNA-seq/Lists/AllOrtho_Noise_Cluster_1.txt")) +
  theme(legend.position = "none")
M1C2GO <- GOPlot(GOEnrich(gene2go = "DEGAnalysis/Pfam/Ortho.gene2go.tsv",
                      GOIs="DEGAnalysis/RNA-seq/Lists/AllOrtho_Noise_Cluster_2.txt")) +
  theme(legend.position = "none")
M2C1GO <- GOPlot(GOEnrich(gene2go = "DEGAnalysis/Pfam/Ortho.gene2go.tsv",
                          GOIs="DEGAnalysis/RNA-seq/Lists/AllOrtho_DEGByFruit_Cluster_1.txt")) +
  theme(legend.position = "none")
M2C2GO <- GOPlot(GOEnrich(gene2go = "DEGAnalysis/Pfam/Ortho.gene2go.tsv",
                      GOIs="DEGAnalysis/RNA-seq/Lists/AllOrtho_DEGByFruit_Cluster_2.txt")) +
  theme(legend.position = "none")
M2C14GO <- GOPlot(GOEnrich(gene2go = "DEGAnalysis/Pfam/Ortho.gene2go.tsv",
                      GOIs="DEGAnalysis/RNA-seq/Lists/AllOrtho_DEGByFruit_Cluster_14.txt")) +
  theme(legend.position = "none")
M2C19GO <- GOPlot(GOEnrich(gene2go = "DEGAnalysis/Pfam/Ortho.gene2go.tsv",
                       GOIs="DEGAnalysis/RNA-seq/Lists/AllOrtho_DEGByFruit_Cluster_19.txt")) +
  theme(legend.position = "none")

Model1_Title <- ggdraw() + 
  draw_label("Model 1",
             fontface = 'bold',
             hjust = 0.5,
             size=18)
Model2_Title <- ggdraw() + 
  draw_label("Model 2",
             fontface = 'bold',
             hjust = 0.5,
             size=18) +
  theme(plot.margin = margin(0, 0, 0, 7))
Model3_Title <- ggdraw() + 
  draw_label("Model 3",
             fontface = 'bold',
             hjust = 0.5,
             size=18) +
  theme(plot.margin = margin(0, 0, 0, 7))

Cluster1_Title <- ggdraw() + 
  draw_label("Cluster 1",             
             fontface = 'bold',
             hjust = 0.5,
             size=14) +
  theme(plot.margin = margin(0, 0, 0, 7))
Cluster2_Title <- ggdraw() + 
  draw_label("Cluster 2",             
             fontface = 'bold',
             hjust = 0.5,
             size=14) +
  theme(plot.margin = margin(0, 0, 0, 7))
Cluster14_Title <- ggdraw() + 
  draw_label("Cluster 14",             
             fontface = 'bold',
             hjust = 0.5,
             size=14) +
  theme(plot.margin = margin(0, 0, 0, 7))
Cluster19_Title <- ggdraw() + 
  draw_label("Cluster 19",             
             fontface = 'bold',
             hjust = 0.5,
             size=14) +
  theme(plot.margin = margin(0, 0, 0, 7))



# Assemble Figure 2 ---------------------------------------------------

# By row
TopRow <- plot_grid(Model1_Title,
                    Model2_Title,
                    M1GO,
                    M2GO,
                    Cluster1_Title,
                    Cluster1_Title,
                    nrow = 3,
                    labels=c("","","A","G","",""),
                    rel_heights = c(0.1,1,.07),
                    align="h")

MiddleRow1 <- plot_grid(M1C1_Plot,
                        M1C1GO,
                        M2C1_Plot,
                        M2C1GO,
                        nrow=1,
                        rel_widths = c(2,3,2,3),
                        labels=c("B","C","H","I"))
MiddleTitle <- plot_grid(Cluster2_Title,
                         Cluster14_Title,
                         nrow=1)
MiddleRow2 <- plot_grid(M1C2_Plot,
                        M1C2GO,
                        M2C14_Plot,
                        M2C14GO,
                        nrow=1,
                        rel_widths = c(2,3,2,3),
                        labels=c("D","E","J","K"))
BottomTitle <- plot_grid(Model3_Title,
                         Cluster19_Title,
                         nrow=1)
BottomPlots <- plot_grid(M3GO,
                         M2C19_Plot,
                         M2C19GO,
                         nrow=1,
                         rel_widths = c(5,2,3),
                         labels=c("F","L","M"))
BottomRow <- plot_grid(BottomTitle,
                       BottomPlots,
                       nrow=2,
                       rel_heights = c(0.1,1))
FinalPlot <- plot_grid(TopRow,
                       MiddleRow1,
                       MiddleTitle,
                       MiddleRow2,
                       BottomRow,
                       nrow=5,
                       rel_heights = c(1,1,0.1,1,1))

ggsave2("Figures/Figure2.pdf", width=20, height=20)


# Assembly Figure S1 ------------------------------------------------------

Figure_S1 <- plot_grid(M2C5_Plot,
                       M2C6_Plot,
                       M2C18_Plot,
                       nrow=1,
                       labels=c("A","B","C"))
ggsave2("Figures/FigureS1.pdf", width=20, height=6)


# Extenisve GO Plotting Functions -----------------------------------------
# Function to do the enrichment
GOEnrich <- function(gene2go="",
                     GOIs="",
                     GOCategory="BP",
                     NumCategories=20) {
  # Clean the lists of GO terms from the bash script and convert to a named list object
  GO <- read.table(gene2go, stringsAsFactors = F)
  GO <- separate_rows(as.data.frame(GO[,c(1,2)]), 2, sep="\\|")
  GO <- GO %>% 
    distinct() %>% 
    group_by(V1) %>% 
    mutate(V2 = paste0(V2, collapse = "|")) %>%
    distinct()
  GO <- setNames(strsplit(x = unlist(GO$V2), split = "\\|"), GO$V1)
  # Read in the genes of interest to be tested for GO enrichment
  GOI <- scan(file=GOIs,
              what = list(""))[[1]]
  GOI <- factor(as.integer(names(GO) %in% GOI))
  names(GOI) <- names(GO)
  # Create a TopGO object with the necessary data, also get a summary with
  GOData <- new("topGOdata",
                description = "Slyc Cluster 1",
                ontology = GOCategory,
                allGenes = GOI,
                nodeSize = 5,
                annot = annFUN.gene2GO,
                gene2GO=GO)
  # Do the enrichment test
  GOResults <- runTest(GOData,
                       algorithm="weight01",
                       statistic = "fisher")
  # Summarize the test with a table
  GOTable <- GenTable(GOData,
                      Fisher = GOResults,
                      orderBy = "Fisher",
                      ranksOf = "Fisher",
                      topNodes = NumCategories)
  # Summarize with a prettier graph from https://www.biostars.org/p/350710/
  GOTable$Fisher <- gsub("^< ", "", GOTable$Fisher) #remove too low pvals
  GoGraph <- GOTable[as.numeric(GOTable$Fisher)<0.05, c("GO.ID", "Term", "Fisher", "Significant")]
  if(dim(GoGraph)[1]==0){
    return(NULL) #stop empty graphs
  }
  GoGraph$Term <- gsub(" [a-z]*\\.\\.\\.$", "", GoGraph$Term) #clean elipses
  GoGraph$Term <- gsub("\\.\\.\\.$", "", GoGraph$Term)
  GoGraph$Term <- substring(GoGraph$Term,1,26)
  GoGraph$Term <- paste(GoGraph$GO.ID, GoGraph$Term, sep=", ")
  GoGraph$Term <- factor(GoGraph$Term, levels=rev(GoGraph$Term))
  GoGraph$Fisher <- as.numeric(GoGraph$Fisher)
  return(GoGraph)
}

# Function to make a cute plot
GOPlot <- function(GoGraph=X,
                   Title="") 
{
  require(ggplot2)
  GoPlot <- ggplot(GoGraph, aes(x=Term, y=-log10(Fisher), fill=Significant)) +
    stat_summary(geom = "bar", fun = mean, position = "dodge") +
    xlab(element_blank()) +
    ylab("Log Fold Enrichment") +
    scale_fill_gradientn(colours = wes_palette("Royal1", 2, type="discrete"),
                         limits=c(1,130),
                         breaks=c(1,130)) +
    ggtitle(Title) +
    scale_y_continuous(breaks=round(seq(0, max(-log10(GoGraph$Fisher),3)), 1)) +
    #theme_bw(base_size=12) +
    theme(
      panel.grid = element_blank(),
      legend.position=c(0.86,.5),
      legend.background=element_blank(),
      legend.key=element_blank(),     #removes the border
      legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
      #legend.text=element_text(size=18),  #Text size
      legend.title=element_blank(),
      plot.title=element_text(angle=0, face="bold", vjust=1),
      axis.text.x=element_text(angle=0, hjust=0.5),
      axis.text.y=element_text(angle=0, vjust=0.5),
      axis.title=element_text(hjust=0.2),
      #title=element_text(size=18)
      ) +
    guides(fill=guide_colorbar(ticks=FALSE, label.position = 'left')) +
    coord_flip()
  print(GoPlot)
}



  
  