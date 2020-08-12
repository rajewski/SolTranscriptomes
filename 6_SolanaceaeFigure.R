# To make the figure for the pre-ripening transcriptome data
library("DEGreport")
library("ggplot2")
library("cowplot")
theme_set(theme_cowplot())
library("topGO")
library("wesanderson")
library("viridis")
library("tidyr")
library("dplyr")
library("Rgraphviz")
library("scales")
library("lemon")
source("X_Functions.R")

# Full Model All Cluster Plots -----------------------------------------------------------
Cluster_SolOrtho_Noise <- readRDS("DEGAnalysis/RNA-seq/Cluster_SolOrtho_Noise.rds")
Cluster_SolOrtho_DEGByFruit <- readRDS("DEGAnalysis/RNA-seq/Cluster_SolOrtho_DEGByFruit.rds")
Cluster_SolOrtho_DEGBySpecies <- readRDS("DEGAnalysis/RNA-seq/Cluster_SolOrtho_DEGBySpecies.rds")

pal <- viridis(3, option="D")

Model1_plot <- ggplot(Cluster_SolOrtho_Noise$normalized, 
                            aes(x=Stage, y=value, col=Species, fill=Species)) +
  labs(y="Z-score") +
  scale_fill_manual(values=pal[2]) +
  scale_color_manual(values=pal[2]) +
  facet_rep_wrap(~cluster,
                 nrow = 2) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust=0.5),
        strip.background = element_rect(fill="#FFFFFF"),
        legend.position = "none") +
  geom_violin(position=position_dodge(width=0), alpha=0.5) +
  stat_summary(fun=mean, geom="line", aes(group=Species))

Model2_plot <- ggplot(Cluster_SolOrtho_DEGByFruit$normalized, 
                     aes(x=Stage, y=value, col=Fruit, fill=Fruit)) +
  labs(y="Z-score") +
  scale_fill_manual(values=pal[1:2]) +
  scale_color_manual(values=pal[1:2]) +
  facet_rep_wrap(~cluster,
                 nrow = 4) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust=0.5),
        strip.background = element_rect(fill="#FFFFFF")) +
  geom_violin(position=position_dodge(width=0), alpha=0.5) +
  stat_summary(fun=mean, geom="line", aes(group=Fruit))

Model3_plot <- ggplot(Cluster_SolOrtho_DEGBySpecies$normalized, 
                     aes(x=Stage, y=value, col=Species, fill=Species)) +
  labs(y="Z-score") +
  scale_fill_manual(values=pal) +
  scale_color_manual(values=pal) +
  facet_rep_wrap(~cluster,
                 nrow = 4) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust=0.5),
        strip.background = element_rect(fill="#FFFFFF")) +
  geom_violin(position=position_dodge(width=0), alpha=0.5) +
  stat_summary(fun=mean, geom="line", aes(group=Species))


# GO Tables to Pick Clusters ----------------------------------------------

# Make a list of the GO Tables
Model1Tables <- list()
for (i in levels(as.factor(Cluster_SolOrtho_Noise$normalized$cluster))) {
  tmpList <- as.factor(as.numeric(Cluster_SolOrtho_Noise$normalized$genes %in% Cluster_SolOrtho_Noise$normalized$genes[Cluster_SolOrtho_Noise$normalized$cluster==i]))
  names(tmpList) <- Cluster_SolOrtho_Noise$normalized$genes
  Model1Tables[[i]] <- GOEnrich(gene2go = "DEGAnalysis/Pfam/Ortho.gene2go.tsv",
                                  GOIs=tmpList)
}
capture.output(Model1Tables, file="DEGAnalysis/RNA-seq/SolOrtho_Noise_GOTables.txt")

Model2Tables <- list()
for (i in levels(as.factor(Cluster_SolOrtho_DEGByFruit$normalized$cluster))) {
  tmpList <- as.factor(as.numeric(Cluster_SolOrtho_DEGByFruit$normalized$genes %in% Cluster_SolOrtho_DEGByFruit$normalized$genes[Cluster_SolOrtho_DEGByFruit$normalized$cluster==i]))
  names(tmpList) <- Cluster_SolOrtho_DEGByFruit$normalized$genes
  Model2Tables[[i]] <- GOEnrich(gene2go = "DEGAnalysis/Pfam/Ortho.gene2go.tsv",
                                GOIs=tmpList)
}
capture.output(Model2Tables, file="DEGAnalysis/RNA-seq/SolOrtho_DEGByFruit_GOTables.txt")

Model3Tables <- list()
for (i in levels(as.factor(Cluster_SolOrtho_DEGBySpecies$normalized$cluster))) {
  tmpList <- as.factor(as.numeric(Cluster_SolOrtho_DEGBySpecies$normalized$genes %in% Cluster_SolOrtho_DEGBySpecies$normalized$genes[Cluster_SolOrtho_DEGByFruit$normalized$cluster==i]))
  if (length(levels(tmpList))==1) {
    next
  }
  names(tmpList) <- Cluster_SolOrtho_DEGBySpecies$normalized$genes
  Model3Tables[[i]] <- GOEnrich(gene2go = "DEGAnalysis/Pfam/Ortho.gene2go.tsv",
                                GOIs=tmpList)
}
capture.output(Model3Tables, file="DEGAnalysis/RNA-seq/SolOrtho_DEGBySpecies_GOTables.txt")

# Needs update after cluster re-org ---------------------------------------


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




  
  