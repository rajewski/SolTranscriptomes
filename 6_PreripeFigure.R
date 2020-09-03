# To make the figure for the pre-ripening transcriptome data
library("DEGreport")
library("ggplot2")
library("ggplotify")
library("cowplot")
theme_set(theme_cowplot())
library("topGO")
library("tidyr")
library("dplyr")
library("Rgraphviz")
library("scales")
library("lemon")
library("UpSetR")
source("X_Functions.R")

# Read in Data and set global stuff ---------------------------------------
Cluster_Fruit <- readRDS("DEGAnalysis/RNA-seq/Cluster_FiveOrtho_Unripe_DEGByFruit.rds")
Cluster_Species <- readRDS("DEGAnalysis/RNA-seq/Cluster_FiveOrtho_Unripe_DEGBySpecies.rds")
DDS_Noise <- readRDS("DEGAnalysis/RNA-seq/DDS_FiveOrtho_Unripe_Noise.rds")
DDS_Fruit <- readRDS("DEGAnalysis/RNA-seq/DDS_FiveOrtho_Unripe_DEGByFruit.rds")
DDS_Species <- readRDS("DEGAnalysis/RNA-seq/DDS_FiveOrtho_Unripe_DEGBySpecies.rds")

# Set relabeling for clusters
# By fruit
M2_Labs <- paste("Cluster", 1:length(unique(Cluster_Fruit$normalized$cluster)))
names(M2_Labs) <- sort(unique(Cluster_Fruit$normalized$cluster))
# By Species
M3_Labs <- paste("Cluster", 1:length(unique(Cluster_Species$normalized$cluster)))
names(M3_Labs) <- sort(unique(Cluster_Species$normalized$cluster))

# Clusters --------------------------------------------------------------
# Clusters of 5-species data by fruit type
ggplot(Cluster_Fruit$normalized, 
                      aes(x=Stage, y=value, col=Fruit, fill=Fruit)) +
  labs(y="Z-score of Expression") +
  scale_fill_manual(values=palfill[c(1,4)]) +
  scale_color_manual(values=palline[c(1,4)]) +
  facet_rep_wrap(~cluster,
                 labeller=labeller(cluster=M2_Labs),
                 nrow = 2) +
  scale_x_discrete(labels=c("1" = "1", "2" = "2", "3" = "3")) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust=0.5),
        strip.background = element_rect(fill="#FFFFFF"),
        legend.position = c(.95, 0), legend.justification = c(1, -1)) +
  geom_violin(position=position_dodge(width=0), alpha=0.5) +
  stat_summary(fun=mean, geom="line", aes(group=Fruit))
ggsave("DEGAnalysis/RNA-seq/Plots/ClusterProfiles_FiveOrtho_Fruit.pdf", height=7, width=13)

# Clusters of 5-species data by species
ggplot(Cluster_Species$normalized, 
                      aes(x=Stage, y=value, col=Species, fill=Species)) +
  labs(y="Z-score of Expression") +
  facet_rep_wrap(~cluster,
                 labeller=labeller(cluster=M3_Labs),
                 nrow = 3) +
  scale_fill_manual(values=palfill[c(3,7,6,2,5)]) +
  scale_color_manual(values=palline[c(3,7,6,2,5)]) +
  scale_x_discrete(labels=c("1" = "1", "2" = "2", "3" = "3")) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust=0.5),
        strip.background = element_rect(fill="#FFFFFF"),
        legend.position = c(.95, 0), legend.justification = c(1, -1)) +
  geom_violin(position=position_dodge(width=0), alpha=1) +
  stat_summary(fun=mean, geom="line", aes(group=Species))
ggsave("DEGAnalysis/RNA-seq/Plots/ClusterProfiles_FiveOrtho_Species.pdf", height=7, width=13)






Cluster_AllOrtho_Noise <- readRDS("DEGAnalysis/RNA-seq/Cluster_AllOrtho_Noise.rds")
#write.table(unique(Cluster_AllOrtho_Noise$normalized$genes),quote=F, col.names = F, row.names = F, file="DEGAnalysis/RNA-seq/Lists/AllOrtho_Noise.txt")

Cluster_AllOrtho_DEGByFruit <- readRDS("DEGAnalysis/RNA-seq/Cluster_AllOrtho_DEGByFruit.rds")
#write.table(unique(Cluster_AllOrtho_DEGByFruit$normalized$genes),quote=F, col.names = F, row.names = F, file="DEGAnalysis/RNA-seq/Lists/AllOrtho_DEGByFruit.txt")

Cluster_AllOrtho_DEGBySpecies <- readRDS("DEGAnalysis/RNA-seq/Cluster_AllOrtho_DEGBySpecies.rds")
#write.table(unique(Cluster_AllOrtho_DEGBySpecies$normalized$genes),quote=F, col.names = F, row.names = F, file="DEGAnalysis/RNA-seq/Lists/AllOrtho_DEGBySpecies.txt")

# 4 Species, Stage 1-Breaker Plots ----------------------------------------------------------
# Labeller object to rename clusters
M1_Labs <- c("1" = "Cluster 1", "2" = "Cluster 2", "3" = "Cluster 3", "5"="Cluster 4")
Model1_plot <- ggplot(Cluster_AllOrtho_Noise$normalized, 
                            aes(x=Stage, y=value, col=Species, fill=Species)) +
  labs(y="Z-score of Expression") +
  scale_fill_manual(values="#9B9B9B") +
  scale_color_manual(values="#000000") +
  facet_rep_wrap(~cluster,
                 labeller=labeller(cluster=M1_Labs),
                 nrow = 2) +
  scale_x_discrete(labels=c("1" = "1", "2" = "2", "3" = "3", "3.5" = "Br")) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust=0.5),
        strip.background = element_rect(fill="#FFFFFF"),
        legend.position = "none") +
  geom_violin(position=position_dodge(width=0), alpha=0.5) +
  stat_summary(fun=mean, geom="line", aes(group=Species))
Model1_plot
ggsave2("DEGAnalysis/RNA-seq/Plots/ClusterProfiles_AllOrtho_Noise.pdf", height=7, width=13)

M2_Labs <- paste("Cluster", 1:length(unique(Cluster_AllOrtho_DEGByFruit$normalized$cluster)))
names(M2_Labs) <- sort(unique(Cluster_AllOrtho_DEGByFruit$normalized$cluster))
Model2_plot <- ggplot(Cluster_AllOrtho_DEGByFruit$normalized, 
                     aes(x=Stage, y=value, col=Fruit, fill=Fruit)) +
  labs(y="Z-score of Expression") +
  scale_fill_manual(values=palw[3:4]) +
  scale_color_manual(values=palw[3:4]) +
  facet_rep_wrap(~cluster,
                 labeller=labeller(cluster=M2_Labs),
                 nrow = 3) +
  scale_x_discrete(labels=c("1" = "1", "2" = "2", "3" = "3", "3.5" = "Br")) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust=0.5),
        strip.background = element_rect(fill="#FFFFFF"),
        legend.position = c(.95, 0), legend.justification = c(1, -1)) +
  geom_violin(position=position_dodge(width=0), alpha=0.5) +
  stat_summary(fun=mean, geom="line", aes(group=Fruit))
Model2_plot
ggsave("DEGAnalysis/RNA-seq/Plots/ClusterProfiles_AllOrtho_Fruit.pdf", height=7, width=13)

M3_Labs <- paste("Cluster", 1:length(unique(Cluster_AllOrtho_DEGBySpecies$normalized$cluster)))
names(M3_Labs) <- sort(unique(Cluster_AllOrtho_DEGBySpecies$normalized$cluster))
Model3_plot <- ggplot(Cluster_AllOrtho_DEGBySpecies$normalized, 
                     aes(x=Stage, y=value, col=Species, fill=Species)) +
  labs(y="Z-score of Expression") +
  scale_fill_manual(values=palw[c(2,1,3,4)]) +
  scale_color_manual(values=palw[c(2,1,3,4)]) +
  facet_rep_wrap(~cluster,
                 labeller=labeller(cluster=M3_Labs),
                 nrow = 5) +
  scale_x_discrete(labels=c("1" = "1", "2" = "2", "3" = "3", "3.5" = "Br")) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust=0.5),
        strip.background = element_rect(fill="#FFFFFF"),
        legend.position = c(1, 0), legend.justification = c(1.5, -0.3)) +
  geom_violin(position=position_dodge(width=0), alpha=0.5) +
  stat_summary(fun=mean, geom="line", aes(group=Species))
Model3Legend <- get_legend(Model3_plot + theme(legend.position = c(0.5,0.8)))
Model3Legend_Alone <- get_legend(Model3_plot +theme(legend.position=c(1.2,0.8)))
Model3_plot
ggsave2("DEGAnalysis/RNA-seq/Plots/ClusterProfiles_AllOrtho_Species.pdf", height=7, width=13)

# Broken up plot to wrap around GO plot
Model3_BrokenPlot_1 <- ggplot(subset(Cluster_AllOrtho_DEGBySpecies$normalized, 
                             cluster %in% c(1,4,7,8,10,11,13,18,19)), 
                      aes(x=Stage, y=value, col=Species, fill=Species)) +
  labs(y=NULL, x=NULL) +
  scale_fill_manual(values=palw[c(2,1,3,4)]) +
  scale_color_manual(values=palw[c(2,1,3,4)]) +
  facet_rep_wrap(~cluster,
                 labeller=labeller(cluster=M3_Labs),
                 nrow = 3) +
  scale_x_discrete(labels=c("1" = "1", "2" = "2", "3" = "3", "3.5" = "Br")) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust=0.5),
        strip.background = element_rect(fill="#FFFFFF"),
        legend.position = "none",
        axis.text.x=element_blank()) +
  geom_violin(position=position_dodge(width=0), alpha=0.5) +
  stat_summary(fun=mean, geom="line", aes(group=Species))

Model3_BrokenPlot_2 <- ggplot(subset(Cluster_AllOrtho_DEGBySpecies$normalized, 
                                     cluster %notin% c(1,4,7,8,10,11,13,18,19)), 
                              aes(x=Stage, y=value, col=Species, fill=Species)) +
  labs(y="Z-score of Expression") +
  scale_fill_manual(values=palw[c(2,1,3,4)]) +
  scale_color_manual(values=palw[c(2,1,3,4)]) +
  facet_rep_wrap(~cluster,
                 labeller=labeller(cluster=M3_Labs),
                 nrow = 3) +
  scale_x_discrete(labels=c("1" = "1", "2" = "2", "3" = "3", "3.5" = "Br")) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust=0.5),
        strip.background = element_rect(fill="#FFFFFF"),
        legend.position = "none") +
  geom_violin(position=position_dodge(width=0), alpha=0.5) +
  stat_summary(fun=mean, geom="line", aes(group=Species))

# Expression Range by Cluster ---------------------------------------------
M2Aggr <- aggregate(Cluster_AllOrtho_DEGByFruit$normalized,
                    by=list(Cluster_AllOrtho_DEGByFruit$normalized$Fruit,
                            Cluster_AllOrtho_DEGByFruit$normalized$cluster,
                            Cluster_AllOrtho_DEGByFruit$normalized$Stage), FUN=mean)[,c(1,2,3,6)]
M2Aggr <- M2Aggr %>% spread(Group.3, value)
M2Aggr$Range <- apply(M2Aggr[,3:6], 1, FUN=max) - apply(M2Aggr[,3:6], 1, FUN=min)
M2Flx <- M2Aggr[M2Aggr$Range>=1,]
M3Aggr <- aggregate(Cluster_AllOrtho_DEGBySpecies$normalized,
                    by=list(Cluster_AllOrtho_DEGBySpecies$normalized$Species,
                            Cluster_AllOrtho_DEGBySpecies$normalized$cluster,
                            Cluster_AllOrtho_DEGBySpecies$normalized$Stage), FUN=mean)[,c(1,2,3,6)]
M3Aggr <- M3Aggr %>% spread(Group.3, value)
M3Aggr$Range <- apply(M3Aggr[,3:6], 1, FUN=max) - apply(M3Aggr[,3:6], 1, FUN=min)
M3Flx <- M3Aggr[M3Aggr$Range>=1,]
M3_Labs[names(M3_Labs) %in% M3Flx$Group.2[M3Flx$Group.1=="Tobacco"]]

# Orthogene Overlap among Models ------------------------------------------
M1orthos <- list(as.data.frame(unique(Cluster_AllOrtho_Noise$normalized$genes)))
M2orthos <- list(as.data.frame(unique(Cluster_AllOrtho_DEGByFruit$normalized$genes)))
M3orthos <- list(as.data.frame(unique(Cluster_AllOrtho_DEGBySpecies$normalized$genes)))
List_Models <- Map(list,
                   M1orthos, M2orthos, M3orthos)
Gene_List <- as.data.frame(unique(unlist(List_Models)))
for (GeneTable in c(M1orthos, M2orthos, M3orthos)) {
  Common<-as.data.frame(Gene_List[,1] %in% as.data.frame(GeneTable)[,1])
  Common[,1][Common[,1]==TRUE]<- 1
  Gene_List<-cbind(Gene_List,Common)
}
colnames(Gene_List) <- c("ID","Model 1", "Model 2", "Model 3")
OverlapUpset <- upset(Gene_List[c(1,4,3,2)], 
      sets = colnames(Gene_List)[c(4,3,2)],
      mainbar.y.label = "Number of DEGs", 
      sets.x.label = NULL, 
      keep.order = T,
      text.scale = 1.5,
      point.size = 3, 
      line.size = 1,
      nintersects = 100)

#Get the genes present in all models
AlwaysDE <- as.character(Gene_List$ID[apply(Gene_List[,c(2:4)],1, FUN=min)==1])
orthogenes <- read.table("Orthofinder/OrthoFinder/Results_May17/Orthogroups/Orthogroups.tsv", stringsAsFactors = F, sep="\t", header = T)
AlwaysDE <- orthogenes[orthogenes$Orthogroup %in% AlwaysDE,]
write.table(AlwaysDE$Arabidopsis, file="DEGAnalysis/RNA-seq/Lists/SevenOrthos.tsv", row.names = F, col.names = F, quote = F)

# 4 Species, Stage 1-Breaker GO Plots ----------------------------------------------
# Make a list of the GO Tables
Model1Tables <- list()
for (i in levels(as.factor(Cluster_AllOrtho_Noise$normalized$cluster))) {
  tmpList <- as.factor(as.numeric(Cluster_AllOrtho_Noise$normalized$genes %in% Cluster_AllOrtho_Noise$normalized$genes[Cluster_AllOrtho_Noise$normalized$cluster==i]))
  names(tmpList) <- Cluster_AllOrtho_Noise$normalized$genes
  Model1Tables[[i]] <- GOEnrich(gene2go = "DEGAnalysis/Pfam/Ortho.gene2go.tsv",
                                  GOIs=tmpList)
}
capture.output(Model1Tables, file="DEGAnalysis/RNA-seq/AllOrtho_Noise_GOTables.txt")
# All the genes as a cohort
Model1_AllTable <- GOEnrich(gene2go = "DEGAnalysis/Pfam/Ortho.gene2go.tsv",
         GOIs="DEGAnalysis/RNA-seq/Lists/AllOrtho_Noise.txt")
capture.output(Model1_AllTable, file="DEGAnalysis/RNA-seq/AllOrtho_Noise_AllGOTable.txt")
Model1_AllPlot <- GOPlot(Model1_AllTable, Title="Overall GO", LegendLimit = 36) +
  theme(legend.position = "none")

Model2Tables <- list()
for (i in levels(as.factor(Cluster_AllOrtho_DEGByFruit$normalized$cluster))) {
  tmpList <- as.factor(as.numeric(Cluster_AllOrtho_DEGByFruit$normalized$genes %in% Cluster_AllOrtho_DEGByFruit$normalized$genes[Cluster_AllOrtho_DEGByFruit$normalized$cluster==i]))
  names(tmpList) <- Cluster_AllOrtho_DEGByFruit$normalized$genes
  Model2Tables[[i]] <- GOEnrich(gene2go = "DEGAnalysis/Pfam/Ortho.gene2go.tsv",
                                GOIs=tmpList)
}
names(Model2Tables) <- M2_Labs[names(Model2Tables)]
capture.output(Model2Tables, file="DEGAnalysis/RNA-seq/AllOrtho_DEGByFruit_GOTables.txt")

M2C4_Plot <- GOPlot(Model2Tables[['Cluster 4']], Title = "Cluster 4 GO", LegendLimit = 36) +
  theme(legend.position = "none")
M2C12_Plot <- GOPlot(Model2Tables[['Cluster 11']], Title = "Cluster 11 GO", LegendLimit = 36) +
  theme(legend.position = "none")
M2C18_Plot <- GOPlot(Model2Tables[['Cluster 16']], Title = "Cluster 16 GO", LegendLimit = 36) +
  theme(legend.position = "none")

# All the genes as a cohort
Model2_AllTable <- GOEnrich(gene2go = "DEGAnalysis/Pfam/Ortho.gene2go.tsv",
                            GOIs="DEGAnalysis/RNA-seq/Lists/AllOrtho_DEGByFruit.txt")
Model2_AllPlot <- GOPlot(Model2_AllTable, Title="Overall GO", LegendLimit = 36) +
  theme(legend.position="none")
capture.output(Model2_AllTable, file="DEGAnalysis/RNA-seq/AllOrtho_DEGByFruit_AllGOTable.txt")

Model3Tables <- list()
for (i in levels(as.factor(Cluster_AllOrtho_DEGBySpecies$normalized$cluster))) {
  tmpList <- as.factor(as.numeric(Cluster_AllOrtho_DEGBySpecies$normalized$genes %in% Cluster_AllOrtho_DEGBySpecies$normalized$genes[Cluster_AllOrtho_DEGByFruit$normalized$cluster==i]))
  if (length(levels(tmpList))==1) {
    next
  }
  names(tmpList) <- Cluster_AllOrtho_DEGBySpecies$normalized$genes
  Model3Tables[[i]] <- GOEnrich(gene2go = "DEGAnalysis/Pfam/Ortho.gene2go.tsv",
                                GOIs=tmpList)
}
names(Model3Tables) <- M3_Labs[names(Model3Tables)]
capture.output(Model3Tables, file="DEGAnalysis/RNA-seq/AllOrtho_DEGBySpecies_GOTables.txt")
# All the genes as a cohort
Model3_AllTable <- GOEnrich(gene2go = "DEGAnalysis/Pfam/Ortho.gene2go.tsv",
                            GOIs="DEGAnalysis/RNA-seq/Lists/AllOrtho_DEGBySpecies.txt",
                            NumCategories = 50)
Model3_AllPlot <- GOPlot(Model3_AllTable, Title="Overall GO", LegendLimit = 36) +
  theme(legend.position=c(0.86,0.3))
GO_Legend <- get_legend(Model3_AllPlot + theme(legend.position=c(0,0)))
capture.output(Model3_AllTable, file="DEGAnalysis/RNA-seq/AllOrtho_DEGBySpecies_AllGOTable.txt")

# The seven genes in all models
SevenTables <- GOEnrich(gene2go = "DEGAnalysis/Pfam/TAIR10.gene2go.tsv",
                              GOIs="DEGAnalysis/RNA-seq/Lists/SevenOrthos.tsv")


# Assemble the Figure -----------------------------------------------------

#Make Titles
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
Overlap_Title <- ggdraw() + 
  draw_label("All Models",
             fontface = 'bold',
             hjust = 0.5,
             size=18) +
  theme(plot.margin = margin(0, 0, 0, 7))

#Model 1 Panel
Panel_M1 <- plot_grid(Model1_Title,
                      Model1_AllPlot,
                      Model1_plot,
                      nrow = 3,
                      labels=c("","A", "B"),
                      rel_heights = c(0.1,.4,1),
                      axis="r",
                      align="h")
# Save Panel A as its own figure
#Panel_M1
#ggsave2("Figures/Figure 2 Model 1.pdf", height=10, width=11)

Panel_M2R1 <- plot_grid(Model2_AllPlot,
                      M2C4_Plot,
                      M2C12_Plot,
                      M2C18_Plot,
                      nrow = 1,
                      labels=c("C","E", "F", "G"), #For Model 2 with others
                      #labels=c("A","C","D","E"), # For Model 2 alone
                      axis="r",
                      align="h")
Panel_M2 <- plot_grid(Model2_Title,
                      Panel_M2R1,
                      Model2_plot,
                      nrow = 3,
                      labels=c("","", "D"), #For Model 2 with others
                      #labels=c("", "", "B"), # For Model 2 alone
                      rel_heights = c(0.1,1,2),
                      axis="r",
                      align="h")
# Save Model 2 as its own figure
Panel_M2
ggsave2("Figures/Figure 2 Model 2.pdf", height=11, width=22)

Panel_Top <- plot_grid(Panel_M1,
                       Panel_M2,
                       nrow=1,
                       rel_widths = c(1,2))
Panel_M3_1 <- plot_grid(Model3_AllPlot,
                      Model3_BrokenPlot_1,
                      ncol=2,
                      rel_widths = c(3.9,3),
                      labels=c("A", "B")) # For Model 3 alone
                      #labels=c("H", "I")) # For Model 3 with others
Panel_M3_2 <- plot_grid(Model3_Title,
                        Panel_M3_1,
                        Model3_BrokenPlot_2,
                        nrow=3,
                        rel_heights = c(0.1,1,1))

# Save Model 3 as its own figure
# Panel_M3_Alone <- plot_grid(Panel_M3_2,
#                             Model3Legend_Alone,
#                             ncol=2,
#                             rel_widths = c(11.6,1))
# Panel_M3_Alone
# ggsave2("Figures/Figure 2 Model 3.pdf", width=24, height=11)

Panel_Legend <- plot_grid(Model3Legend,
                          align="h",
                          axis="l")
Panel_Upset1 <- plot_grid(Panel_Legend,
                         OverlapUpset$Main_bar, 
                         OverlapUpset$Sizes, 
                         OverlapUpset$Matrix,
                         labels=c("","J","", ""),
                         nrow=2, 
                         align='v', 
                         axis="t",
                         rel_heights = c(3,1),
                         rel_widths = c(2,3))
Panel_Upset2 <- plot_grid(Overlap_Title,
                         Panel_Upset1,
                         nrow=2,
                         rel_heights = c(0.1,1))
# Save Upset Plot as its own Figure
Panel_Upset1_Alone <- plot_grid(NULL,
                                OverlapUpset$Main_bar, 
                                OverlapUpset$Sizes, 
                                OverlapUpset$Matrix,
                                nrow=2, 
                                align='v', 
                                axis="t",
                                rel_heights = c(3,1),
                                rel_widths = c(2,3))
Panel_Upset2_Alone <- plot_grid(Overlap_Title,
                                Panel_Upset1_Alone,
                                nrow=2,
                                rel_heights = c(0.1,1))
Panel_Upset2_Alone
ggsave2("Figures/Figure 2 All Models.pdf", height=6, width=10)
Panel_Bottom <- plot_grid(Panel_M3_2,
                          Panel_Upset2,
                          ncol = 2,
                          rel_widths = c(2,1))
FullFigure <- plot_grid(Panel_Top,
                        Panel_Bottom,
                        nrow = 2)
FullFigure
ggsave2(file="Figures/Preripe.pdf",plot=FullFigure, width=32, height=22)

  
  