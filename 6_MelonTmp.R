library("ggplot2")
library("ggplotify")
library("cowplot")
theme_set(theme_cowplot())
library("scales")
library("lemon")
source("X_Functions.R")

#Trying out the 5 species data, might add into the preripe figure
Cluster_Melon <- readRDS("DEGAnalysis/RNA-seq/Cluster_Melon.rds")
Cluster_Fruit <- readRDS("DEGAnalysis/RNA-seq/Cluster_FiveOrtho_Unripe_DEGByFruit.rds")
Cluster_Species <- readRDS("DEGAnalysis/RNA-seq/Cluster_FiveOrtho_Unripe_DEGBySpecies.rds")
DDS_FiveOrtho_Noise <- readRDS("DEGAnalysis/RNA-seq/DDS_FiveOrtho_Unripe_Noise.rds")
DDS_FiveOrtho_Fruit <- readRDS("DEGAnalysis/RNA-seq/DDS_FiveOrtho_Unripe_DEGByFruit.rds")
DDS_FiveOrtho_Species <- readRDS("DEGAnalysis/RNA-seq/DDS_FiveOrtho_Unripe_DEGBySpecies.rds")

# Set relabeling for clusters
# Melon-only
Melon_Labs <- paste("Cluster", 1:length(unique(Cluster_Melon$normalized$cluster)))
names(Melon_Labs) <- sort(unique(Cluster_Melon$normalized$cluster))
# By fruit
M2_Labs <- paste("Cluster", 1:length(unique(Cluster_Fruit$normalized$cluster)))
names(M2_Labs) <- sort(unique(Cluster_Fruit$normalized$cluster))
# By Species
M3_Labs <- paste("Cluster", 1:length(unique(Cluster_Species$normalized$cluster)))
names(M3_Labs) <- sort(unique(Cluster_Species$normalized$cluster))

# Clusters ----------------------------------------------------------------
# Clusters of Melon-only data
Melon_plot <- ggplot(Cluster_Melon$normalized, 
                      aes(x=DAP, y=value, col=Fruit, fill=Fruit)) +
  labs(y="Z-score of Expression") +
  scale_fill_manual(values=palw[3]) +
  scale_color_manual(values=palw[3]) +
  facet_rep_wrap(~cluster,
                 labeller=labeller(cluster=Melon_Labs),
                 nrow = 3) +
  scale_x_discrete(labels=c("1" = "1", "2" = "2", "3" = "3")) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust=0.5),
        strip.background = element_rect(fill="#FFFFFF"),
        legend.position = c(.95, 0), legend.justification = c(1, -1)) +
  geom_violin(position=position_dodge(width=0), alpha=0.5) +
  stat_summary(fun=mean, geom="line", aes(group=Fruit))
Melon_plot
ggsave("DEGAnalysis/RNA-seq/Plots/ClusterProfiles_Melon.pdf", height=7, width=13)

# Clusters of 5-species data by fruit type
Model2_plot <- ggplot(Cluster_Fruit$normalized, 
                      aes(x=Stage, y=value, col=Fruit, fill=Fruit)) +
  labs(y="Z-score of Expression") +
  scale_fill_manual(values=palw[3:4]) +
  scale_color_manual(values=palw[3:4]) +
  facet_rep_wrap(~cluster,
                 labeller=labeller(cluster=M2_Labs),
                 nrow = 3) +
  scale_x_discrete(labels=c("1" = "1", "2" = "2", "3" = "3")) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust=0.5),
        strip.background = element_rect(fill="#FFFFFF"),
        legend.position = c(.95, 0), legend.justification = c(1, -1)) +
  geom_violin(position=position_dodge(width=0), alpha=0.5) +
  stat_summary(fun=mean, geom="line", aes(group=Fruit))
Model2_plot
ggsave("DEGAnalysis/RNA-seq/Plots/ClusterProfiles_FiveOrtho_Fruit.pdf", height=7, width=13)

# Clusters of 5-species data by species
Model3_plot <- ggplot(Cluster_Species$normalized, 
                      aes(x=Stage, y=value, col=Species, fill=Species)) +
  labs(y="Z-score of Expression") +
  facet_rep_wrap(~cluster,
                 labeller=labeller(cluster=M3_Labs),
                 nrow = 3) +
  scale_x_discrete(labels=c("1" = "1", "2" = "2", "3" = "3")) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust=0.5),
        strip.background = element_rect(fill="#FFFFFF"),
        legend.position = c(.95, 0), legend.justification = c(1, -1)) +
  geom_violin(position=position_dodge(width=0), alpha=0.5) +
  stat_summary(fun=mean, geom="line", aes(group=Species))
Model3_plot
ggsave("DEGAnalysis/RNA-seq/Plots/ClusterProfiles_FiveOrtho_Species.pdf", height=7, width=13)

# GO Tables ---------------------------------------------------------------
# GO tables for melon-only data
MelonTables <- list()
for (i in levels(as.factor(Cluster_Melon$normalized$cluster))) {
  tmpList <- as.factor(as.numeric(Cluster_Melon$normalized$genes %in% Cluster_Melon$normalized$genes[Cluster_Melon$normalized$cluster==i]))
  names(tmpList) <- strtrim(Cluster_Melon$normalized$genes,12)
  
  MelonTables[[i]] <- GOEnrich(gene2go = "DEGAnalysis/Pfam/Melon.gene2go.tsv",
                                GOIs=tmpList)
}
names(MelonTables) <- Melon_Labs[names(MelonTables)]
capture.output(MelonTables, file="DEGAnalysis/RNA-seq/Melon_GOTables.txt")

# GO Tables for 5-species data noise, all genes
tmpList <- results(DDS_FiveOrtho_Noise)
tmpList <- subset(tmpList, padj<0.01)
tmpList <- as.factor(as.numeric(rownames(DDS_FiveOrtho_Noise) %in% rownames(tmpList)))
names(tmpList) <- rownames(DDS_FiveOrtho_Noise)
FiveOrtho_Noise_Table <- GOEnrich(gene2go = "DEGAnalysis/Pfam/Ortho.831All.gene2go.tsv",GOIs=tmpList)
FiveOrtho_Noise_Plot <- GOPlot(FiveOrtho_Noise_Table, Title="Overall GO")
capture.output(FiveOrtho_Noise_Table, file="DEGAnalysis/RNA-seq/FiveOrtho_Noise_AllGOTable.txt")

# GO Tables for 5-species data by fruit type, all genes
tmpList <- results(DDS_FiveOrtho_Fruit)
tmpList <- subset(tmpList, padj<0.01)
tmpList <- as.factor(as.numeric(rownames(DDS_FiveOrtho_Fruit) %in% rownames(tmpList)))
names(tmpList) <- rownames(DDS_FiveOrtho_Fruit)
FiveOrtho_Fruit_Table <- GOEnrich(gene2go = "DEGAnalysis/Pfam/Ortho.831All.gene2go.tsv",
                                  GOIs=tmpList)
FiveOrtho_Fruit_Plot <- GOPlot(FiveOrtho_Fruit_Table, Title="Overall GO")
capture.output(FiveOrtho_Fruit_Table, file="DEGAnalysis/RNA-seq/FiveOrtho_Fruit_AllGOTable.txt")

# GO tables for 5-species data by fruit type, clustered
Model2Tables <- list()
for (i in levels(as.factor(Cluster_Fruit$normalized$cluster))) {
  tmpList <- as.factor(as.numeric(Cluster_Fruit$normalized$genes %in% Cluster_Fruit$normalized$genes[Cluster_Fruit$normalized$cluster==i]))
  names(tmpList) <- Cluster_Fruit$normalized$genes
  Model2Tables[[i]] <- GOEnrich(gene2go = "DEGAnalysis/Pfam/Ortho.831All.gene2go.tsv",
                                GOIs=tmpList)
}
names(Model2Tables) <- M2_Labs[names(Model2Tables)]
capture.output(Model2Tables, file="DEGAnalysis/RNA-seq/Five_Fruit_GOTables.txt")

# GO Tables for 5-species data by species, all genes
tmpList <- results(DDS_FiveOrtho_Species)
tmpList <- subset(tmpList, padj<0.01)
tmpList <- as.factor(as.numeric(rownames(DDS_FiveOrtho_Species) %in% rownames(tmpList)))
names(tmpList) <- rownames(DDS_FiveOrtho_Species)
FiveOrtho_Species_Table <- GOEnrich(gene2go = "DEGAnalysis/Pfam/Ortho.831All.gene2go.tsv",GOIs=tmpList) # Is a null set

# GO tables for 5-species data by species, clustered
Model3Tables <- list()
for (i in levels(as.factor(Cluster_Species$normalized$cluster))) {
  tmpList <- as.factor(as.numeric(Cluster_Species$normalized$genes %in% Cluster_Species$normalized$genes[Cluster_Species$normalized$cluster==i]))
  names(tmpList) <- Cluster_Species$normalized$genes
  Model3Tables[[i]] <- GOEnrich(gene2go = "DEGAnalysis/Pfam/Ortho.831All.gene2go.tsv",
                                GOIs=tmpList)
}
names(Model3Tables) <- M3_Labs[names(Model3Tables)]
capture.output(Model3Tables, file="DEGAnalysis/RNA-seq/Five_Species_GOTables.txt")


# PCA ---------------------------------------------------------------------
RLD_FiveOrtho_Noise <- rlog(DDS_FiveOrtho_Noise, blind=FALSE)
plotPCA(RLD_FiveOrtho_Noise, intgroup = c("Species"))

RLD_FiveOrtho_Fruit <- rlog(DDS_FiveOrtho_Fruit, blind=FALSE)
plotPCA(RLD_FiveOrtho_Fruit, intgroup = c("Species"))

RLD_FiveOrtho_Species <- rlog(DDS_FiveOrtho_Species, blind=FALSE)
plotPCA(RLD_FiveOrtho_Species, intgroup = c("Species"))
