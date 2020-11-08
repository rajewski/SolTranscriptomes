library("ggplot2")
library("ggplotify")
library("cowplot")
theme_set(theme_cowplot())
library("scales")
library("lemon")
library("DESeq2")
source("X_Functions.R")

# Set up ------------------------------------------------------------------
#Trying out the 5 species data, might add into the preripe figure
Cluster_Noise <- readRDS("DEGAnalysis/RNA-seq/Cluster_FiveOrtho_Noise.rds")
Cluster_Fruit <- readRDS("DEGAnalysis/RNA-seq/Cluster_FiveOrtho_Fruit.rds")
Cluster_Species <- readRDS("DEGAnalysis/RNA-seq/Cluster_FiveOrtho_Species.rds")
DDS_Noise <- readRDS("DEGAnalysis/RNA-seq/DDS_FiveOrtho_Noise.rds")
DDS_Fruit <- readRDS("DEGAnalysis/RNA-seq/DDS_FiveOrtho_Fruit.rds")
DDS_Species <- readRDS("DEGAnalysis/RNA-seq/DDS_FiveOrtho_Species.rds")

# Set relabeling for clusters
# Stage Names
Stage_Labs <- c("1"="1", "2"="2", "3"="3", "3.5"="Br", "4"="RR")
# Melon-only
M1_Labs <- paste("Cluster", 1:length(unique(Cluster_Noise$normalized$cluster)))
names(M1_Labs) <- sort(unique(Cluster_Noise$normalized$cluster))
# By fruit
M2_Labs <- paste("Cluster", 1:length(unique(Cluster_Fruit$normalized$cluster)))
names(M2_Labs) <- sort(unique(Cluster_Fruit$normalized$cluster))
# By Species
M3_Labs <- paste("Cluster", 1:length(unique(Cluster_Species$normalized$cluster)))
names(M3_Labs) <- sort(unique(Cluster_Species$normalized$cluster))

# Clusters ----------------------------------------------------------------
# Clusters of 5-species data by fruit type
Model1_plot <- ggplot(Cluster_Noise$normalized, 
                      aes(x=Stage, y=value)) +
  labs(y="Z-score of Expression") +
  facet_rep_wrap(~cluster,
                 labeller=labeller(cluster=M1_Labs),
                 nrow = 1) +
  scale_x_discrete(labels=Stage_Labs) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust=0.5),
        strip.background = element_rect(fill="#FFFFFF"),
        legend.position = c(.95, 0), legend.justification = c(1, -1)) +
  geom_violin(position=position_dodge(width=0), alpha=0.5) +
  stat_summary(fun=mean, geom="line", aes(group=Fruit))
Model1_plot
ggsave("DEGAnalysis/RNA-seq/Plots/ClusterProfiles_FiveOrtho_Noise.pdf", height=7, width=13)

# Clusters of 5-species data by fruit type
Model2_plot <- ggplot(Cluster_Fruit$normalized, 
                      aes(x=Stage, y=value, col=Fruit, fill=Fruit)) +
  labs(y="Z-score of Expression") +
  scale_fill_manual(values=palw[3:4]) +
  scale_color_manual(values=palw[3:4]) +
  facet_rep_wrap(~cluster,
                 labeller=labeller(cluster=M2_Labs),
                 nrow = 3) +
  scale_x_discrete(labels=Stage_Labs) +
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
  scale_x_discrete(labels=Stage_Labs) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust=0.5),
        strip.background = element_rect(fill="#FFFFFF"),
        legend.position = c(.95, 0), legend.justification = c(1, -1)) +
  geom_violin(position=position_dodge(width=0), alpha=0.5) +
  stat_summary(fun=mean, geom="line", aes(group=Species))
Model3_plot
ggsave("DEGAnalysis/RNA-seq/Plots/ClusterProfiles_FiveOrtho_Species.pdf", height=7, width=13)

# GO Tables ---------------------------------------------------------------

# GO Tables for 5-species data noise, all genes
tmpList <- subset(results(DDS_Noise), padj<0.01)
tmpList <- as.factor(as.numeric(rownames(DDS_Noise) %in% rownames(tmpList)))
names(tmpList) <- rownames(DDS_Noise)
FiveOrtho_Noise_Table <- GOEnrich(gene2go = "DEGAnalysis/Pfam/Ortho.831All.gene2go.tsv",GOIs=tmpList)
FiveOrtho_Noise_Plot <- GOPlot(FiveOrtho_Noise_Table, Title="Model 1 GO Enrichment")
capture.output(FiveOrtho_Noise_Table, file="DEGAnalysis/RNA-seq/FiveOrtho_Noise_AllGOTable.txt")

# GO Tables for 5-species data by fruit type, all genes
tmpList <- subset(results(DDS_Fruit), padj<0.01)
tmpList <- as.factor(as.numeric(rownames(DDS_Fruit) %in% rownames(tmpList)))
names(tmpList) <- rownames(DDS_Fruit)
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
tmpList <- subset(results(DDS_Species), padj<0.01)
tmpList <- as.factor(as.numeric(rownames(DDS_Species) %in% rownames(tmpList)))
names(tmpList) <- rownames(DDS_Species)
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
ResSig_Noise <- subset(results(DDS_Noise), padj<=0.01)
Subset_Noise <- DDS_Noise[rownames(DDS_Noise) %in% rownames(ResSig_Noise)]
RLD_Noise <- rlog(Subset_Noise, blind=FALSE)
RLD_Noise$Stage[RLD_Noise$Stage=="3.5"] <- "Br"
RLD_Noise$Stage <- as.factor(RLD_Noise$Stage)
plotPCAmod(RLD_Noise, xPC = 1, yPC=2)

Res_Fruit <- results(DDS_Fruit)
ResSig_Fruit <- subset(Res_Fruit, padj<=0.01)
Subset_Fruit <- DDS_Fruit[rownames(DDS_Fruit) %in% rownames(ResSig_Fruit)]
RLD_Fruit <- rlog(Subset_Fruit, blind=FALSE)
RLD_Fruit$Stage[RLD_Fruit$Stage=="3.5"] <- "Br"
RLD_Fruit$Stage <- as.factor(RLD_Fruit$Stage)
plotPCAmod(RLD_Fruit, xPC=1, yPC=2)

Res_Species <- results(DDS_Species)
ResSig_Species <- subset(Res_Species, padj<=0.01)
Subset_Species <- DDS_Species[rownames(DDS_Species) %in% rownames(ResSig_Species)]
RLD_Species <- rlog(Subset_Species, blind=FALSE)
RLD_Species$Stage[RLD_Species$Stage=="3.5"] <- "Br"
RLD_Species$Stage <- as.factor(RLD_Species$Stage)
plotPCAmod(RLD_Species, xPC=1, yPC=3)



