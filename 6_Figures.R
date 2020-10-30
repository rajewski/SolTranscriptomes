source("X_Functions.R")
library("ggplot2")
library("cowplot")
theme_set(theme_cowplot())
library("lemon")
library("patchwork")


# Read Data In ------------------------------------------------------------
Cluster_Nobt <- readRDS("DEGAnalysis/RNA-seq/Cluster_Nobt.rds")
Labs_NobtCluster <- ClusterLabs(Cluster_Nobt)
Labs_NobtStage <- c("0"="1", "3"="2", "6"="3", "11"="Br")

Cluster_Solanum <- readRDS("DEGAnalysis/RNA-seq/Cluster_Solanum_3DF.rds")
Labs_SolanumCluster <- ClusterLabs(Cluster_Solanum)
Cluster_Solanum_Noise <-readRDS("DEGAnalysis/RNA-seq/Cluster_Solanum_3DF_Noise.rds")
Labs_SolanumNoiseCluster <- ClusterLabs(Cluster_Solanum_Noise)
Labs_SolanumStage <- c("1"="1", "3"="2", "15"="3", "35"="Br", "45"="RR")

# Tobacco Clusters -----------------------------------------------------------
# consider using lapply to make a list of the individual (unfaceted) plots?
C1 <- ggplot(Cluster_Nobt$normalized, 
                      aes(x=DAP, y=value, col=Species, fill=Species)) +
  labs(y="Z-score of Expression",
       x="Stage") +
  facet_rep_wrap(~cluster,
                 labeller=labeller(cluster=Labs_NobtCluster),
                 nrow = 2) +
  scale_fill_manual(values=palfill[2]) +
  scale_color_manual(values=palline[2]) +
  scale_x_discrete(labels=Labs_NobtStage) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust=0.5),
        strip.background = element_rect(fill="#FFFFFF"),
        legend.position = "none",
        legend.justification = c(1, -1)) +
  geom_violin(position=position_dodge(width=0), alpha=0.5) +
  stat_summary(fun=mean, geom="line", aes(group=Fruit))

C2 <- list()
C2 <- lapply(seq_along(unique(Cluster_Nobt$normalized$cluster)),
             function(i) ggplot(Cluster_Nobt$normalized[Cluster_Nobt$normalized$cluster==i,], 
             aes(x=DAP, y=value, col=Species, fill=Species)) +
               labs(y="Z-score of Expression",
                    x="Stage") +
             scale_fill_manual(values=palfill[1]) +
             scale_color_manual(values=palline[1]) +
             scale_x_discrete(labels=Labs_NobtStage) +
             theme(plot.title = element_text(hjust = 0.5),
                   plot.subtitle = element_text(hjust=0.5),
                   strip.background = element_rect(fill="#FFFFFF"),
                   legend.position = "none",
                   legend.justification = c(1, -1)) +
             geom_violin(position=position_dodge(width=0), alpha=0.5) +
             stat_summary(fun=mean, geom="line", aes(group=Fruit))+
             ggtitle(paste("Cluster", i, "Expression")))


# Tobacco GO Plots --------------------------------------------------------
# Do the GO Enrichment
Tables_Nobt <- list()
for (i in levels(as.factor(Cluster_Nobt$normalized$cluster))) {
  tmpList <- as.factor(as.numeric(Cluster_Nobt$normalized$genes %in% Cluster_Nobt$normalized$genes[Cluster_Nobt$normalized$cluster==i]))
  names(tmpList) <- Cluster_Nobt$normalized$genes
  Tables_Nobt[[i]] <- GOEnrich(gene2go = "DEGAnalysis/Pfam/Nobt.gene2go.tsv",
                                GOIs=tmpList,
                               NumCategories = 10)
}
names(Tables_Nobt) <- Labs_NobtCluster[names(Tables_Nobt)]
# Create a list of the GO Plots
GO_Nobt <- list()
GO_Nobt <- lapply(seq_along(Tables_Nobt), 
                  function(i) GOPlot(Tables_Nobt[[i]], Title = paste("Cluster", i, "GO Enrichment"))+
                    theme(legend.position = "none"))


# Tobacco Gene Plots ------------------------------------------------------
# List important genes to plot
impt_genes <- c("Solyc02g077920.4.1"="CNR",
                "Solyc10g006880.3.1"="NOR",
                "Solyc05g012020.4.1"="MADS-RIN",
                "Solyc05g056620.2.1"="MC",
                "Solyc11g010570.2.1"="J",
                "Solyc12g038510.2.1"="J2",
                "Solyc03g114840.3.1"="EJ2/MADS1",
                "Solyc01g093965.2.1"="TM3",
                "Solyc01g092950.3.1"="STM3",
                "Solyc05g015750.3.1"="TM5",
                "Solyc02g089200.4.1"="TM29",
                "Solyc03g006830.3.1"="FYFL",
                "Solyc06g069430.3.1"="FUL1",
                "Solyc03g114830.3.1"="FUL2",
                "Solyc02g065730.2.1"="MBP10",
                "Solyc02g089210.4.1"="MBP20",
                "Solyc09g075440.4.1"="NR",
                "Solyc07g055920.4.1"="TAGL1",
                "Solyc01g095080.3.1"="ACS2",
                "Solyc05g050010.3.1"="ACS4",
                "Solyc07g049530.3.1"="ACO1",
                "Solyc12g005940.2.1"="ACO2",
                "Solyc07g049550.3.1"="ACO3",
                "Solyc02g081190.4.1"="ACO4",
                "Solyc07g026650.3.1"="ACO5",
                "Solyc02g036350.3.1"="ACO6",
                "Solyc06g060070.3.1"="ACO7",
                "Solyc02g071730.4.1"="TAG1",
                "Solyc11g028020.3.1"="AGL11",
                "Solyc06g064840.4.1"="MBP3",
                "Solyc03g031860.3.1"="PSY1",
                "Solyc03g111690.4.1"="PL1",
                "Solyc09g010210.3.1"="Cel2",
                "Solyc10g080210.2.1"="PGA2A",
                "Solyc06g051800.3.1"="EXP1",
                "Solyc03g098240.3.1"="GAD1",
                "Solyc11g011920.2.1"="GAD2",
                "Solyc01g005000.3.1"="GAD3",
                "Solyc09g091510.3.1"="CHS-1",
                "Solyc05g053550.3.1"="CHS-2",
                "Solyc07g052960.3.1"="GRAS38",
                "Solyc10g005060.4.1"="CTOMT1",
                "Solyc01g006540.4.1"="TOMLOXC",
                "Solyc06g073720.3.1"="EIL1",
                "Solyc01g009170.4.1"="EIL2",
                "Solyc06g073730.2.1"="EIL4",
                "Solyc07g064300.3.1"="Solyc07g064300.3.1",
                "Solyc03g117740.3.1"="Solyc03g117740.3.1",
                "Solyc05g005150.1.1"="Solyc05g005150.1.1",
                "Solyc04g072038.1.1"="Solyc04g072038.1.1",
                "Solyc06g065310.3.1"="Solyc06g065310.3.1",
                "Solyc08g066860.2.1"="Solyc08g066860.2.1",
                "Solyc07g005840.2.1"="SlCel3")
# Get Nobt orthologs
orthogroups <- read.table("Orthofinder/OrthoFinder/Results_Oct29/Orthogroups/Orthogroups.tsv",
                          sep="\t",
                          stringsAsFactors = F,
                          header=T)[,c(2:3)]
orthogroups <- orthogroups %>% filter_all(all_vars(!grepl(',',.))) #Remove multiples
Gene_Nobt <- as.data.frame(impt_genes)
tmp <- merge(Gene_Nobt, orthogroups, by.x=0, by.y="Solanum")
Gene_Nobt <- as.character(tmp$impt_genes)
names(Gene_Nobt) <- tmp$Nicotiana
Gene_Nobt <- c(Gene_Nobt, c("NIOBTv3_g28929-D2.t1"="NoFUL1","NIOBTv3_g39464.t1"="NoFUL2","NIOBTv3_g07845.t1"="NoMBP10","NIOBT_gMBP20.t1"="NoMBP20"))

# Make a new dataframe with just those genes
Subset_Nobt <- subset(Cluster_Nobt$normalized, genes %in% names(Gene_Nobt))
Subset_Nobt$Abbr <- Gene_Nobt[Subset_Nobt$genes]

# Plot the overaged expression of these genes
G1 <- list()
G1 <- lapply(seq_along(unique(Subset_Nobt$Abbr)),
             function(i) ggplot(Subset_Nobt[Subset_Nobt$Abbr==unique(Subset_Nobt$Abbr)[i],],
            aes(x=DAP, y=value, col=Species, fill=Species, group=Species)) +
  labs(y="Z-score of Expression",
       x="Stage") +
  scale_fill_manual(values=palfill[2]) +
  scale_color_manual(values=palline[2]) +
  scale_x_discrete(labels=Labs_NobtStage) +
  theme(plot.title = element_text(hjust = 0.5, face="italic"),
        plot.subtitle = element_text(hjust=0.5),
        strip.background = element_rect(fill="#FFFFFF"),
        strip.text.x = element_text(face="italic"),
        legend.position = "none") +
  geom_line(position=position_dodge(width=0)) +
    ggtitle(unique(Subset_Nobt$Abbr)[i]))

# Tobacco Figures ----------------------------------------------------------

((C2[[1]] | GO_Nobt[[1]]) / (C2[[2]] | GO_Nobt[[2]]) / (C2[[3]] | GO_Nobt[[3]]) /
 (C2[[4]] | GO_Nobt[[4]]) / (C2[[5]] | GO_Nobt[[5]]) / (C2[[6]] | GO_Nobt[[6]])) + 
  plot_annotation(tag_levels = "A")
ggsave("Figures/Tobacco_Clusters.pdf", height=20, width=15)


# Solanum Clusters ------------------------------------------------------
# For common patterns
C3 <- list()
C3 <- lapply(seq_along(unique(Cluster_Solanum_Noise$normalized$cluster)),
             function(i) ggplot(Cluster_Solanum_Noise$normalized[Cluster_Solanum_Noise$normalized$cluster==i,], 
                                aes(x=DAP, y=value, col=Species, fill=Species)) +
               labs(y="Z-score of Expression",
                    x="Stage") +
               scale_fill_manual(values=palfill[4]) +
               scale_color_manual(values=palline[4]) +
               scale_x_discrete(labels=Labs_SolanumStage) +
               theme(plot.title = element_text(hjust = 0.5),
                     plot.subtitle = element_text(hjust=0.5),
                     strip.background = element_rect(fill="#FFFFFF"),
                     legend.justification = c(1, -1)) +
               geom_violin(position=position_dodge(width=0), alpha=0.5) +
               stat_summary(fun=mean, geom="line", aes(group=Fruit))+
               ggtitle(paste("Cluster", i, "Expression")))
# For divergent patterns
C4 <- list()
C4 <- lapply(seq_along(unique(Cluster_Solanum$normalized$cluster)),
             function(i) ggplot(Cluster_Solanum$normalized[Cluster_Solanum$normalized$cluster==i,], 
                                aes(x=DAP, y=value, col=Species, fill=Species)) +
               labs(y="Z-score of Expression",
                    x="Stage") +
               scale_fill_manual(values=palfill[c(5,6)]) +
               scale_color_manual(values=palline[c(5,6)]) +
               scale_x_discrete(labels=Labs_SolanumStage) +
               theme(plot.title = element_text(hjust = 0.5),
                     plot.subtitle = element_text(hjust=0.5),
                     strip.background = element_rect(fill="#FFFFFF")) +
               geom_violin(position=position_dodge(width=0), alpha=0.5) +
               stat_summary(fun=mean, geom="line", aes(group=Species))+
               ggtitle(paste("Cluster", i, "Expression")))

# Solanum GO Plots --------------------------------------------------------
# For Common patterns
Tables_SolanumNoise <- list()
for (i in levels(as.factor(Cluster_Solanum_Noise$normalized$cluster))) {
  tmpList <- as.factor(as.numeric(Cluster_Solanum_Noise$normalized$genes %in% Cluster_Solanum_Noise$normalized$genes[Cluster_Solanum_Noise$normalized$cluster==i]))
  names(tmpList) <- Cluster_Solanum_Noise$normalized$genes
  Tables_SolanumNoise[[i]] <- GOEnrich(gene2go = "DEGAnalysis/Pfam/Slyc.gene2go.tsv",
                                  GOIs=tmpList,
                                  NumCategories = 10)
}
names(Tables_SolanumNoise) <- Labs_SolanumNoiseCluster[names(Tables_SolanumNoise)]
# Create a list of the GO Plots
GO_SolanumNoise <- list()
GO_SolanumNoise <- lapply(seq_along(Tables_SolanumNoise), 
                     function(i) GOPlot(Tables_SolanumNoise[[i]], Title = paste("Cluster", i, "GO Enrichment"))+
                       theme(legend.position = "none"))

# For divergent patterns
Tables_Solanum <- list()
for (i in levels(as.factor(Cluster_Solanum$normalized$cluster))) {
  tmpList <- as.factor(as.numeric(Cluster_Solanum$normalized$genes %in% Cluster_Solanum$normalized$genes[Cluster_Solanum$normalized$cluster==i]))
  names(tmpList) <- Cluster_Solanum$normalized$genes
  Tables_Solanum[[i]] <- GOEnrich(gene2go = "DEGAnalysis/Pfam/Slyc.gene2go.tsv",
                               GOIs=tmpList,
                               NumCategories = 10)
}
names(Tables_Solanum) <- Labs_SolanumCluster[names(Tables_Solanum)]
# Create a list of the GO Plots
GO_Solanum <- list()
GO_Solanum <- lapply(seq_along(Tables_Solanum), 
                  function(i) GOPlot(Tables_Solanum[[i]], Title = paste("Cluster", i, "GO Enrichment"))+
                    theme(legend.position = "none"))

# Solanum Gene Plots -----------------------------------------------------------

# Make a new dataframe with just those genes
Subset_Solanum <- subset(Cluster_Solanum$normalized, genes %in% names(impt_genes))
Subset_Solanum$Abbr <- impt_genes[Subset_Solanum$genes]
Subset_SolanumNoise <- subset(Cluster_Solanum_Noise$normalized, genes %in% names(impt_genes))
Subset_SolanumNoise$Abbr <- impt_genes[Subset_SolanumNoise$genes]

# Plot the overaged expression of these genes
G2 <- list()
G2 <- lapply(seq_along(unique(Subset_Solanum$Abbr)),
             function(i) ggplot(Subset_Solanum[Subset_Solanum$Abbr==unique(Subset_Solanum$Abbr)[i],], 
             aes(x=DAP, y=value, col=Species, fill=Species, group=Species, lty=Species)) +
  labs(y="Z-score of Expression",
       x="Stage") +
  scale_color_manual(values=palline[c(5,6)]) +
  scale_x_discrete(labels=Labs_SolanumStage) +
  theme(plot.title = element_text(hjust = 0.5,face="italic"),
        plot.subtitle = element_text(hjust=0.5),
        strip.background = element_rect(fill="#FFFFFF"),
        strip.text.x = element_text(face="italic")) +
  geom_line(position=position_dodge(width=0)) +
    ggtitle(unique(Subset_Solanum$Abbr)[i]))

G3 <- list()
G3 <- lapply(seq_along(unique(Subset_SolanumNoise$Abbr)),
             function(i) ggplot(Subset_SolanumNoise[Subset_SolanumNoise$Abbr==unique(Subset_SolanumNoise$Abbr)[i],], 
             aes(x=DAP, y=value, col=Species, fill=Species, group=Species)) +
  labs(y="Z-score of Expression",
       x="Stage") +
  scale_color_manual(values=palline[4]) +
  scale_x_discrete(labels=Labs_SolanumStage) +
  theme(plot.title = element_text(hjust = 0.5,face="italic"),
        plot.subtitle = element_text(hjust=0.5),
        strip.background = element_rect(fill="#FFFFFF"),
        strip.text.x = element_text(face="italic"),
        legend.position = "none") +
  geom_line(position=position_dodge(width=0))+
    ggtitle(unique(Subset_SolanumNoise$Abbr)[i]))

# Gene Figure -------------------------------------------------------------

(G1[[1]] | G2[[1]] | G3[[1]]) +
   plot_annotation(tag_levels = "A")


# Five-Species Venn Diagram -----------------------------------------------


# Five-Species GO Plots ---------------------------------------------------


# Five-Species PCA --------------------------------------------------------


# Five-Species Clusters ---------------------------------------------------


