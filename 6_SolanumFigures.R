# To make the figure for the pre-ripening transcriptome data
library("DEGreport")
library("ggplot2")
library("cowplot")
theme_set(theme_cowplot())
library("topGO")
library("tidyr")
library("dplyr")
library("Rgraphviz")
library("scales")
library("lemon")
source("X_Functions.R")

# Read in the datasets ----------------------------------------------------
Cluster_Solanum <- readRDS("DEGAnalysis/RNA-seq/Cluster_Solanum_3DF.rds")
Cluster_Solanum_Noise <-readRDS("DEGAnalysis/RNA-seq/Cluster_Solanum_3DF_Noise.rds")


# Cluster Profile Plots ---------------------------------------------------
# Plot all clusters for the pimp vs AC model
SA_Labs <- paste("Cluster", 1:length(unique(Cluster_Solanum$normalized$cluster)))
names(SA_Labs) <- sort(unique(Cluster_Solanum$normalized$cluster))

SolanumAll_plot <- ggplot(Cluster_Solanum$normalized, aes(x=DAP, y=value, col=Species, fill=Species)) +
  labs(y="Z-score of Expression") +
  scale_fill_manual(values=palw[c(1,4)]) +
  scale_color_manual(values=palw[c(1,4)]) +
  facet_rep_wrap(~cluster,
             labeller=labeller(cluster=SA_Labs),
             nrow = 3) +
  scale_x_discrete(labels=c("1" = "1", "3" = "3", "15" = "15", "35" = "Br", "45"="RR")) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust=0.5),
        strip.background = element_rect(fill="#FFFFFF")) +
  geom_violin(position=position_dodge(width=0), alpha=0.5) +
  stat_summary(fun=mean, geom="line", aes(group=Species))
SolanumAll_plot
ggsave("DEGAnalysis/RNA-seq/Plots/ClusterProfiles_Solanum_3DF.pdf", height=8, width=20)

# Plot all clusters for the Tom v Noise model
SN_Labs <- paste("Cluster", 1:length(unique(Cluster_Solanum_Noise$normalized$cluster)))
names(SN_Labs) <- sort(unique(Cluster_Solanum_Noise$normalized$cluster))
SolanumNoise_plot <- ggplot(Cluster_Solanum_Noise$normalized, 
                            aes(x=DAP, y=value, col=Species, fill=Species)) +
  labs(y="Z-score of Expression") +
  scale_fill_manual(values=palw[4]) +
  scale_color_manual(values=palw[4]) +
  facet_rep_wrap(~cluster,
                 nrow = 4,
                 labeller=labeller(cluster=SN_Labs)) +
  scale_x_discrete(labels=c("1" = "1", "3" = "3", "15" = "15", "35" = "Br", "45"="RR")) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust=0.5),
        strip.background = element_rect(fill="#FFFFFF"),
        legend.position = "none") +
  geom_violin(position=position_dodge(width=0), alpha=0.5) +
  stat_summary(fun=mean, geom="line", aes(group=Species))
SolanumNoise_plot
ggsave("DEGAnalysis/RNA-seq/Plots/ClusterProfiles_Solanum_3DF_Noise.pdf", height=8, width=20)

# Make GO Plots and Lists -----------------------------------------------------------
# Cohort of all genes in first model
GOPlot(GOEnrich(gene2go = "DEGAnalysis/Pfam/Slyc.gene2go.tsv",
                         GOIs="DEGAnalysis/RNA-seq/Lists/Solanum_3DF_AllGenes.txt"),
                Title = "Overall GO Enrichment")
ggsave2("DEGAnalysis/RNA-seq/Plots/Solanum_3DF_AllGO.pdf", height=8, width=11 )
# Cohort of all genes in second model
GOPlot(GOEnrich(gene2go = "DEGAnalysis/Pfam/Slyc.gene2go.tsv",
                         GOIs="DEGAnalysis/RNA-seq/Lists/Solanum_3DF_Noise_AllGenes.txt"),
                Title = "Overall GO Enrichment")
ggsave2("DEGAnalysis/RNA-seq/Plots/Solanum_3DF_Noise_AllGO.pdf", height=8, width=11 )

# Get table of GOs by cluster
SolTables <- list()
for (i in levels(as.factor(Cluster_Solanum$normalized$cluster))) {
  tmpList <- as.factor(as.numeric(Cluster_Solanum$normalized$genes %in% Cluster_Solanum$normalized$genes[Cluster_Solanum$normalized$cluster==i]))
  names(tmpList) <- Cluster_Solanum$normalized$genes
  SolTables[[i]] <- GOEnrich(gene2go = "DEGAnalysis/Pfam/Slyc.gene2go.tsv",GOIs=tmpList)}
names(SolTables) <- SA_Labs[names(SolTables)]
capture.output(SolTables, file="DEGAnalysis/RNA-seq/Solanum_3DF_GOTables.txt")

SolNoiseTables <- list()
for (i in levels(as.factor(Cluster_Solanum_Noise$normalized$cluster))) {
  tmpList <- as.factor(as.numeric(Cluster_Solanum_Noise$normalized$genes %in% Cluster_Solanum_Noise$normalized$genes[Cluster_Solanum_Noise$normalized$cluster==i]))
  names(tmpList) <- Cluster_Solanum_Noise$normalized$genes
  SolNoiseTables[[i]] <- GOEnrich(gene2go = "DEGAnalysis/Pfam/Slyc.gene2go.tsv",
                               GOIs=tmpList)}
names(SolNoiseTables) <- SN_Labs[names(SolNoiseTables)]
capture.output(SolNoiseTables, file="DEGAnalysis/RNA-seq/Solanum_3DF_Noise_GOTables.txt")

# Plot Individual Genes --------------------------------------------------------
# named vector of important genes
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

# make a new dataframe to get the ordering of gene names right
Subset_Solanum_Noise <- subset(Cluster_Solanum_Noise$normalized, genes %in% names(impt_genes))
Subset_Solanum_Noise$Abbr <- impt_genes[Subset_Solanum_Noise$genes]

Subset_Solanum <- subset(Cluster_Solanum$normalized, genes %in% names(impt_genes))
Subset_Solanum$Abbr <- impt_genes[Subset_Solanum$genes]

# Plot of differential expression between the two species
Indiv_Plot_1 <- ggplot(Subset_Solanum,
                     aes(x=DAP, y=value, col=Species, fill=Species, group=Species)) +
  labs(y="Z-score of Expression") +
  scale_fill_manual(values=palw[c(1,4)]) +
  scale_color_manual(values=palw[c(1,4)]) +
  facet_wrap(~Abbr) +
  scale_x_discrete(labels=c("1" = "1", "3" = "3", "15" = "15", "35" = "Br", "45"="RR")) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust=0.5),
        strip.background = element_rect(fill="#FFFFFF"),
        strip.text.x = element_text(face="italic")) +
  geom_line(position=position_dodge(width=0))

# Plot of  differential expression in common 
Indiv_Plot_2 <- ggplot(Subset_Solanum_Noise,
                       aes(x=DAP, y=value, col=Species, fill=Species, group=Species)) +
  labs(y="Z-score of Expression") +
  scale_fill_manual(values=palw[4]) +
  scale_color_manual(values=palw[4]) +
  facet_rep_wrap(~Abbr,
                 nrow=4) +
  scale_x_discrete(labels=c("1" = "1", "3" = "3", "15" = "15", "35" = "Br", "45"="RR")) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust=0.5),
        legend.position = "none",
        strip.background = element_rect(fill="#FFFFFF"),
        strip.text.x = element_text(face="italic")) +
  geom_line(position=position_dodge(width=0))
