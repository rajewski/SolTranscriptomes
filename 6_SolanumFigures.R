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
library("lemon")

#Need the plotting functions from the other figure script 6_Figures
Cluster_Solanum <- readRDS("DEGAnalysis/RNA-seq/Cluster_Solanum_3DF.rds")
Cluster_Solanum_Noise <-readRDS("DEGAnalysis/RNA-seq/Cluster_Solanum_3DF_Noise.rds")
Solanum <- Cluster_Solanum$normalized
palsol <- wes_palette("Zissou1", 2, type="continuous")

# Plot all clusters for the pimp vs AC model
SolanumAll_plot <- ggplot(Solanum, aes(x=DAP, y=value, col=Species, fill=Species)) +
  labs(y="Z-score") +
  scale_fill_manual(values=palsol) +
  scale_color_manual(values=palsol) +
  facet_rep_wrap(~cluster,
             labeller=label_both,
             nrow = 3) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust=0.5),
        strip.background = element_rect(fill="#FFFFFF")) +
  geom_violin(position=position_dodge(width=0), alpha=0.5) +
  stat_summary(fun=mean, geom="line", aes(group=Species))
ggsave("DEGAnalysis/RNA-seq/Plots/ClusterProfiles_Solanum_3DF.pdf", height=8, width=20)

SolGO <- GOPlot(GOEnrich(gene2go = "DEGAnalysis/Pfam/Slyc.gene2go.tsv",
                        GOIs="DEGAnalysis/RNA-seq/Lists/Solanum_3DF_AllGenes.txt"),
               Title = "Overall GO Enrichment") +
  theme(legend.position = "none")

# Plot all clusters for the Tom v Noise model
SolanumNoise_plot <- ggplot(Cluster_Solanum_Noise$normalized, 
                            aes(x=DAP, y=value, col=Species, fill=Species)) +
  labs(y="Z-score") +
  scale_fill_manual(values=palsol[2]) +
  scale_color_manual(values=palsol[2]) +
  facet_rep_wrap(~cluster,
                 nrow = 4) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust=0.5),
        strip.background = element_rect(fill="#FFFFFF"),
        legend.position = "none") +
  geom_violin(position=position_dodge(width=0), alpha=0.5) +
  stat_summary(fun=mean, geom="line", aes(group=Species))

ggsave("DEGAnalysis/RNA-seq/Plots/ClusterProfiles_Solanum_3DF_Noise.pdf", height=8, width=20)

# Individual Genes --------------------------------------------------------
# named vector of important genes
impt_genes <- c("CNR","NOR","RIN","MC","J","FUL1","FUL2","MBP10","MBP20","NR","TAGL1", "ACS2", "ACS4", "ACO1", "ACO2", "ACO3", "ACO4", "ACO5", "ACO6", "ACO7", "TAG1", "AGL11", "MBP3", "PSY1", "PL1", "Cel2", "PGA2A", "EXP1", "GAD1", "GAD2", "GAD3", "CHS-1", "CHS-2", "GRAS38", "CTOMT1", "TOMLOXC" )
names(impt_genes) <- c("Solyc02g077920.4.1","Solyc10g006880.3.1","Solyc05g012020.4.1","Solyc05g056620.2.1","Solyc11g010570.2.1","Solyc06g069430.3.1","Solyc03g114830.3.1","Solyc02g065730.2.1","Solyc02g089210.4.1","Solyc09g075440.4.1","Solyc07g055920.4.1", "Solyc01g095080.3.1", "Solyc05g050010.3.1", "Solyc07g049530.3.1", "Solyc12g005940.2.1", "Solyc07g049550.3.1", "Solyc02g081190.4.1", "Solyc07g026650.3.1", "Solyc02g036350.3.1", "Solyc06g060070.3.1", "Solyc02g071730.4.1", "Solyc11g028020.3.1", "Solyc06g064840.4.1" ,"Solyc03g031860.3.1", "Solyc03g111690.4.1", "Solyc09g010210.3.1", "Solyc10g080210.2.1", "Solyc06g051800.3.1", "Solyc03g098240.3.1", "Solyc11g011920.2.1", "Solyc01g005000.3.1", "Solyc09g091510.3.1", "Solyc05g053550.3.1", "Solyc07g052960.3.1", "Solyc10g005060.4.1", "Solyc01g006540.4.1")

# make a new dataframe to get the ordering of gene names right
Subset_Solanum_Noise <- subset(Cluster_Solanum_Noise$normalized, genes %in% names(impt_genes))
Subset_Solanum_Noise$Abbr <- impt_genes[Subset_Solanum_Noise$genes]

Subset_Solanum <- subset(Cluster_Solanum$normalized, genes %in% names(impt_genes))
Subset_Solanum$Abbr <- impt_genes[Subset_Solanum$genes]

# Plot of important genes with differential expression between the two species
Indiv_Plot_1 <- ggplot(Subset_Solanum,
                     aes(x=DAP, y=value, col=Species, fill=Species, group=Species)) +
  labs(y="Z-score") +
  scale_fill_manual(values=palsol) +
  scale_color_manual(values=palsol) +
  facet_wrap(~Abbr) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust=0.5),
        strip.background = element_rect(fill="#FFFFFF"),
        strip.text.x = element_text(face="italic")) +
  geom_line(position=position_dodge(width=0))

# Plot of important genes with differentail expression in common among the solanum species
Indiv_Plot_2 <- ggplot(Subset_Solanum_Noise,
                       aes(x=DAP, y=value, col=Species, fill=Species, group=Species)) +
  labs(y="Z-score") +
  scale_fill_manual(values=palsol[2]) +
  scale_color_manual(values=palsol[2]) +
  facet_rep_wrap(~Abbr,
             nrow=4,
             ) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust=0.5),
        legend.position = "bottomright",
        strip.background = element_rect(fill="#FFFFFF"),
        strip.text.x = element_text(face="italic")) +
  geom_line(position=position_dodge(width=0))
