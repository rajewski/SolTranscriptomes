source("X_Functions.R")
library("DESeq2")
library("ggplot2")
library("cowplot")
theme_set(theme_cowplot())
library("lemon")
library("VennDiagram")
library("dplyr")
library("reshape")
library("patchwork")

# Read Data In ------------------------------------------------------------
Cluster_Nobt <- readRDS("DEGAnalysis/RNA-seq/Cluster_Nobt.rds")
Labs_NobtCluster <- ClusterLabs(Cluster_Nobt)

DDS_Nobt <- readRDS("DEGAnalysis/RNA-seq/DDS_NobtSE.rds")
Labs_NobtStage <- c("0"="1", "3"="2", "6"="3", "11"="Br")

Cluster_Solanum <- readRDS("DEGAnalysis/RNA-seq/Cluster_Solanum_3DF.rds")
Labs_SolanumCluster <- ClusterLabs(Cluster_Solanum)
Cluster_Solanum_Noise <-readRDS("DEGAnalysis/RNA-seq/Cluster_Solanum_3DF_Noise.rds")
Labs_SolanumNoiseCluster <- ClusterLabs(Cluster_Solanum_Noise)

DDS_SolanumNoise <- readRDS("DEGAnalysis/RNA-seq/DDS_Solanum_3DF_Noise.rds")
DDS_Solanum <- readRDS("DEGAnalysis/RNA-seq/DDS_Solanum_3DF.rds")
Labs_SolanumStage <- c("1"="1", "3"="2", "15"="3", "35"="Br", "45"="RR")

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
                "Solyc07g052960.3.1"="SlFSR",
                "Solyc10g005060.4.1"="CTOMT1",
                "Solyc01g006540.4.1"="TOMLOXC",
                "Solyc06g073720.3.1"="EIL1",
                "Solyc01g009170.4.1"="EIL2",
                "Solyc06g073730.2.1"="EIL4",
                "Solyc07g064300.3.1"="Solyc07g064300.3.1",
                "Solyc03g117740.3.1"="Solyc03g117740.3.1",
                "Solyc05g005150.1.1"="SlKFB",
                "Solyc04g072038.1.1"="Solyc04g072038.1.1",
                "Solyc06g065310.3.1"="Solyc06g065310.3.1",
                "Solyc08g066860.2.1"="SlFKD",
                "Solyc07g005840.2.1"="Cel3")

Limit_genes <- c(-2.5,2.5)

Cluster_Noise <- readRDS("DEGAnalysis/RNA-seq/Cluster_FiveOrtho_Noise.rds")
Cluster_Fruit <- readRDS("DEGAnalysis/RNA-seq/Cluster_FiveOrtho_Fruit.rds")
DDS_Noise <- readRDS("DEGAnalysis/RNA-seq/DDS_FiveOrtho_Noise.rds")
DDS_Fruit <- readRDS("DEGAnalysis/RNA-seq/DDS_FiveOrtho_Fruit.rds")

# Tobacco Clusters -----------------------------------------------------------
C1 <- list()
C1 <- lapply(seq_along(unique(Cluster_Nobt$normalized$cluster)),
             function(i) ggplot(Cluster_Nobt$normalized[Cluster_Nobt$normalized$cluster==i,], 
                                aes(x=DAP, y=value, col=Species, fill=Species)) +
               labs(y="Z-score of Expression",
                    x="Stage") +
               scale_fill_manual(values=palw2[1]) +
               scale_color_manual(values=palw2[1]) +
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
                  function(i) GOPlot(Tables_Nobt[[i]], 
                                     Title = paste("Cluster", i, "GO Enrichment"),
                                     colorHex = palw2[1]) +
                    theme(legend.position = "none"))


# Tobacco Gene Plots ------------------------------------------------------
#Get Nobt orthologs
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
rm(tmp)

#Rename and hardcode genes
Gene_Nobt <- c('NIOBTv3_g02352.t1'='NoACO6', 
               'NIOBTv3_g22632-D2.t1'='NoAG', 
               'NIOBTv3_g27953.t1'='NoCNR', 
               'NIOBTv3_g13660.t1'='NoACO4', 
               'NIOBTv3_g14235.t1'='NoSEP1', 
               'NIOBTv3_g10096.t1'='NoFYFL', 
               'NIOBTv3_g17569.t1'='NoPSY1', 
               'NIOBTv3_g11084.t1'='NoGAD1', 
               'NIOBTv3_g22270.t1'='NIOBTv3_g22270.t1', 
               'NIOBTv3_g10008.t1'='NIOBTv3_g10008.t1', 
               'NIOBTv3_g35607.t1'='NoKFB', 
               'NIOBTv3_g18077.t1'='NoMC', 
               'NIOBTv3_g12238.t1'='NIOBTv3_g12238.t1', 
               'NIOBTv3_g12440.t1'='NoCel3', 
               'NIOBTv3_g38689.t1'='NoACO5', 
               'NIOBTv3_g37943.t1'='NoGRAS38', 
               'NIOBTv3_g13969.t1'='NoSHP', 
               'NIOBTv3_g11662.t1'='NIOBTv3_g11662.t1', 
               'NIOBTv3_g12218.t1'='NoFKD', 
               'NIOBTv3_g19880.t1'='NoCel2', 
               'NIOBTv3_g10291.t1'='NoNR', 
               'NIOBTv3_g14436.t1'='NoAGL11', 
               'NIOBTv3_g28929-D2.t1'='NoFUL1', 
               'NIOBTv3_g39464.t1'='NoFUL2', 
               'NIOBTv3_g07845.t1'='NoMBP10', 
               'NIOBT_gMBP20.t1'='NoMBP20')

#Test if I should plot a gene at all
Test_Nobt <- as.data.frame(results(DDS_Nobt))
Test_Nobt <- Test_Nobt[row.names(Test_Nobt) %in% names(Gene_Nobt),]
Test_Nobt$Abbr <- Gene_Nobt[row.names(Test_Nobt)]
Test_Nobt <- Test_Nobt[order(Test_Nobt$Abbr),]
Test_Nobt$Plot <- TRUE
Test_Nobt$Plot[Test_Nobt$padj>0.01 | is.na(Test_Nobt$padj)] <- FALSE
Test_Nobt <- Test_Nobt[,-c(1:5)]


# Make a new dataframe with just those genes
# Subset to just the important genes
Subset_Nobt <- as.data.frame(assay(subset(DDS_Nobt, rownames(DDS_Nobt) %in% names(Gene_Nobt))))
Subset_Nobt$Gene <- rownames(Subset_Nobt)
Subset_Nobt$Abbr <- Gene_Nobt[Subset_Nobt$Gene]
Subset_Nobt <- Subset_Nobt[order(Subset_Nobt$Abbr),]
Subset_Nobt <- melt(Subset_Nobt, id.vars = c("Gene", "Abbr"))
Subset_Nobt$DAP <- as.factor(rep(colData(DDS_Nobt)$DAP,
                                 each=length(unique(Subset_Nobt$Gene))))

G1 <- list()
G1 <- lapply(seq_along(unique(Subset_Nobt$Abbr)),
             function(i) ggplot(Subset_Nobt[Subset_Nobt$Abbr==unique(Subset_Nobt$Abbr)[i],], 
                                aes(x=DAP, y=value, col=Abbr, fill=Abbr)) +
               labs(y="Normalized Counts",
                    x="Stage") +
               scale_color_manual(values=palw2[1]) +
               scale_fill_manual(values=palw2[1]) +
               scale_x_discrete(labels=Labs_NobtStage) +
               theme(plot.title = element_text(hjust = 0.5,
                                               face = ifelse(Test_Nobt[i,"Plot"], 'bold.italic', 'italic')),
                     plot.subtitle = element_text(hjust=0.5),
                     strip.background = element_rect(fill="#FFFFFF"),
                     strip.text.x = element_text(face="italic"),
                     legend.position = "none") +
               geom_violin(position=position_dodge(width=0), alpha=0.5) +
               stat_summary(fun=mean, geom="line", group="Abbr")+
               ggtitle(unique(Subset_Nobt$Abbr)[i]))

# Tobacco Figures ----------------------------------------------------------
(G1[[5]] + G1[[6]] + G1[[7]] + G1[[8]] + G1[[9]] + G1[[10]] + 
   G1[[11]] + G1[[12]] + G1[[13]] + G1[[14]] + G1[[15]] + 
  G1[[17]] + G1[[18]] + G1[[19]] +
   G1[[21]] + G1[[22]] + G1[[23]] + G1[[24]] + G1[[25]] + 
   G1[[26]]) +
  plot_annotation(tag_levels = "A")
ggsave2("Figures/Tobacco_Genes.pdf", height=15, width=15)

((C1[[1]] | GO_Nobt[[1]]) / (C1[[2]] | GO_Nobt[[2]]) / (C1[[3]] | GO_Nobt[[3]]) /
    (C1[[4]] | GO_Nobt[[4]]) / (C1[[5]] | GO_Nobt[[5]]) / (C1[[6]] | GO_Nobt[[6]])) + 
  plot_annotation(tag_levels = "A")
ggsave("Figures/Tobacco_Clusters.pdf", height=20, width=15)

# Solanum Clusters ------------------------------------------------------
# For common patterns
C2 <- list()
C2 <- lapply(seq_along(unique(Cluster_Solanum_Noise$normalized$cluster)),
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
C3 <- list()
C3 <- lapply(seq_along(unique(Cluster_Solanum$normalized$cluster)),
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



# Solanum Gene Plots DDS --------------------------------------------------
# Determine if a plot should be combined by species or now
Test_Solanum <- merge(as.data.frame(results(DDS_Solanum)),
                      as.data.frame(results(DDS_SolanumNoise)),
                      by=0)
Test_Solanum <- Test_Solanum[Test_Solanum$Row.names %in% names(impt_genes),c(1,7,13)]
Test_Solanum$Abbr <- impt_genes[Test_Solanum$Row.names]
Test_Solanum <- Test_Solanum[order(Test_Solanum$Abbr),]
Test_Solanum$padj.x[Test_Solanum$padj.x>0.01] <- 1
Test_Solanum$padj.y[Test_Solanum$padj.y>0.01] <- 1
Test_Solanum$Choose[Test_Solanum$padj.x < Test_Solanum$padj.y] <- "Separate"
Test_Solanum$Choose[Test_Solanum$padj.y < Test_Solanum$padj.x] <- "Together"
Test_Solanum$Choose[is.na(Test_Solanum$Choose)] <- "Neither"

# Subset to just the important genes
Subset_SolanumNoise <- as.data.frame(assay(subset(DDS_SolanumNoise, rownames(DDS_SolanumNoise) %in% names(impt_genes))))
Subset_SolanumNoise$Gene <- rownames(Subset_SolanumNoise)
Subset_SolanumNoise$Abbr <- impt_genes[Subset_SolanumNoise$Gene]
Subset_SolanumNoise <- Subset_SolanumNoise[order(Subset_SolanumNoise$Abbr),]
Subset_SolanumNoise <- melt(Subset_SolanumNoise, id.vars = c("Gene", "Abbr"))
Subset_SolanumNoise$Species <- rep(colData(DDS_SolanumNoise)$Species,
                                   each=length(unique(Subset_SolanumNoise$Gene)))
Subset_SolanumNoise$DAP <- as.factor(rep(colData(DDS_SolanumNoise)$DAP,
                               each=length(unique(Subset_SolanumNoise$Gene))))

G2 <- list()
G2 <- lapply(seq_along(unique(Subset_SolanumNoise$Abbr)),
             function(i) ggplot(Subset_SolanumNoise[Subset_SolanumNoise$Abbr==unique(Subset_SolanumNoise$Abbr)[i],], 
                                aes(x=DAP, y=value, col=Abbr, fill=Abbr)) +
               labs(y="Normalized Counts",
                    x="Stage") +
               #ylim(Limit_genes) +
               scale_color_manual(values=palw2[3]) +
               scale_fill_manual(values=palw2[3]) +
               scale_x_discrete(labels=Labs_SolanumStage) +
               theme(plot.title = element_text(hjust = 0.5,
                                  face = ifelse(Test_Solanum[i,"Choose"]=="Together", 'bold.italic', 'italic')),
                     plot.subtitle = element_text(hjust=0.5),
                     strip.background = element_rect(fill="#FFFFFF"),
                     strip.text.x = element_text(face="italic"),
                     legend.position = "none") +
               geom_violin(position=position_dodge(width=0), alpha=0.5) +
               stat_summary(fun=mean, geom="line", group="Abbr")+
               ggtitle(unique(Subset_SolanumNoise$Abbr)[i]))

# Subset to just the important genes
Subset_Solanum <- as.data.frame(assay(subset(DDS_Solanum, rownames(DDS_Solanum) %in% names(impt_genes))))
Subset_Solanum$Gene <- rownames(Subset_Solanum)
Subset_Solanum$Abbr <- impt_genes[Subset_Solanum$Gene]
Subset_Solanum <- Subset_Solanum[order(Subset_Solanum$Abbr),]
Subset_Solanum <- melt(Subset_Solanum, id.vars = c("Gene", "Abbr"))
Subset_Solanum$Species <- rep(colData(DDS_Solanum)$Species,
                              each=length(unique(Subset_Solanum$Gene)))
Subset_Solanum$DAP <- as.factor(rep(colData(DDS_Solanum)$DAP,
                                    each=length(unique(Subset_Solanum$Gene))))

G3 <- list()
G3 <- lapply(seq_along(unique(Subset_Solanum$Abbr)),
             function(i) ggplot(Subset_Solanum[Subset_Solanum$Abbr==unique(Subset_Solanum$Abbr)[i],], 
                                aes(x=DAP, y=value, col=Species, fill=Species)) +
               labs(y="Normalized Counts",
                    x="Stage") +
               scale_color_manual(values=palw2[c(2,3)],
                                  labels=c("Wild", "Cultivated or Both")) +
               scale_fill_manual(values=palw2[c(2,3)],
                                 labels=c("Wild", "Cultivated or Both")) +
               scale_x_discrete(labels=Labs_SolanumStage) +
               theme(plot.title = element_text(hjust = 0.5,
                                               face = ifelse(Test_Solanum[i,"Choose"]=="Separate", 'bold.italic', 'italic')),
                     plot.subtitle = element_text(hjust=0.5),
                     strip.background = element_rect(fill="#FFFFFF"),
                     strip.text.x = element_text(face="italic")) +
               geom_violin(position=position_dodge(width=0), alpha=0.5) +
               stat_summary(fun=mean, geom="line", aes(group=Species)) +
               #stat_summary(fun.data=mean_se, geom="errorbar", aes(group=Species)) +
               ggtitle(unique(Subset_Solanum$Abbr)[i]))


# Gene Figure -------------------------------------------------------------
#Ethylene Biosynth and perception
(G2[[1]] + G2[[2]] + G2[[3]] + G2[[4]] + G2[[5]] + G3[[6]] + G2[[7]] + G2[[8]] + 
   G2[[9]] + G2[[17]] +  G2[[18]] + G2[[19]] + G2[[36]] + guide_area()) +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides="collect")
ggsave2("Figures/Tomato_Ethylene_Genes.pdf",
        height=12,
        width=12)

#Flavor and Pigment
(G2[[13]] + G2[[14]] + G2[[15]] + G2[[16]] + G2[[25]] + G2[[26]] + G2[[27]] +
   G2[[39]] +G2[[53]])  +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides="collect") 
ggsave2("Figures/Tomato_Pigment_Genes.pdf",
        height=10,
        width=10)

#Developmental TFs
(G2[[10]] + G2[[20]] + G2[[22]] + G2[[23]] + G2[[24]] + G2[[28]] +
    G2[[29]] + G2[[30]] + G2[[31]] + G2[[32]] + G2[[33]] + G3[[34]] +
    G2[[35]] + G2[[36]] + G2[[47]] + G3[[48]] + G3[[49]] + G2[[50]] + 
    G2[[51]] + G2[[52]] + guide_area())  +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides="collect",
              nrow=3) 
ggsave2("Figures/Tomato_TF_Genes.pdf",
        height=10,
        width=21)

#Cell Wall
(G2[[11]] + G2[[12]] + G2[[21]] + G2[[37]] + G2[[38]] + G2[[39]] + G2[[40]])
ggsave2("Figures/Tomato_Cell_Genes.pdf",
        height=10,
        width=10)

#Misc
(G2[[41]] + G2[[42]] + G3[[43]] + G2[[44]] + G2[[45]] + G2[[46]]) +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides="collect")
ggsave2("Figures/Tomato_Misc_Genes.pdf", height=30, width=30)

 
# Five-Species Venn Diagram -----------------------------------------------
VennOrthos <- read.table("Orthofinder/OrthoFinder/Results_Aug31/Orthogroups/Orthogroups.tsv",
                         sep="\t",
                         stringsAsFactors = F,
                         header=T)
#All Orthogroups
# Remove empties for each species
Set1 <- VennOrthos[,c(1:2)] %>% filter_all(all_vars(!grepl("^$",.)))
Set2 <- VennOrthos[,c(1,3)] %>% filter_all(all_vars(!grepl("^$",.)))
Set3 <- VennOrthos[,c(1,4)] %>% filter_all(all_vars(!grepl("^$",.)))
Set4 <- VennOrthos[,c(1,5)] %>% filter_all(all_vars(!grepl("^$",.)))
venn.diagram(
  x = list(Set1[,1], Set2[,1], Set3[,1], Set4[,1]),
  category.names = c("Arabidopsis","Cucumis", "Nicotiana", "Solanum"),
  filename = "Figures/AllOrthogroups.png",
  imagetype = "png",
  main="Orthogenes",
  main.fontfamily = "sans",
  main.fontface = "bold",
  cat.fontface="italic",
  cat.fontfamily="sans",
  fontfamily="sans",
  scaled=T,
  euler.d=T,
  output=F,
  lwd=2,
  lty=1)
# All Species, single copy
SingleSet1 <- Set1 %>% filter_all(all_vars(!grepl(',',.)))
SingleSet2 <- Set2 %>% filter_all(all_vars(!grepl(',',.)))
SingleSet3 <- Set3 %>% filter_all(all_vars(!grepl(',',.)))
SingleSet4 <- Set4 %>% filter_all(all_vars(!grepl(',',.)))
venn.diagram(
  x = list(SingleSet1[,1], SingleSet2[,1], SingleSet3[,1], SingleSet4[,1]),
  category.names = c("Arabidopsis","Cucumis", "Nicotiana", "Solanum"),
  filename = "Figures/SingleCopyOrthogroups.png",
  imagetype = "png",
  main="Shared Single-Copy Orthogenes",
  main.fontfamily = "sans",
  main.fontface = "bold",
  cat.fontface="italic",
  cat.fontfamily="sans",
  fontfamily="sans",
  scaled=T,
  euler.d=T,
  output=F,
  lwd=2,
  lty=1)

# Five-Species GO Plots ---------------------------------------------------
# First get all the clusters when grouped by fruit
Tables_Ortho <- list()
Tables_Ortho <- lapply(seq_along(unique(Cluster_Fruit$normalized$cluster)), 
                       function(i)  GOEnrich(gene2go = "DEGAnalysis/Pfam/Ortho.831All.gene2go.tsv",
                                             GOIs=setNames(as.factor(as.numeric(Cluster_Fruit$normalized$genes %in% Cluster_Fruit$normalized$genes[Cluster_Fruit$normalized$cluster==unique(Cluster_Fruit$normalized$cluster)[i]])),
                                                           Cluster_Fruit$normalized$genes)))
# Then add in Overalls
Tables_Ortho[["Model 1"]] <- GOEnrich(gene2go = "DEGAnalysis/Pfam/Ortho.831All.gene2go.tsv",
                                     GOIs=setNames(as.factor(as.numeric(rownames(DDS_Noise) %in% rownames(subset(results(DDS_Noise), padj<0.01)))), rownames(DDS_Noise)))

Tables_Ortho[["Model 2"]] <- GOEnrich(gene2go = "DEGAnalysis/Pfam/Ortho.831All.gene2go.tsv",
                                     GOIs=setNames(as.factor(as.numeric(rownames(DDS_Fruit) %in% rownames(subset(results(DDS_Fruit), padj<0.01)))), rownames(DDS_Fruit)))
                                     
# Save the GO Plots
GO_Ortho <- list()
GO_Ortho <- lapply(seq_along(Tables_Ortho), 
                     function(i) GOPlot(Tables_Ortho[[i]],
                                        Title = ifelse(names(Tables_Ortho[i])=='',
                                                             paste("Cluster", i, "GO Enrichment"),
                                                             paste(names(Tables_Ortho[i]), "GO Enrichment")),
                                        colorHex = "#9EBE91CC")+
                       theme(legend.position = "none"))

# Five-Species PCA --------------------------------------------------------
Subset_Noise <- DDS_Noise[rownames(DDS_Noise) %in% rownames(subset(results(DDS_Noise), padj<=0.01))]
RLD_Noise <- rlog(Subset_Noise, blind=FALSE)
RLD_Noise$Stage[RLD_Noise$Stage=="3.5"] <- "Br"
RLD_Noise$Stage <- as.factor(RLD_Noise$Stage)
RLD_Noise$Species <- factor(RLD_Fruit$Species, levels=c("Arabidopsis", "Melon", "Tobacco", "Tomato", "Pimpinellifolium"))
# Make three PCA plots
PCA1 <- lapply(list(c(1,2),c(2,3),c(1,3)),
               function(i) plotPCAmod(RLD_Noise,xPC=i[1], yPC=i[2]))

Subset_Fruit <- DDS_Fruit[rownames(DDS_Fruit) %in% rownames(subset(results(DDS_Fruit), padj<=0.01))]
RLD_Fruit <- rlog(Subset_Fruit, blind=FALSE)
RLD_Fruit$Stage[RLD_Fruit$Stage=="3.5"] <- "Br"
RLD_Fruit$Stage <- as.factor(RLD_Fruit$Stage)
RLD_Fruit$Species <- factor(RLD_Fruit$Species, levels=c("Arabidopsis", "Melon", "Tobacco", "Tomato", "Pimpinellifolium"))
PCA2 <-  lapply(list(c(1,2),c(2,3),c(1,3)),
                function(i) plotPCAmod(RLD_Fruit,
                                       xPC=i[1], yPC=i[2],
                                       ntop=dim(RLD_Fruit)[1]))

# Five-Species Clusters ---------------------------------------------------

#From 6_Dry_v_Fleshy.R

