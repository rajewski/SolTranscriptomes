source("X_Functions.R")
library("DESeq2")
library("ggplot2")
library("cowplot")
theme_set(theme_cowplot())
library("VennDiagram")
library("dplyr")
library("reshape")
library("patchwork")
library("gt")

# Read Data In ------------------------------------------------------------
Cluster_Nobt <- readRDS("DEGAnalysis/RNA-seq/Cluster_Nobt.rds")
DDS_Nobt <- readRDS("DEGAnalysis/RNA-seq/DDS_NobtSE.rds")
Labs_NobtStage <- c("0"="1", "3"="2", "6"="3", "11"="Tr")

Cluster_Solanum <- readRDS("DEGAnalysis/RNA-seq/Cluster_Solanum_3DF.rds")
Cluster_Solanum_Noise <-readRDS("DEGAnalysis/RNA-seq/Cluster_Solanum_3DF_Noise.rds")
DDS_SolanumNoise <- readRDS("DEGAnalysis/RNA-seq/DDS_Solanum_3DF_Noise.rds")
DDS_Solanum <- readRDS("DEGAnalysis/RNA-seq/DDS_Solanum_3DF.rds")
Labs_SolanumStage <- c("1"="1", "3"="2", "15"="3", "35"="Br", "45"="RR")

Cluster_Solanaceae_Noise <- readRDS("DEGAnalysis/RNA-seq/Cluster_Solanaceae_Noise.rds")
Cluster_Solanaceae <- readRDS("DEGAnalysis/RNA-seq/Cluster_Solanaceae.rds")
DDS_Solanaceae_Noise <- readRDS("DEGAnalysis/RNA-seq/DDS_Solanaceae_Noise.rds")
DDS_Solanaceae <- readRDS("DEGAnalysis/RNA-seq/DDS_Solanaceae.rds")

Cluster_Noise <- readRDS("DEGAnalysis/RNA-seq/Cluster_FiveOrtho_Noise.rds")
Cluster_Fruit <- readRDS("DEGAnalysis/RNA-seq/Cluster_FiveOrtho_Fruit.rds")
Labs_Stage <- c("1"="1","2"="2", "3"="3", "3.5"="Tr")
DDS_Noise <- readRDS("DEGAnalysis/RNA-seq/DDS_FiveOrtho_Noise.rds")
DDS_Fruit <- readRDS("DEGAnalysis/RNA-seq/DDS_FiveOrtho_Fruit.rds")

# Data Source and Stage Table ---------------------------------------------
#PATH=/bigdata/littlab/arajewski/Datura/software/phantomjs-2.1.1-linux-x86_64/bin:$PATH
Design <- read.csv("ExternalData/Design.csv", stringsAsFactors = F)[,-1]
row.names(Design) <- c("Stage 1",
                       "Stage 2",
                       "Stage 3",
                       "Transition",
                       "Stage 4",
                       "Bioproject Accession")

Design %>%
  gt(rownames_to_stub = TRUE) %>%
  tab_header(title="Description of Developmental Stages") %>%
  fmt_missing(
    columns = ends_with("Day")) %>%
  cols_merge(
    columns = starts_with("Nobtusifolia"),
    hide_columns = vars("NobtusifoliaDay"),
    pattern = "{1}</br>({2})") %>%
  cols_merge(
    columns = starts_with("Athaliana"),
    hide_columns = vars("AthalianaDay"),
    pattern = "{1}</br>({2})") %>%
  cols_merge(
    columns = starts_with("Solanum"),
    hide_columns = vars("SolanumDay"),
    pattern = "{1}</br>({2})") %>%
  cols_merge(
    columns = starts_with("Cmelo"),
    hide_columns = vars("CmeloDay"),
    pattern = "{1}</br>({2})") %>%
  cols_label(
    NobtusifoliaDesc = md("**Desert Tobacco**"),
    AthalianaDesc = md("***A.&nbsp;thaliana***"),
    SolanumDesc = md("**Tomato**"),
    CmeloDesc = md("**Melon**")) %>%
  tab_footnote(
    footnote = "Not sampled",
    locations = list(
      cells_body(columns=starts_with(c("Nobt", "Atha")),
                 rows=5),
      cells_body(columns=starts_with("Cmel"),
                 rows=1))) %>%
  tab_source_note(
    source_note = md("Pab&#243;n-Mora and Litt, 2011; Mizzotti et al, 2018; and Zhang et al, 2016")) %>%
  gtsave(filename = "Design.rtf",
         path = "Tables")

# Simplified Table for defense presentation
Stages <- data.frame(t(Design)[c(4,2,8,6,6),],
                     stringsAsFactors = F)
row.names(Stages) <- c("Arabidopsis", "Desert Tobacco", "Melon", "Tomato (Cultivated)", "Tomato (Wild)")
Stages[,1:5] <- Stages[,1:5] %>% mutate_all(~replace(., !is.na(.), "&#10003;"))
Stages[,1:5] <- Stages[,1:5] %>% mutate_all(~replace(., is.na(.), "&#8709;"))
Stages$Type <- c("Dry - Silique", "Dry - Capsule", "Fleshy - Pepo", "Fleshy - Berry", "Fleshy - Berry")

Stages %>%
  gt(
    rownames_to_stub = TRUE) %>%
  fmt_markdown(columns=c(2:6)) %>%
  cols_label(
    Stage.1 = md("**Stage&nbsp;1**"),
    Stage.2 = md("**Stage&nbsp;2**"),
    Stage.3 = md("**Stage&nbsp;3**"),
    Transition = md("**Transition**"),
    Stage.4 = md("**Stage&nbsp;4**"),
    Type = md("**Fruit Type**"),
    Bioproject.Accession = md("**Source**")) %>%
  cols_move_to_end(vars(Bioproject.Accession)) %>%
  gtsave(path = "Defense",
         filename = "Stages.rtf")
rm(Stages)

# Important Genes for Plotting --------------------------------------------
#unimportant genes
#Solyc06g068440.4.1 - CCR1 - Not DE
#Solyc03g117600.3.1 - HCT - conserved 
#Solyc09g007920.4.1 - PAL - conserved
#Solyc12g042460.2.1 - 4CL - not DE
#Solyc06g068650.4.1 - 4CL - conserved
#Solyc03g117870.3.1 - 4CL - conserved
#Solyc01g096670.4.1 - C3H -  conserved
#Solyc06g150137.1.1 - C4H - conserved
#Solyc02g084570.4.1 - F5H - conserved
#Solyc01g107910.4.1 - CCOAOMT6 - not DE
#Solyc10g050160.2.1 - CCOAOMT5 - not DE
#Solyc02g093270.4.1 - CCOAOMT - conserved
#Solyc03g080180.4.1 - SlCOMT - conserved
# Manually curated lsit of impt solaum genes
impt_genes <- c("Solyc02g077920.4.1"="SPL-CNR",
                "Solyc10g006880.3.1"="NAC-NOR",
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
                "Solyc09g075440.4.1"="NR/ETR3",
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

#Get Nobt orthologs
orthogroups <- read.table("Orthofinder/OrthoFinder/Results_Oct29/Orthogroups/Orthogroups.tsv",
                          sep="\t",
                          stringsAsFactors = F,
                          header=T)
orthogroups <- orthogroups %>% filter_all(all_vars(!grepl(',',.))) #Remove multiples
All_Genes <- cbind(Solanum=names(impt_genes), Solanum_Abbr=impt_genes)
All_Genes <- merge(All_Genes, orthogroups, by.x="Solanum", by.y="Solanum", all.x=TRUE)
#Sys.setlocale(category = "LC_ALL", locale ="en_US.UTF-8")
All_Genes <- All_Genes[order(All_Genes$Solanum_Abbr),]
row.names(All_Genes) <- NULL
rm(orthogroups)

# Manually add in genes
tmp <- c("NIOBTv3_g17210.t1" = "Solyc06g051800.3.1",
         "NIOBTv3_g28929-D2.t1" = "Solyc06g069430.3.1",
         "NIOBTv3_g39464.t1" = "Solyc03g114830.3.1",
         "NIOBTv3_g15806.t1" = "Solyc12g038510.2.1",
         "NIOBTv3_g07845.t1" = "Solyc02g065730.2.1",
         "NIOBT_gMBP20.t1" = "Solyc02g089210.4.1",
         "NIOBTv3_g08302.t1" = "Solyc10g006880.3.1")
All_Genes[All_Genes$Solanum %in%tmp,"Nicotiana"] <- names(tmp)
rm(tmp)
tmp <- c('NoACO4'="NIOBTv3_g13660.t1",
         'NoACO5'="NIOBTv3_g38689.t1",
         'NoACO6'="NIOBTv3_g02352.t1",
         'NoAGL11'="NIOBTv3_g14436.t1",
         'NoCel2'="NIOBTv3_g19880.t1",
         'NoCel3'="NIOBTv3_g12440.t1",
         'NoEXP'="NIOBTv3_g17210.t1",
         'NoFUL1'="NIOBTv3_g28929-D2.t1",
         'NoFUL2'="NIOBTv3_g39464.t1",
         'NoFYFL'="NIOBTv3_g10096.t1",
         'NoGAD1'="NIOBTv3_g11084.t1",
         'NoJ2'="NIOBTv3_g15806.t1",
         'NoMBP10'="NIOBTv3_g07845.t1",
         'NoMBP20'="NIOBT_gMBP20.t1",
         'NoMC'="NIOBTv3_g18077.t1",
         'NoNOR'="NIOBTv3_g08302.t1",
         'NoETR3'="NIOBTv3_g10291.t1",
         'NoPSY1'="NIOBTv3_g17569.t1",
         'NoFKD'="NIOBTv3_g12218.t1",
         'NoFSR'="NIOBTv3_g37943.t1",
         'NoKFB'="NIOBTv3_g35607.t1",
         'NIOBTv3_g22270.t1'="NIOBTv3_g22270.t1",
         'NIOBTv3_g10008.t1'="NIOBTv3_g10008.t1",
         'NIOBTv3_g12238.t1'="NIOBTv3_g12238.t1", 
         'NIOBTv3_g11662.t1'="NIOBTv3_g11662.t1",
         'NoSPL-CNR'="NIOBTv3_g27953.t1",
         'NoAG'="NIOBTv3_g22632-D2.t1",
         'NoSHP'="NIOBTv3_g13969.t1",
         'NoSEP1'="NIOBTv3_g14235.t1")
All_Genes[All_Genes$Nicotiana %in% tmp,"Nicotiana_Abbr"] <- names(tmp)
rm(tmp)
# Create a dummy column for whether this gene is actually used
All_Genes$Type <- TRUE
All_Genes$Type[39:41] <- FALSE
write.table(All_Genes,
            file="Tables/Orthologs.csv",
            row.names = F,
            quote = F,
            sep = ",")

# # Check if we can do ortholog stuff for 5-species data
# deleteme <- read.table("Orthofinder/OrthoFinder/Results_Aug31/Orthogroups/Orthogroups.tsv",
#                   sep="\t",
#                   stringsAsFactors = F,
#                   header=T) %>% 
#   filter_all(all_vars(!grepl(',',.))) %>%
#   filter_all(all_vars(!grepl("^$",.)))
# deleteme <- merge(deleteme,
#              All_Genes[All_Genes$Type,],
#              by.x="Solanum",
#              by.y="Solanum")

# Ortholog Table ----------------------------------------------------------
# Make table for Nobt to Slyc orthos
All_Genes[All_Genes$Type,-c(3,6)] %>%
  gt( ) %>%
  tab_header(
    title = "Orthology and Abbreviations") %>%
  tab_style(
    style=list(cell_text(style="italic")),
    locations=cells_body()) %>%
  cols_align(
    align="left") %>%
  cols_label(
    Solanum_Abbr = md("Gene Name"),
    Solanum = md("*Solanum* Gene ID"),
    Nicotiana_Abbr = md("*Nicotiana* Abbreviation"),
    Nicotiana = md("*Nicotiana* Ortholog ID")) %>%
  cols_move_to_start(
    columns=c(2,1,3,4)) %>%
  fmt_missing(
    columns=c(1:4)) %>%
  tab_footnote(
    footnote = "Only one-to-one and many-to-one orthologs",
    locations = cells_column_labels(3)) %>%
  gtsave(filename = "Orthologs.rtf",
         path = "Tables")

#Insert a list of the 78 conserved genes and find a way to include their TAIR abbreviations
orthogroups_five <- read.table("Orthofinder/OrthoFinder/Results_Aug31/Orthogroups/Orthogroups.tsv",
                               sep="\t",
                               stringsAsFactors = F,
                               header=T)
orthogroups_noise <- orthogroups_five[orthogroups_five$Orthogroup %in% rownames(subset(results(DDS_Noise), padj<0.01)),]
orthogroups_noise <- orthogroups_noise[order(orthogroups_noise$Arabidopsis),]
write.table(orthogroups_noise,
            file="Tables/Orthologs_Conserved.csv",
            row.names = F,
            quote = F,
            sep = ",")

orthogroups_noise[,-1] %>%
  gt( ) %>%
  tab_header(
    title = "Genes with Conserved Expression") %>%
  tab_style(
    style=list(cell_text(style="italic")),
    locations=cells_body()
  ) %>%
  cols_align(
    align="left"
  ) %>%
  cols_label(
    Arabidopsis = md("*A.&nbsp;thaliana*"),
    Cucumis = md("Melon"),
    Nicotiana = md("Desert Tobacco"),
    Solanum = md("Tomato")) %>%
  gtsave(filename = "Orthologs_Conserved.png", 
         path = "/bigdata/littlab/arajewski/FULTranscriptomes/Tables")

orthogroups_fruit <- orthogroups_five[orthogroups_five$Orthogroup %in% rownames(subset(results(DDS_Fruit), padj<0.01)),]
write.table(orthogroups_fruit,
            file="Tables/Orthologs_Divergent.csv",
            row.names = F,
            quote = F,
            sep = ",")


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
Tables_Nobt <- lapply(seq_along(unique(Cluster_Nobt$normalized$cluster)), 
                      function(i)  GOEnrich(gene2go = "DEGAnalysis/Pfam/Nobt.gene2go.tsv",
                                            GOIs=setNames(as.factor(as.numeric(Cluster_Nobt$normalized$genes %in% subset(Cluster_Nobt$normalized, cluster==unique(Cluster_Nobt$normalized$cluster)[i])$genes)),
                                                          Cluster_Nobt$normalized$genes),
                                            NumCategories = 10))

Tables_Nobt[["Overall"]] <- GOEnrich(gene2go = "DEGAnalysis/Pfam/Nobt.gene2go.tsv",
                                             GOIs=setNames(as.factor(as.numeric(rownames(DDS_Nobt) %in% rownames(subset(results(DDS_Nobt), padj<0.01)))), rownames(DDS_Nobt)))

# Create a list of the GO Plots
GO_Nobt <- list()
GO_Nobt <- lapply(seq_along(Tables_Nobt), 
                  function(i) GOPlot(Tables_Nobt[[i]],
                                     Title = ifelse(names(Tables_Nobt[i])=='',
                                                    paste("Cluster", i, "GO Enrichment"),
                                                    paste(names(Tables_Nobt[i]), "GO Enrichment")),
                                     colorHex = palw2[1]) +
                    theme(legend.position = c(.9,.3)))


# Tobacco Gene Plots ------------------------------------------------------

#Test if I should plot a gene at all
Test_Nobt <- as.data.frame(results(DDS_Nobt))
Test_Nobt <- Test_Nobt[row.names(Test_Nobt) %in% All_Genes$Nicotiana,]
Test_Nobt <- merge(Test_Nobt, 
                        All_Genes[,c("Nicotiana", "Nicotiana_Abbr")],
                        by.x=0,
                        by.y="Nicotiana")
Test_Nobt <- Test_Nobt[order(Test_Nobt$Nicotiana_Abbr),]
Test_Nobt$Plot <- TRUE
Test_Nobt$Plot[Test_Nobt$padj>0.01 | is.na(Test_Nobt$padj)] <- FALSE
Test_Nobt <- Test_Nobt[,-c(2:6)]


# Make a new dataframe with just those genes
# Subset to just the important genes
Subset_Nobt <- as.data.frame(assay(subset(DDS_Nobt, rownames(DDS_Nobt) %in% All_Genes$Nicotiana)))
Subset_Nobt <- merge(Subset_Nobt,
                     All_Genes[,c("Nicotiana", "Nicotiana_Abbr")],
                     by.x=0,
                     by.y="Nicotiana")
Subset_Nobt <- Subset_Nobt %>% dplyr::rename(Gene=Row.names, Abbr=Nicotiana_Abbr)
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
names(G1) <- unique(Subset_Nobt$Abbr)

# Tobacco Figures ----------------------------------------------------------
### Figure 3
(wrap_elements(GO_Nobt[[7]] & theme(legend.position = c(0.9,0.2))) /
   (C1[[1]] | C1[[2]]) /
   (C1[[3]] | C1[[4]]) /
   (C1[[5]] | C1[[6]])) +
  plot_annotation(tag_levels = "A") +
  plot_layout(nrow=4,
              heights = c(2.5,1,1,1))
ggsave2("Figures/Figure 3.pdf",
        height=13,
        width=9)

# Cluster GO supplement
(GO_Nobt[[1]] +
    GO_Nobt[[2]] +
    GO_Nobt[[3]] + 
    GO_Nobt[[4]] +
    GO_Nobt[[5]] +
    GO_Nobt[[6]]) *
  theme(legend.position = c(0.9,0.4)) + 
  plot_annotation(tag_levels = "A") +
  plot_layout(nrow=3,
              byrow = TRUE)
ggsave("Figures/Suppl_Tobacco_GO.pdf", 
       height=9, 
       width=15)


 (G1[[5]] + G1[[6]] + G1[[7]] + G1[[8]] + G1[[9]] +
   G1[[10]] + G1[[11]] + G1[[12]] + G1[[13]] + G1[[16]] +
   G1[[17]] + G1[[19]] + G1[[20]] +  G1[[23]] + G1[[24]] + 
   G1[[25]] +  G1[[26]] + G1[[27]] + G1[[28]] + G1[[29]] ) +
  plot_annotation(tag_levels = "A") +
  plot_layout(nrow=4)
ggsave2("Figures/Tobacco_Genes.pdf",
        height=10,
        width=15)


# Solanum Clusters ------------------------------------------------------
# For common patterns
C2 <- list()
C2 <- lapply(seq_along(unique(Cluster_Solanum_Noise$normalized$cluster)),
             function(i) ggplot(subset(Cluster_Solanum_Noise$normalized, 
                                       cluster==unique(Cluster_Solanum_Noise$normalized$cluster)[i]), 
                                aes(x=DAP, y=value, col=Species, fill=Species)) +
               labs(y="Z-score of Expression",
                    x="Stage") +
               scale_fill_manual(values=palw2[3]) +
               scale_color_manual(values=palw2[3]) +
               scale_x_discrete(labels=Labs_SolanumStage) +
               theme(plot.title = element_text(hjust = 0.5),
                     plot.subtitle = element_text(hjust=0.5),
                     strip.background = element_rect(fill="#FFFFFF"),
                     legend.position = "none") +
               geom_violin(position=position_dodge(width=0), alpha=0.5) +
               stat_summary(fun=mean, geom="line", aes(group=Fruit))+
               ggtitle(paste("Cluster", i, "Expression")))
# For divergent patterns
C3 <- list()
C3 <- lapply(seq_along(unique(Cluster_Solanum$normalized$cluster)),
             function(i) ggplot(subset(Cluster_Solanum$normalized, 
                                       cluster==unique(Cluster_Solanum$normalized$cluster)[i]), 
                                aes(x=DAP, y=value, col=Species, fill=Species)) +
               labs(y="Z-score of Expression",
                    x="Stage") +
               scale_fill_manual(values=palw2[c(2,3)]) +
               scale_color_manual(values=palw2[c(2,3)]) +
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
Tables_SolanumNoise <- lapply(seq_along(unique(Cluster_Solanum_Noise$normalized$cluster)), 
                             function(i)  GOEnrich(gene2go = "DEGAnalysis/Pfam/Slyc.gene2go.tsv",
                                                   GOIs=setNames(as.factor(as.numeric(Cluster_Solanum_Noise$normalized$genes %in% subset(Cluster_Solanum_Noise$normalized, cluster==unique(Cluster_Solanum_Noise$normalized$cluster)[i])$genes)),
                                                                 Cluster_Solanum_Noise$normalized$genes)))
# Include the overall, unclustered genes
Tables_SolanumNoise[["Conserved Genes"]] <- GOEnrich(gene2go = "DEGAnalysis/Pfam/Slyc.gene2go.tsv",
                                            GOIs=setNames(as.factor(as.numeric(rownames(DDS_SolanumNoise) %in% rownames(subset(results(DDS_SolanumNoise), padj<0.01)))),
                                                          rownames(DDS_SolanumNoise)))

# Create a list of the GO Plots
GO_SolanumNoise <- list()
GO_SolanumNoise <- lapply(seq_along(Tables_SolanumNoise), 
       function(i) GOPlot(Tables_SolanumNoise[[i]],
                Title = ifelse(names(Tables_SolanumNoise[i])=='',
                               paste("Cluster", i, "GO Enrichment"),
                               paste(names(Tables_SolanumNoise[i]), "GO Enrichment")),
                colorHex = "#B40F20CC",
                LegendLimit = max(sapply(Tables_SolanumNoise, function(x) max(x$Significant)))) +
                           theme(legend.position = c(.9,.2)))


# For divergent patterns
Tables_Solanum <- list()
Tables_Solanum <- lapply(seq_along(unique(Cluster_Solanum$normalized$cluster)), 
       function(i)  GOEnrich(gene2go = "DEGAnalysis/Pfam/Slyc.gene2go.tsv",
                             GOIs=setNames(as.factor(as.numeric(Cluster_Solanum$normalized$genes %in% subset(Cluster_Solanum$normalized, cluster==unique(Cluster_Solanum$normalized$cluster)[i])$genes)),
                                           Cluster_Solanum$normalized$genes)))
# Include the overall, unclustered genes
Tables_Solanum[["Divergent Genes"]] <- GOEnrich(gene2go = "DEGAnalysis/Pfam/Slyc.gene2go.tsv",
                                        GOIs=setNames(as.factor(as.numeric(rownames(DDS_Solanum) %in% rownames(subset(results(DDS_Solanum), padj<0.01)))), rownames(DDS_Solanum)))


# Create a list of the GO Plots
GO_Solanum <- list()
GO_Solanum <- lapply(seq_along(Tables_Solanum), 
                 function(i) GOPlot(Tables_Solanum[[i]],
                        Title = ifelse(names(Tables_Solanum[i])=='',
                                       paste("Cluster", i, "GO Enrichment"),
                                       paste(names(Tables_Solanum[i]), "GO Enrichment")),
                        colorHex = "#B40F20CC",
                        LegendLimit = max(sapply(Tables_Solanum, function(x) max(x$Significant)))) +
                   theme(legend.position = c(.9,.2)))


# Solanum Genes --------------------------------------------------
# Determine if a plot should be combined by species or now
Test_Solanum <- merge(as.data.frame(results(DDS_Solanum)),
                         as.data.frame(results(DDS_SolanumNoise)),
                         by=0)
Test_Solanum <- Test_Solanum[Test_Solanum$Row.names %in% All_Genes$Solanum,]
Test_Solanum <- merge(Test_Solanum, 
                   All_Genes[,c("Solanum", "Solanum_Abbr")],
                   by.x="Row.names",
                   by.y="Solanum")
Test_Solanum <- Test_Solanum[,-c(2:6,8:12)]
Test_Solanum$Choose[Test_Solanum$padj.x<0.01] <- "Separate"
Test_Solanum$Choose[Test_Solanum$padj.y<0.01 & is.na(Test_Solanum$Choose)] <- "Together"
Test_Solanum$Choose[is.na(Test_Solanum$Choose)] <- "Neither"
Test_Solanum <- Test_Solanum[order(Test_Solanum$Solanum_Abbr),]
row.names(Test_Solanum) <- NULL

# Subset to just the important genes
Subset_SolanumNoise <- as.data.frame(assay(subset(DDS_SolanumNoise, rownames(DDS_SolanumNoise) %in% All_Genes$Solanum)))
Subset_SolanumNoise <- merge(Subset_SolanumNoise,
                     All_Genes[,c("Solanum", "Solanum_Abbr")],
                     by.x=0,
                     by.y="Solanum")
Subset_SolanumNoise <- Subset_SolanumNoise %>% dplyr::rename(Gene=Row.names, Abbr=Solanum_Abbr)
Subset_SolanumNoise <- Subset_SolanumNoise[order(Subset_SolanumNoise$Abbr),]
Subset_SolanumNoise <- melt(Subset_SolanumNoise, id.vars = c("Gene", "Abbr"))
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
                                  face = ifelse(Test_Solanum[i,"Choose"]=="Together",
                                                'bold.italic',
                                                'italic')),
                     plot.subtitle = element_text(hjust=0.5),
                     strip.background = element_rect(fill="#FFFFFF"),
                     strip.text.x = element_text(face="italic"),
                     legend.position = "none") +
               geom_violin(position=position_dodge(width=0), alpha=0.5) +
               stat_summary(fun=mean, geom="line", group="Abbr")+
               ggtitle(unique(Subset_SolanumNoise$Abbr)[i]))
names(G2) <- unique(Subset_SolanumNoise$Abbr)

# Subset to just the important genes
Subset_Solanum <- as.data.frame(assay(subset(DDS_Solanum, rownames(DDS_Solanum) %in% All_Genes$Solanum)))
Subset_Solanum <- merge(Subset_Solanum,
                             All_Genes[,c("Solanum", "Solanum_Abbr")],
                             by.x=0,
                             by.y="Solanum")
Subset_Solanum <- Subset_Solanum %>% dplyr::rename(Gene=Row.names, Abbr=Solanum_Abbr)
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
names(G3) <- unique(Subset_Solanum$Abbr)


# Solanum Figures -------------------------------------------------------------
# Conserved Clusters
((C2[[1]] | GO_SolanumNoise[[1]]) /  
   (C2[[2]] | GO_SolanumNoise[[2]]) /  
   (C2[[3]] | GO_SolanumNoise[[3]]) /  
   (C2[[4]] | GO_SolanumNoise[[4]]) /  
   (C2[[5]] | GO_SolanumNoise[[5]]) /  
   (C2[[6]] | GO_SolanumNoise[[6]]) /  
   (C2[[7]] | GO_SolanumNoise[[7]]) /  
   (C2[[8]] | GO_SolanumNoise[[8]]) /  
   (C2[[9]] | GO_SolanumNoise[[9]]) /  
   (C2[[10]] | GO_SolanumNoise[[10]]) /  
   (C2[[11]] | GO_SolanumNoise[[11]]) /  
   (C2[[12]] | GO_SolanumNoise[[12]]) /  
   (C2[[13]] | GO_SolanumNoise[[13]]) /  
   (C2[[14]] | GO_SolanumNoise[[14]]) /  
   (C2[[15]] | GO_SolanumNoise[[15]]) /  
   (C2[[16]] | GO_SolanumNoise[[16]]) /  
   (C2[[17]] | GO_SolanumNoise[[17]]) /  
   (C2[[18]] | GO_SolanumNoise[[18]]) /  
   (C2[[19]] | GO_SolanumNoise[[19]]) /  
   (C2[[20]] | GO_SolanumNoise[[20]]) /
   guide_area()) +
  plot_layout(guides="collect")
ggsave2("Figures/Suppl_TomatoNoise_Clusters.pdf",
        height=70,
        width=10,
        limitsize = FALSE)

#Divergent Clusters
((C3[[1]] | GO_Solanum[[1]]) / 
    (C3[[2]] | GO_Solanum[[2]]) / 
    (C3[[3]] | GO_Solanum[[3]]) / 
    (C3[[4]] | GO_Solanum[[4]]) / 
    (C3[[5]] | GO_Solanum[[5]]) / 
    (C3[[6]] | GO_Solanum[[6]]) / 
    (C3[[7]] | GO_Solanum[[7]]) / 
    (C3[[8]] | GO_Solanum[[8]]) / 
    (C3[[9]] | GO_Solanum[[9]]) / 
    (C3[[10]] | GO_Solanum[[10]]) /     
    (C3[[11]] | GO_Solanum[[11]]) / 
    (C3[[12]] | GO_Solanum[[12]]) / 
    (C3[[13]] | GO_Solanum[[13]]) / 
    (C3[[14]] | GO_Solanum[[14]]) / 
    (C3[[15]] | GO_Solanum[[15]]) +
    guide_area()) +
  plot_layout(guides="collect")
ggsave2("Figures/Suppl_Tomato_Clusters.pdf",
        height=49,
        width=10)

### Figure 1
# Overall GO and selected Clusters
# align better with wrap_elements, and overwrite the scale for GO Plots
wrap_elements(full=GO_SolanumNoise[[21]] & 
     scale_fill_gradientn(colours = c("#87868140", "#B40F20CC"), 
                          limits=c(min(Tables_SolanumNoise[[21]]$Significant),
                                   max(Tables_SolanumNoise[[21]]$Significant)),
                          breaks=c(min(Tables_SolanumNoise[[21]]$Significant),
                                   max(Tables_SolanumNoise[[21]]$Significant)))) /
  (C2[[4]] | GO_SolanumNoise[[4]] &
     scale_fill_gradientn(colours = c("#87868140", "#B40F20CC"),
                          limits=c(min(Tables_SolanumNoise[[4]]$Significant),
                                   max(Tables_SolanumNoise[[4]]$Significant)),
                          breaks=c(min(Tables_SolanumNoise[[4]]$Significant),
                                   max(Tables_SolanumNoise[[4]]$Significant)))) /
  (C2[[10]] | GO_SolanumNoise[[10]] &
     scale_fill_gradientn(colours = c("#87868140", "#B40F20CC"), 
                          limits=c(min(Tables_SolanumNoise[[10]]$Significant),
                                   max(Tables_SolanumNoise[[10]]$Significant)),
                          breaks=c(min(Tables_SolanumNoise[[10]]$Significant),
                                   max(Tables_SolanumNoise[[10]]$Significant)))) /
  wrap_elements(full=(GO_Solanum[[16]] & 
     scale_fill_gradientn(colours = c("#87868140", "#B40F20CC"), 
                          limits=c(min(Tables_Solanum[[16]]$Significant),
                                   max(Tables_Solanum[[16]]$Significant)),
                          breaks=c(min(Tables_Solanum[[16]]$Significant),
                                   max(Tables_Solanum[[16]]$Significant))))) +
  plot_annotation(tag_levels = "A") + 
  plot_layout(widths=1, heights=1)
ggsave2("Figures/Figure 1.pdf",
        height=20,
        width=15)

#Structural Genes
(G2[["ACO1"]] + 
    G2[["ACO2"]] + 
    G3[["ACO3"]] + 
    G2[["ACO4"]] + 
    G2[["ACO5"]] + 
    G3[["ACO6"]] + 
    G3[["ACO7"]] + 
    G2[["ACS2"]] +
    G2[["ACS4"]] + 
    G3[["EIL1"]] +  
    G3[["EIL2"]] + 
    G2[["EIL4"]] +
    G2[["NR/ETR3"]] + plot_spacer() +plot_spacer() +
    G2[["CHS-1"]] + 
    G2[["CHS-2"]] + 
    G2[["CTOMT1"]] + 
    G3[["GAD1"]] + 
    G2[["GAD2"]] + 
    G2[["GAD3"]] +
    G2[["PSY1"]] + 
    G3[["TOMLOXC"]] +
    guide_area()) +
  plot_layout(guides="collect",
              nrow=5) +
  plot_annotation(tag_levels = "A")
ggsave2("Figures/Suppl_Tomato_Structural.pdf",
        height=14,
        width=14)


#Developmental TFs
(G3[["AGL11"]] + 
    G3[["EJ2/MADS1"]] + 
    G2[["FUL1"]] + 
    G3[["FUL2"]] + 
    G2[["FYFL"]] + 
    G3[["J"]] +
    G3[["J2"]] + 
    G2[["MADS-RIN"]] + 
    G3[["MBP3"]] + 
    G3[["MBP10"]] + 
    G3[["MBP20"]] + 
    G2[["NAC-NOR"]] + 
    G3[["SPL-CNR"]] + 
    G3[["STM3"]] + 
    G3[["TAG1"]] + 
    G3[["TAGL1"]] + 
    G2[["TM3"]] + 
    G3[["TM5"]] + 
    guide_area())  +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides="collect",
              nrow=4) 
ggsave2("Figures/Suppl_Tomato_Regulatory.pdf",
        height=11.2,
        width=14)

### Figure 2
# Selected Divergent Genes
((G3[["ACO6"]] +
    G3[["TOMLOXC"]] +
    G3[["GAD1"]] +
    G3[["MBP3"]] + 
    G3[["SPL-CNR"]] + 
    G3[["TAG1"]] + 
    G3[["TAGL1"]] +
    guide_area()) *
    scale_color_manual(values=palw2[c(2,3)],labels=c("Wild", "Cultivated")) *
    scale_fill_manual(values=palw2[c(2,3)],labels=c("Wild", "Cultivated"))) +
  plot_annotation(tag_levels = "A") + 
  plot_layout(guides="collect")
ggsave2("Figures/Figure 2.pdf",
        height=8.4,
        width=8.4)

 
# Solanaceae Clusters -----------------------------------------------------
C6 <- list()
C6 <- lapply(seq_along(unique(Cluster_Solanaceae_Noise$normalized$cluster)),
             function(i) ggplot(subset(Cluster_Solanaceae_Noise$normalized, 
                                       cluster==unique(Cluster_Solanaceae_Noise$normalized$cluster)[i]), 
                                aes(x=Stage, y=value, col=Species, fill=Species)) +
               labs(y="Z-score of Expression",
                    x="Stage") +
               scale_fill_manual(values=palw2[1]) +
               scale_color_manual(values=palw2[3]) +
               scale_x_discrete(labels=Labs_Stage) +
               theme(plot.title = element_text(hjust = 0.5),
                     plot.subtitle = element_text(hjust=0.5),
                     strip.background = element_rect(fill="#FFFFFF"),
                     legend.position = "none") +
               geom_violin(position=position_dodge(width=0), alpha=0.5) +
               stat_summary(fun=mean, geom="line", aes(group=Fruit))+
               ggtitle(paste("Cluster", i, "Expression")))
# For divergent patterns
C7 <- list()
C7 <- lapply(seq_along(unique(Cluster_Solanaceae$normalized$cluster)),
             function(i) ggplot(subset(Cluster_Solanaceae$normalized, 
                                       cluster==unique(Cluster_Solanaceae$normalized$cluster)[i]), 
                                aes(x=Stage, y=value, col=Fruit, fill=Fruit)) +
               labs(y="Z-score of Expression",
                    x="Stage") +
               scale_fill_manual(values=palw2[c(1,3)]) +
               scale_color_manual(values=palw2[c(1,3)]) +
               scale_x_discrete(labels=Labs_Stage) +
               theme(plot.title = element_text(hjust = 0.5),
                     plot.subtitle = element_text(hjust=0.5),
                     strip.background = element_rect(fill="#FFFFFF")) +
               geom_violin(position=position_dodge(width=0), alpha=0.5) +
               stat_summary(fun=mean, geom="line", aes(group=Fruit))+
               ggtitle(paste("Cluster", i, "Expression")))

# Solanaceae GO Plots --------------------------------------------------------
# For Common patterns
Tables_SolanaceaeNoise <- list()
Tables_SolanaceaeNoise <- lapply(seq_along(unique(Cluster_Solanaceae_Noise$normalized$cluster)), 
                              function(i)  GOEnrich(gene2go = "DEGAnalysis/Pfam/Ortho.1029Sol.gene2go.tsv",
                                                    GOIs=setNames(as.factor(as.numeric(Cluster_Solanaceae_Noise$normalized$genes %in% subset(Cluster_Solanaceae_Noise$normalized, cluster==unique(Cluster_Solanaceae_Noise$normalized$cluster)[i])$genes)),
                                                                  Cluster_Solanaceae_Noise$normalized$genes)))

# Include the overall, unclustered genes
Tables_SolanaceaeNoise[["Conserved Genes"]] <- GOEnrich(gene2go = "DEGAnalysis/Pfam/Ortho.1029Sol.gene2go.tsv",
                                                     GOIs=setNames(as.factor(as.numeric(rownames(DDS_Solanaceae_Noise) %in% rownames(subset(results(DDS_Solanaceae_Noise), padj<0.01)))), rownames(DDS_Solanaceae_Noise)))

# Create a list of the GO Plots
GO_SolanaceaeNoise <- list()
GO_SolanaceaeNoise <- lapply(seq_along(Tables_SolanaceaeNoise), 
                          function(i) GOPlot(Tables_SolanaceaeNoise[[i]],
                                             Title = ifelse(names(Tables_SolanaceaeNoise[i])=='',
                                                            paste("Cluster", i, "GO Enrichment"),
                                                            paste(names(Tables_SolanaceaeNoise[i]), "GO Enrichment")),
                                             colorHex = "#B40F20CC") +
                            theme(legend.position = c(.9,.2)))


# For divergent patterns
Tables_Solanaceae <- list()
Tables_Solanaceae <- lapply(seq_along(unique(Cluster_Solanaceae$normalized$cluster)), 
                         function(i)  GOEnrich(gene2go = "DEGAnalysis/Pfam/Ortho.1029Sol.gene2go.tsv",
                                               GOIs=setNames(as.factor(as.numeric(Cluster_Solanaceae$normalized$genes %in% subset(Cluster_Solanaceae$normalized, cluster==unique(Cluster_Solanaceae$normalized$cluster)[i])$genes)),
                                                             Cluster_Solanaceae$normalized$genes)))
# Include the overall, unclustered genes
Tables_Solanaceae[["Divergent Genes"]] <- GOEnrich(gene2go = "DEGAnalysis/Pfam/Ortho.1029Sol.gene2go.tsv",
                                                GOIs=setNames(as.factor(as.numeric(rownames(DDS_Solanaceae) %in% rownames(subset(results(DDS_Solanaceae), padj<0.01)))), rownames(DDS_Solanaceae)))

# Create a list of the GO Plots
GO_Solanaceae <- list()
GO_Solanaceae <- lapply(seq_along(Tables_Solanaceae)[-c(which(sapply(Tables_Solanaceae,is.null)))], 
                     function(i) GOPlot(Tables_Solanaceae[[i]],
                                        Title = ifelse(names(Tables_Solanaceae[i])=='',
                                                       paste("Cluster", i, "GO Enrichment"),
                                                       paste(names(Tables_Solanaceae[i]), "GO Enrichment")),
                                        colorHex = "#B40F20CC") +
                       theme(legend.position = c(.9,.2)))



# Solanaceae Genes --------------------------------------------------------
# Determine if a plot should be combined by species or now
Test_Solanaceae <- merge(as.data.frame(results(DDS_Solanaceae)),
                      as.data.frame(results(DDS_Solanaceae_Noise)),
                      by=0)
Test_Solanaceae <- Test_Solanaceae[Test_Solanaceae$Row.names %in% All_Genes$Orthogroup,]
Test_Solanaceae <- merge(Test_Solanaceae, 
                      All_Genes[,c("Orthogroup", "Solanum_Abbr")],
                      by.x="Row.names",
                      by.y="Orthogroup")
Test_Solanaceae <- Test_Solanaceae[,-c(2:6,8:12)]
Test_Solanaceae$Choose[Test_Solanaceae$padj.x<0.01] <- "Separate"
Test_Solanaceae$Choose[Test_Solanaceae$padj.y<0.01 & is.na(Test_Solanaceae$Choose)] <- "Together"
Test_Solanaceae$Choose[is.na(Test_Solanaceae$Choose)] <- "Neither"
Test_Solanaceae <- Test_Solanaceae[order(Test_Solanaceae$Solanum_Abbr),]
row.names(Test_Solanaceae) <- NULL

# Subset to just the important genes
Subset_Solanaceae_Noise <- as.data.frame(assay(subset(DDS_Solanaceae_Noise, rownames(DDS_Solanaceae_Noise) %in% All_Genes$Orthogroup)))
Subset_Solanaceae_Noise <- merge(Subset_Solanaceae_Noise,
                             All_Genes[,c("Orthogroup", "Solanum_Abbr")],
                             by.x=0,
                             by.y="Orthogroup")
Subset_Solanaceae_Noise <- Subset_Solanaceae_Noise %>% dplyr::rename(Gene=Row.names, Abbr=Solanum_Abbr)
Subset_Solanaceae_Noise <- Subset_Solanaceae_Noise[order(Subset_Solanaceae_Noise$Abbr),]
Subset_Solanaceae_Noise <- melt(Subset_Solanaceae_Noise, id.vars = c("Gene", "Abbr"))
Subset_Solanaceae_Noise$Stage <- as.factor(rep(colData(DDS_Solanaceae_Noise)$Stage,
                                         each=length(unique(Subset_Solanaceae_Noise$Gene))))

# Genes with conserved patterns (there really arent any)
G4 <- list()
G4 <- lapply(seq_along(unique(Subset_Solanaceae_Noise$Abbr)),
             function(i) ggplot(Subset_Solanaceae_Noise[Subset_Solanaceae_Noise$Abbr==unique(Subset_Solanaceae_Noise$Abbr)[i],], 
                                aes(x=Stage, y=value, col=Abbr, fill=Abbr)) +
               labs(y="Normalized Counts",
                    x="Stage") +
               scale_color_manual(values=palw2[3]) +
               scale_fill_manual(values=palw2[1]) +
               scale_x_discrete(labels=Labs_Stage) +
               theme(plot.title = element_text(hjust = 0.5,
                                               face = ifelse(Test_Solanaceae[i,"Choose"]=="Together",
                                                             'bold.italic',
                                                             'italic')),
                     plot.subtitle = element_text(hjust=0.5),
                     strip.background = element_rect(fill="#FFFFFF"),
                     strip.text.x = element_text(face="italic"),
                     legend.position = "none") +
               geom_violin(position=position_dodge(width=0), alpha=0.5) +
               stat_summary(fun=mean, geom="line", group="Abbr")+
               ggtitle(unique(Subset_Solanaceae_Noise$Abbr)[i]))
names(G4) <- unique(Subset_Solanaceae_Noise$Abbr)

# Subset to just the important genes
Subset_Solanaceae <- as.data.frame(assay(subset(DDS_Solanaceae, rownames(DDS_Solanaceae) %in% 
                                                  All_Genes$Orthogroup)))
Subset_Solanaceae <- merge(Subset_Solanaceae,
                        All_Genes[,c("Orthogroup", "Solanum_Abbr")],
                        by.x=0,
                        by.y="Orthogroup")
Subset_Solanaceae <- Subset_Solanaceae %>% dplyr::rename(Gene=Row.names, Abbr=Solanum_Abbr)
Subset_Solanaceae <- Subset_Solanaceae[order(Subset_Solanaceae$Abbr),]
Subset_Solanaceae <- melt(Subset_Solanaceae, id.vars = c("Gene", "Abbr"))
Subset_Solanaceae$Fruit <- rep(colData(DDS_Solanaceae)$Fruit, 
                              each=length(unique(Subset_Solanaceae$Gene)))
Subset_Solanaceae$Stage <- as.factor(rep(colData(DDS_Solanaceae)$Stage,
                                    each=length(unique(Subset_Solanaceae$Gene))))
Subset_Solanaceae$Species <- rep(colData(DDS_Solanaceae)$Species,
                                 each=length(unique(Subset_Solanaceae$Gene)))

#Genes with divergent patterns between fruits
G5 <- list()
G5 <- lapply(seq_along(unique(Subset_Solanaceae$Abbr)),
             function(i) ggplot(Subset_Solanaceae[Subset_Solanaceae$Abbr==unique(Subset_Solanaceae$Abbr)[i],], 
                                aes(x=Stage, y=value, col=Fruit, fill=Fruit)) +
               labs(y="Normalized Counts",
                    x="Stage") +
               scale_color_manual(values=c("Pimpinellifolium"=palw2[2],
                                           "Dry"= palw2[1],
                                           "Fleshy"=palw2[3]),
                                  labels=c("Pimpinellifolium"="Wild Tomato",
                                           "Dry"="Desert Tobacco",
                                           "Fleshy"="Cultivated Tomato"),
                                  name="Species") +
               scale_fill_manual(values=c("Pimpinellifolium"=palw2[2],
                                          "Dry"= palw2[1],
                                          "Fleshy"=palw2[3]),
                                 labels=c("Pimpinellifolium"="Wild Tomato",
                                          "Dry"="Desert Tobacco",
                                          "Fleshy"="Cultivated Tomato"),
                                 name="Species") +
               scale_x_discrete(labels=Labs_Stage) +
               theme(plot.title = element_text(hjust = 0.5,
                                               face = ifelse(Test_Solanaceae[i,"Choose"]=="Separate", 'bold.italic', 'italic')),
                     plot.subtitle = element_text(hjust=0.5),
                     strip.background = element_rect(fill="#FFFFFF"),
                     strip.text.x = element_text(face="italic")) +
               geom_violin(position=position_dodge(width=0), alpha=0.5) +
               stat_summary(fun=mean, geom="line", aes(group=Fruit)) +
               ggtitle(unique(Subset_Solanaceae$Abbr)[i]))
names(G5) <- unique(Subset_Solanaceae$Abbr)

# Same as G5 really, but with separate plots by species to highlight differences if necessary
G6 <- list()
G6 <- lapply(seq_along(unique(Subset_Solanaceae$Abbr)),
             function(i) ggplot(Subset_Solanaceae[Subset_Solanaceae$Abbr==unique(Subset_Solanaceae$Abbr)[i],], 
                                aes(x=Stage, y=value, col=Species, fill=Species)) +
               labs(y="Normalized Counts",
                    x="Stage") +
               scale_color_manual(values=palw2[c(2,1,3)],
                                  labels=c("Tobacco"="Desert Tobacco/Dry",
                                           "Tomato"="Cultivated Tomato/Fleshy",
                                           "Pimpinellifolium"="Wild Tomato"),
                                  name="Species") +
               scale_fill_manual(values=palw2[c(2,1,3)],
                                 labels=c("Tobacco"="Desert Tobacco/Dry",
                                          "Tomato"="Cultivated Tomato/Fleshy",
                                          "Pimpinellifolium"="Wild Tomato"),
                                 name="Species") +
               scale_x_discrete(labels=Labs_Stage) +
               theme(plot.title = element_text(hjust = 0.5,
                                               face = ifelse(Test_Solanaceae[i,"Choose"]=="Separate", 'bold.italic', 'italic')),
                     plot.subtitle = element_text(hjust=0.5),
                     strip.background = element_rect(fill="#FFFFFF"),
                     strip.text.x = element_text(face="italic")) +
               geom_violin(position=position_dodge(width=0), alpha=0.5) +
               stat_summary(fun=mean, geom="line", aes(group=Species)) +
               ggtitle(unique(Subset_Solanaceae$Abbr)[i]))
names(G6) <- unique(Subset_Solanaceae$Abbr)

# Solanaceae Figures ------------------------------------------------------
### Figure 4
(wrap_elements(GO_SolanaceaeNoise[[8]] & theme(legend.position = c(0.9,0.3))) /
   (C6[[3]] + GO_SolanaceaeNoise[[3]]) /
   wrap_elements(GO_Solanaceae[[36]] & theme(legend.position = c(0.9,0.3)))) +
  plot_annotation(tag_levels = "A")
ggsave2("Figures/Figure 4.pdf",
        height=14,
        width=13)

# All Conserved Clusters
((C6[[1]] | GO_SolanaceaeNoise[[1]]) /
    (C6[[2]] | GO_SolanaceaeNoise[[2]]) /
    (C6[[3]] | GO_SolanaceaeNoise[[3]]) /
    (C6[[4]] | GO_SolanaceaeNoise[[4]]) /
    (C6[[5]] | GO_SolanaceaeNoise[[5]]) /
    (C6[[6]] | GO_SolanaceaeNoise[[6]]) /
    (C6[[7]] | GO_SolanaceaeNoise[[7]]))
ggsave2("Figures/Suppl_SolNoise_Clusters.pdf",
        height=27,
        width=12)

# God forgive me my sins for plotting this
((C7[[1]] | GO_Solanaceae[[1]] | C7[[2]] | GO_Solanaceae[[2]]) /
    (C7[[3]] | GO_Solanaceae[[3]] | C7[[4]] | GO_Solanaceae[[4]]) /
    (C7[[5]] | GO_Solanaceae[[5]] | C7[[6]] | GO_Solanaceae[[6]]) /
    (C7[[7]] | GO_Solanaceae[[7]] | C7[[8]] | GO_Solanaceae[[8]]) /
    (C7[[9]] | GO_Solanaceae[[9]] | C7[[10]] | GO_Solanaceae[[10]]) /
    (C7[[11]] | GO_Solanaceae[[11]] | C7[[12]] | GO_Solanaceae[[12]]) /
    (C7[[13]] | GO_Solanaceae[[13]] | C7[[14]] | GO_Solanaceae[[14]]) /
    (C7[[15]] | plot_spacer() | C7[[16]] | GO_Solanaceae[[15]]) /
    (C7[[17]] | GO_Solanaceae[[16]] | C7[[18]] | GO_Solanaceae[[17]]) /
    (C7[[19]] | GO_Solanaceae[[18]] | C7[[20]] | GO_Solanaceae[[19]]) /
    (C7[[21]] | GO_Solanaceae[[20]] | C7[[22]] | GO_Solanaceae[[21]]) /
    (C7[[23]] | GO_Solanaceae[[22]] | C7[[24]] | GO_Solanaceae[[23]]) /
    (C7[[25]] | GO_Solanaceae[[24]] | C7[[26]] | GO_Solanaceae[[25]]) /
    (C7[[27]] | GO_Solanaceae[[26]] | C7[[28]] | GO_Solanaceae[[27]]) /
    (C7[[29]] | GO_Solanaceae[[28]] | C7[[30]] | GO_Solanaceae[[29]]) /
    (C7[[31]] | GO_Solanaceae[[30]] | C7[[32]] | GO_Solanaceae[[31]]) /
    (C7[[33]] | GO_Solanaceae[[32]] | C7[[34]] | GO_Solanaceae[[33]]) /
    (C7[[35]] | GO_Solanaceae[[34]] | C7[[36]] | GO_Solanaceae[[35]]))
ggsave2("Figures/Suppl_Sol_Clusters.pdf",
        height=60,
        width=24,
        limitsize=FALSE)

### Figure 5
((G5[["ACO4"]] & theme(legend.position="none")) +  
    (G5[["ACO5"]] & theme(legend.position="none")) + 
    (G6[["ACO6"]] & theme(legend.position="none")) + 
    (G5[["NR/ETR3"]] & theme(legend.position="none")) + 
    plot_spacer() + 
    plot_spacer() +
    (G5[["AGL11"]] & theme(legend.position="none")) + 
    (G6[["FYFL"]]  & theme(legend.position="none")) +  
    (G5[["SPL-CNR"]]  & theme(legend.position="none")) + 
    (G6[["TAG1"]]) + 
    (G5[["TAGL1"]] & theme(legend.position="none")) +
    guide_area()) + 
  plot_annotation(tag_levels = "A") +
  plot_layout(guides="collect",
              ncol = 3)
ggsave2("Figures/Figure 5.pdf",
        height=11.2,
        width=8.4)


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

# Five-Species Clusters ---------------------------------------------------
# Conserved patterns first
C4 <- list()
C4 <- lapply(seq_along(unique(Cluster_Noise$normalized$cluster)),
             function(i) ggplot(subset(Cluster_Noise$normalized, 
                                       cluster==unique(Cluster_Noise$normalized$cluster)[i]),
                                aes(x=Stage, y=value, col=Species, fill=Species)) +
               labs(y="Z-score of Expression",
                    x="Stage") +
               scale_fill_manual(values="#9EBE91") +
               scale_color_manual(values="#9EBE91") +
               scale_x_discrete(labels=Labs_Stage) +
               theme(plot.title = element_text(hjust = 0.5),
                     plot.subtitle = element_text(hjust=0.5),
                     strip.background = element_rect(fill="#FFFFFF"),
                     legend.position = "none") +
               geom_violin(position=position_dodge(width=0), alpha=0.5) +
               stat_summary(fun=mean, geom="line", aes(group=Fruit))+
               ggtitle(paste("Cluster", i, "Expression")))

# For divergent patterns
C5 <- list()
C5 <- lapply(seq_along(unique(Cluster_Fruit$normalized$cluster)),
             function(i) ggplot(subset(Cluster_Fruit$normalized, 
                                       cluster==unique(Cluster_Fruit$normalized$cluster)[i]), 
                                aes(x=Stage, y=value, col=Fruit, fill=Fruit)) +
               labs(y="Z-score of Expression",
                    x="Stage") +
               scale_fill_manual(values=palw2[c(1,3)]) +
               scale_color_manual(values=palw2[c(1,3)]) +
               scale_x_discrete(labels=Labs_Stage) +
               theme(plot.title = element_text(hjust = 0.5),
                     plot.subtitle = element_text(hjust=0.5),
                     strip.background = element_rect(fill="#FFFFFF")) +
               geom_violin(position=position_dodge(width=0), alpha=0.5) +
               stat_summary(fun=mean, geom="line", aes(group=Species))+
               ggtitle(paste("Cluster", i, "Expression")))


# Five-Species GO Plots ---------------------------------------------------
# For Common patterns
Tables_Ortho_Noise <- list()
Tables_Ortho_Noise <- lapply(seq_along(unique(Cluster_Noise$normalized$cluster)), 
                             function(i)  GOEnrich(gene2go = "DEGAnalysis/Pfam/Ortho.831All.gene2go.tsv",
                                                   GOIs=setNames(as.factor(as.numeric(Cluster_Noise$normalized$genes %in% subset(Cluster_Noise$normalized, cluster==unique(Cluster_Noise$normalized$cluster)[i])$genes)),
                                                                 Cluster_Noise$normalized$genes)))

# For divergent patterns 
Tables_Ortho_Fruit <- list()
Tables_Ortho_Fruit <- lapply(seq_along(unique(Cluster_Fruit$normalized$cluster)), 
                       function(i)  GOEnrich(gene2go = "DEGAnalysis/Pfam/Ortho.831All.gene2go.tsv",
                                             GOIs=setNames(as.factor(as.numeric(Cluster_Fruit$normalized$genes %in% Cluster_Fruit$normalized$genes[Cluster_Fruit$normalized$cluster==unique(Cluster_Fruit$normalized$cluster)[i]])),
                                                           Cluster_Fruit$normalized$genes)))
# Then add in Overalls
Tables_Ortho_Noise[["Model 1"]] <- GOEnrich(gene2go = "DEGAnalysis/Pfam/Ortho.831All.gene2go.tsv",
                                     GOIs=setNames(as.factor(as.numeric(rownames(DDS_Noise) %in% rownames(subset(results(DDS_Noise), padj<0.01)))), rownames(DDS_Noise)))

Tables_Ortho_Fruit[["Model 2"]] <- GOEnrich(gene2go = "DEGAnalysis/Pfam/Ortho.831All.gene2go.tsv",
                                     GOIs=setNames(as.factor(as.numeric(rownames(DDS_Fruit) %in% rownames(subset(results(DDS_Fruit), padj<0.01)))), rownames(DDS_Fruit)))
                                     
# Save the GO Plots
GO_Ortho_Noise <- list()
GO_Ortho_Noise <- lapply(seq_along(Tables_Ortho_Noise), 
                         function(i) GOPlot(Tables_Ortho_Noise[[i]],
                                            Title = ifelse(names(Tables_Ortho_Noise[i])=='',
                                                           paste("Cluster", i, "GO Enrichment"),
                                                           paste(names(Tables_Ortho_Noise[i]), "GO Enrichment")),
                                            colorHex = "#9EBE91CC") +
                           theme(legend.position = c(.9,.2)))
GO_Ortho_Fruit <- list()
GO_Ortho_Fruit <- lapply(seq_along(Tables_Ortho_Fruit), 
             function(i) GOPlot(Tables_Ortho_Fruit[[i]],
                      Title = ifelse(names(Tables_Ortho_Fruit[i])=='',
                             paste("Cluster", i, "GO Enrichment"),
                             paste(names(Tables_Ortho_Fruit[i]), "GO Enrichment")),
                      colorHex = "#9EBE91CC") +
               theme(legend.position = c(.9,.2)))


# Five-Species PCA --------------------------------------------------------
Subset_Noise <- DDS_Noise[rownames(DDS_Noise) %in% rownames(subset(results(DDS_Noise), padj<=0.01))]
RLD_Noise <- rlog(Subset_Noise, blind=FALSE)
RLD_Noise$Stage[RLD_Noise$Stage=="3.5"] <- "Tr"
RLD_Noise$Stage <- as.factor(RLD_Noise$Stage)
RLD_Noise$Species <- factor(RLD_Noise$Species, levels=c("Arabidopsis", "Melon", "Tobacco", "Tomato", "Pimpinellifolium"))
# Make three PCA plots
PCA1 <- lapply(list(c(1,2),c(2,3),c(1,3), c(1,4), c(2,4), c(3,4), c(1,5), c(2,5), c(3,5), c(4,5)),
               function(i) plotPCAmod(RLD_Noise,xPC=i[1], yPC=i[2]))

Subset_Fruit <- DDS_Fruit[rownames(DDS_Fruit) %in% rownames(subset(results(DDS_Fruit), padj<=0.01))]
RLD_Fruit <- rlog(Subset_Fruit, blind=FALSE)
RLD_Fruit$Stage[RLD_Fruit$Stage=="3.5"] <- "Tr"
RLD_Fruit$Stage <- as.factor(RLD_Fruit$Stage)
RLD_Fruit$Species <- factor(RLD_Fruit$Species, levels=c("Arabidopsis", "Melon", "Tobacco", "Tomato", "Pimpinellifolium"))
PCA2 <-  lapply(list(c(1,2),c(2,3),c(1,3), c(1,4), c(2,4), c(3,4), c(1,5), c(2,5), c(3,5), c(4,5)),
                function(i) plotPCAmod(RLD_Fruit,
                                       xPC=i[1], yPC=i[2],
                                       ntop=dim(RLD_Fruit)[1]))


# Five Ortho Figures ------------------------------------------------------

### Figure 7
((PCA1[[1]] + PCA1[[6]] + PCA1[[9]] + guide_area() + plot_layout(guides = "collect")) | 
   GO_Ortho_Noise[[3]]) +
  plot_annotation(tag_levels = "A")
ggsave2("Figures/Figure 7.pdf",
        width=16,
        height=7)

### Figure 8
((PCA2[[1]] + PCA2[[2]] + PCA2[[10]] + guide_area() + plot_layout(guides="collect")) |
    GO_Ortho_Fruit[[9]]) +
  plot_annotation(tag_levels = "A")
ggsave2("Figures/Figure 8.pdf",
        width=16,
        height=7)

# Conserved clusters
(C4[[1]] + GO_Ortho_Noise[[1]] + 
    C4[[2]] + GO_Ortho_Noise[[2]]) +
  plot_annotation(tag_levels = "A") +
  plot_layout(nrow=2)
ggsave2("Figures/Suppl_Model1_Clusters.pdf",
       height=10,
       width=15)

### Figure 9
# Selected Clusters from Fruit ortho data
((C5[[4]] & theme(legend.position = c(.9,.9))) + GO_Ortho_Fruit[[4]] + 
    (C5[[6]] & theme(legend.position = "none")) + GO_Ortho_Fruit[[6]] + 
    (C5[[7]] & theme(legend.position = "none")) + GO_Ortho_Fruit[[7]] + 
    (C5[[8]] & theme(legend.position = "none")) + GO_Ortho_Fruit[[8]]) +
  plot_layout(nrow = 4) +
  plot_annotation(tag_levels = "A")
ggsave2("Figures/Figure 9.pdf",
       height=20,
       width=15)

# All cluster for supplement
((C5[[1]] & theme(legend.position = c(.9,.9))) + GO_Ortho_Fruit[[1]] + 
    (C5[[2]] & theme(legend.position = "none")) + GO_Ortho_Fruit[[2]] + 
    (C5[[3]] & theme(legend.position = "none")) + GO_Ortho_Fruit[[3]] + 
    (C5[[4]] & theme(legend.position = "none")) + GO_Ortho_Fruit[[4]] + 
    (C5[[5]] & theme(legend.position = "none")) + GO_Ortho_Fruit[[5]] + 
    (C5[[6]] & theme(legend.position = "none")) + GO_Ortho_Fruit[[6]] + 
    (C5[[7]] & theme(legend.position = "none")) + GO_Ortho_Fruit[[7]] + 
    (C5[[8]] & theme(legend.position = "none")) + GO_Ortho_Fruit[[8]]) +
  plot_layout(nrow = 8) +
  plot_annotation(tag_levels = "A")
ggsave2("Figures/Suppl_Model2_Clusters.pdf",
        height=40,
        width=15)



# Defense Figures ---------------------------------------------------------

# Solanum Conserved Genes
GO_SolanumNoise[[21]] & 
  scale_fill_gradientn(colours = c("#87868140", "#B40F20CC"), 
                       limits=c(min(Tables_SolanumNoise[[21]]$Significant),
                                max(Tables_SolanumNoise[[21]]$Significant)),
                       breaks=c(min(Tables_SolanumNoise[[21]]$Significant),
                                max(Tables_SolanumNoise[[21]]$Significant)))
ggsave2("Defense/GO_Solanum_Conserved.png",
        height=8,
        width=8)
# Solanum Divergent Genes
GO_Solanum[[16]] & 
  scale_fill_gradientn(colours = c("#87868140", "#B40F20CC"), 
                       limits=c(min(Tables_Solanum[[16]]$Significant),
                                max(Tables_Solanum[[16]]$Significant)),
                       breaks=c(min(Tables_Solanum[[16]]$Significant),
                                max(Tables_Solanum[[16]]$Significant)))
ggsave2("Defense/GO_Solanum_Divergent.png",
        height=5.3,
        width=8)

# Solanaceae Conserved Genes
GO_SolanaceaeNoise[[8]]
ggsave2("Defense/GO_Solanaceae_Conserved.png",
        height=4.3,
        width=8)

# All clusters as an example
(C6[[1]] + 
    C6[[2]] + 
    C6[[3]] + 
    C6[[4]] + 
    C6[[5]] + 
    C6[[6]] + 
    C6[[7]]) *
    ggtitle(element_blank()) +
  plot_layout(nrow=2)
ggsave2("Defense/Cluster_Solanaceae_All.png",
        height=6,
        width=13)


# Solanaceae Conserved Cluster 3 example
(C6[[3]] | GO_SolanaceaeNoise[[3]]) +
  plot_layout(widths=c(1.3,1))
ggsave2("Defense/GO_Solanaceae_Cluster3.png",
        height=6,
        width=14)

# Solanaceae Divergent GO
GO_Solanaceae[[36]] &
  theme(legend.position = c(0.8,0.2))
ggsave2("Defense/GO_Solanaceae_Divergent.png",
        height=6,
        width=8)

# Solanaceae Divergent Cluster 36 example
(C7[[36]] &
    scale_fill_manual(values=palw2[c(1,3)],
                      labels=c("Desert Tobacco", "Tomato"),
                      name="Species") &
    scale_color_manual(values=palw2[c(1,3)],
                       labels=c("Desert Tobacco", "Tomato"),
                       name="Species") &
    ggtitle(element_blank()) &
    theme(legend.position = c(0.55,0.85))) / 
    wrap_elements(GO_Solanaceae[[35]] &
                    ggtitle(element_blank()) &
                    theme(legend.position = c(0.85,0.25)))
ggsave2("Defense/GO_Solanaceae_Cluster36.png",
        height=8,
        width=8)

# Some Solanaceae Genes
(G5[["ACO4"]] & theme(legend.position="none")) + 
    (G5[["ACO5"]] & theme(legend.position="none")) + 
    (G6[["ACO6"]] & 
       theme(plot.title = element_text(face="bold.italic")) &
       scale_color_manual(values=palw2[c(2,1,3)],
                          labels=c("Tobacco"="Desert Tobacco",
                                   "Tomato"="Tomato(es)",
                                   "Pimpinellifolium"="(Wild Tomato)"),
                          name="Species") &
       scale_fill_manual(values=palw2[c(2,1,3)],
                         labels=c("Tobacco"="Desert Tobacco",
                                  "Tomato"="Tomato(es)",
                                  "Pimpinellifolium"="(Wild Tomato)"),
                         name="Species")) + 
    (G6[["TAG1"]] & theme(legend.position="none",
                          plot.title = element_text(face="bold.italic"))) +
    (G5[["TAGL1"]] & theme(legend.position="none")) +
  guide_area() +
  plot_layout(guides="collect")
ggsave2("Defense/Genes_Solanaceae.png",
        height=6,
        width=9)

# All Model 1 Clusters
(C4[[1]] + C4[[2]]) * ggtitle(element_blank())
ggsave2("Defense/Cluster_Model1.png",
        height=3,
        width=6.5)

# All Model 2 Clusters
((C5[[1]] & theme(legend.position="none")) + 
    (C5[[2]] & theme(legend.position="none")) +
    (C5[[3]] & theme(legend.position="none")) +
    (C5[[4]] & theme(legend.position="none")) +
    (C5[[5]] & theme(legend.position="none")) +
    (C5[[6]] & theme(legend.position="none")) +
    (C5[[7]] & theme(legend.position="none")) +
    (C5[[8]]) + 
    guide_area()) * 
  ggtitle(element_blank()) +
  plot_layout(guide="collect",
              nrow=2)
ggsave2("Defense/Cluster_Model2.png",
        height=6,
        width=16)

# Five species conserved
((PCA1[[1]] + PCA1[[6]] + PCA1[[9]] + guide_area() + plot_layout(guides = "collect")) | 
    GO_Ortho_Noise[[3]] & ggtitle(element_blank()))
ggsave2("Defense/PCA_Conserved.png",
        width=16,
        height=7)

# Five Species divergent
((PCA2[[1]] + PCA2[[2]] + PCA2[[10]] + guide_area() + plot_layout(guides="collect")) |
    GO_Ortho_Fruit[[9]] & ggtitle(element_blank()))
ggsave2("Defense/PCA_Divergent.png",
        width=16,
        height=7)

# explain model testing
# generate sample data
nObs=100
test.df <- data.frame(Time=rep(rep(c(1:5),each=nObs/5),2),
                      Case=rep(c("Slow", "Fast"), each=nObs))
# create two patterns
test.df$Counts <- c(sample(90:105, 100, replace = TRUE),
                   10*rep(c(1:5),each=nObs/5) + sample(50:100, 100, replace = TRUE))
test.df$RandCase <- rep(c("Fast1", "Fast2"))

# plot as one group
t1 <- ggplot(test.df, aes(x=Time, y=Counts, col=palw2[1])) +
  scale_color_manual(values="black") +
  geom_point() +
  theme(legend.position = "none")

t2 <- t1 +
  geom_smooth(method="lm",
              se=FALSE)

# plot as two groups based on real difference
t3 <- ggplot(test.df, aes(x=Time, y=Counts, col=Case, group=Case)) +
  scale_color_manual(values=palw2[1:2]) +
  geom_point()

t4 <- t3 +
  geom_smooth(method="lm",
              se=FALSE)

# plot as two groups based on nothing
t5 <- ggplot(test.df, aes(x=Time, y=Counts, col=RandCase, group=RandCase)) +
  scale_color_manual(values=palw2[1:2]) +
  geom_point()

t6 <- t5 +
  geom_smooth(method="lm",
              se=FALSE)

# Combine into one figure
(t1 + 
    t2 + 
    t3 +
    t4 +
    t5 +
    t6) * 
  theme(legend.position = "none") +
  plot_layout(nrow=3)
ggsave2("Defense/Test_Explanation.png",
        height=9,
        width=7)


# Reviewer Responses -------------------------------------------------------------

## Gene Counts
# In review, we were asked to provide gene counts for each figure panel
# I'll include them in the figure legends, but make a list here to organize them
Fig_Gene_Counts <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(Fig_Gene_Counts) <- c("Figure", "Panel", "GeneCount")
# Add in rows as they're computed
# Fig 1
Fig_Gene_Counts[nrow(Fig_Gene_Counts) + 1,] = c(1, "B", sum(as.factor(as.numeric(Cluster_Solanum_Noise$normalized$genes %in% subset(Cluster_Solanum_Noise$normalized, cluster==unique(Cluster_Solanum_Noise$normalized$cluster)[4])$genes))==1))
Fig_Gene_Counts[nrow(Fig_Gene_Counts) + 1,] = c(1, "D", sum(as.factor(as.numeric(Cluster_Solanum_Noise$normalized$genes %in% subset(Cluster_Solanum_Noise$normalized, cluster==unique(Cluster_Solanum_Noise$normalized$cluster)[10])$genes))==1))
# Fig 3
Fig_Gene_Counts[nrow(Fig_Gene_Counts) + 1,] = c(3,"B", sum(as.factor(as.numeric(Cluster_Nobt$normalized$genes %in% subset(Cluster_Nobt$normalized, cluster==unique(Cluster_Nobt$normalized$cluster)[1])$genes))==1))
Fig_Gene_Counts[nrow(Fig_Gene_Counts) + 1,] = c(3,"C", sum(as.factor(as.numeric(Cluster_Nobt$normalized$genes %in% subset(Cluster_Nobt$normalized, cluster==unique(Cluster_Nobt$normalized$cluster)[2])$genes))==1))
Fig_Gene_Counts[nrow(Fig_Gene_Counts) + 1,] = c(3,"D", sum(as.factor(as.numeric(Cluster_Nobt$normalized$genes %in% subset(Cluster_Nobt$normalized, cluster==unique(Cluster_Nobt$normalized$cluster)[3])$genes))==1))
Fig_Gene_Counts[nrow(Fig_Gene_Counts) + 1,] = c(3,"E", sum(as.factor(as.numeric(Cluster_Nobt$normalized$genes %in% subset(Cluster_Nobt$normalized, cluster==unique(Cluster_Nobt$normalized$cluster)[4])$genes))==1))
Fig_Gene_Counts[nrow(Fig_Gene_Counts) + 1,] = c(3,"F", sum(as.factor(as.numeric(Cluster_Nobt$normalized$genes %in% subset(Cluster_Nobt$normalized, cluster==unique(Cluster_Nobt$normalized$cluster)[5])$genes))==1))
Fig_Gene_Counts[nrow(Fig_Gene_Counts) + 1,] = c(3,"G", sum(as.factor(as.numeric(Cluster_Nobt$normalized$genes %in% subset(Cluster_Nobt$normalized, cluster==unique(Cluster_Nobt$normalized$cluster)[6])$genes))==1))
#Fig 4
Fig_Gene_Counts[nrow(Fig_Gene_Counts) + 1,] = c(4,"B", sum(Cluster_Solanaceae_Noise$normalized$cluster==3))
# Fig 9
Fig_Gene_Counts[nrow(Fig_Gene_Counts) + 1,] = c(9,"A", sum(Cluster_Fruit$normalized$cluster==4))
Fig_Gene_Counts[nrow(Fig_Gene_Counts) + 1,] = c(9,"C", sum(Cluster_Fruit$normalized$cluster==6))
Fig_Gene_Counts[nrow(Fig_Gene_Counts) + 1,] = c(9,"E", sum(Cluster_Fruit$normalized$cluster==11)) # Ungapped cluster=7, but must use original cluster number=11
Fig_Gene_Counts[nrow(Fig_Gene_Counts) + 1,] = c(9,"G", sum(Cluster_Fruit$normalized$cluster==12)) # Ungapped cluster=8, but must use original cluster number=12

write.csv(Fig_Gene_Counts,
          "./Tables/SelectedClusterGeneCounts.csv",row.names = FALSE)

## Cluster Genes and Expression
# In review we were asked to provide a supplemental list of genes in each cluster with their expression values
# at each stage for every cluster in the MS. We'll do it, of course, but I'm not certain why they want this data.
# I also don't know the most useful way to provide this data. Counts would need normalization and aggregation,
# which would be a lot of work. Z-scores are not readily interpretable, but are aggregated properly and let 
# you recreate the figures.

Z_Solanum <- Cluster_Solanum_Noise$normalized %>% 
  dplyr::select(Gene=genes, Zscore=value, Cluster=cluster,  DAP) %>%
  dplyr::filter(Cluster %in% c(4,10)) %>% 
  mutate(Stage = recode(DAP, !!!Labs_SolanumStage),
         DAP = NULL) %>% 
  arrange(Cluster, Stage, Gene)

Z_Nobt <- Cluster_Nobt$normalized %>% 
  dplyr::select(Gene=genes, Zscore=value, Cluster=cluster,  DAP) %>%
  mutate(Stage = recode(DAP, !!!Labs_NobtStage),
         DAP = NULL) %>% 
  arrange(Cluster, Stage, Gene)

# get orthos
orthogroups <- read.table("Orthofinder/OrthoFinder/Results_Oct29/Orthogroups/Orthogroups.tsv",
                          sep="\t",
                          stringsAsFactors = F,
                          header=T) %>% 
  filter_all(all_vars(!grepl(',',.)))

Z_Solanaceae <- Cluster_Solanaceae_Noise$normalized %>% 
  dplyr::select(Orthogroup=genes, Zscore=value, Cluster=cluster, Stage) %>% 
  dplyr::filter(Cluster %in% c(3)) %>% 
  mutate(Stage = recode(Stage, !!!Labs_Stage)) %>% 
  left_join(orthogroups, by = "Orthogroup") %>%
  dplyr::select(Solanum_Gene=Solanum, Nicotiana_Gene=Nicotiana, Cluster, Stage, Zscore) %>% 
  arrange(Cluster, Stage, Solanum_Gene)

# Get orthos
tmp_Orthos <- VennOrthos %>% 
  filter_all(all_vars(!grepl("^$",.))) %>% 
  filter_all(all_vars(!grepl(',',.)))

Z_Model2 <- Cluster_Fruit$normalized %>% 
  dplyr::select(Orthogroup=genes, Zscore=value, Cluster=cluster, Fruit_Type=Fruit, Stage) %>% 
  dplyr::filter(Cluster %in% c(4,6,11,12)) %>% 
  mutate(Stage = recode(Stage, !!!Labs_Stage),
         Cluster = case_when(Cluster == 4 ~ 4,
                             Cluster == 6 ~ 6,
                             Cluster == 11 ~ 7,
                             Cluster == 12 ~ 8)) %>% 
  left_join(tmp_Orthos, by = "Orthogroup") %>% 
  dplyr::select(Arabidopsis_Gene=Arabidopsis, Cucumis_Gene=Cucumis,
                Solanum_Gene=Solanum, Nicotiana_Gene=Nicotiana, Cluster, Stage, Zscore) %>% 
  arrange(Cluster, Stage, Arabidopsis_Gene)

# Export to Excel
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "Figure1") 
openxlsx::writeData(wb,
                    sheet = "Figure1",
                    x=Z_Solanum) 
openxlsx::addWorksheet(wb, "Figure3") 
openxlsx::writeData(wb,
                    sheet = "Figure3",
                    x=Z_Nobt) 
openxlsx::addWorksheet(wb, "Figure4") 
openxlsx::writeData(wb,
                    sheet = "Figure4",
                    x=Z_Solanaceae) 
openxlsx::addWorksheet(wb, "Figure9") 
openxlsx::writeData(wb,
                    sheet = "Figure9",
                    x=Z_Model2) 
openxlsx::saveWorkbook(wb,
                       file = file.path("./Tables/SelectedClusterExpression.xlsx"),
                       overwrite = T)
  









