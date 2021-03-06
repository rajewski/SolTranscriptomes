library("GenomicAlignments")
library("Rsamtools")
library("GenomicFeatures")
library("DESeq2")
library("ggplot2")
library("magrittr")
library("dplyr")
library("tibble")
library("tidyr")
library("UpSetR")
source("X_Functions.R")

#Based on the upsetting upset plot of ripening genes, I wanted to confirm that our data lined up a bit with the Mizzotti data, so I made a corresponding upset plot of each one to compare. This will not be in the manuscript, but is a useful quality check so I don't kill myself.


#Do it for Nobt
Expt_Nobt_All <- readRDS("DEGAnalysis/RNA-seq/Expt_NobtSE.rds") #try with everything as SE

colData(Expt_Nobt_All)$DAP <- as.factor(colData(Expt_Nobt_All)$DAP)
DDS_NobtPairwise <- DESeqDataSet(Expt_Nobt_All,
                                 design= ~ DAP)
DDS_NobtPairwise$Accession <- sub("\\.\\d", "", DDS_NobtPairwise$Accession) #converted to char FYI
DDS_NobtPairwise <- collapseReplicates(DDS_NobtPairwise,
                                         DDS_NobtPairwise$Accession)
DDS_NobtPairwise <- estimateSizeFactors(DDS_NobtPairwise)
DDS_NobtPairwise <- DESeq(DDS_NobtPairwise)

# 0 v 3 DPA
Res_DAP_03 <- results(DDS_NobtPairwise, 
              contrast=c("DAP", "0", "3"), 
              alpha=0.05)
Res_DAP_03 <- lfcShrink(DDS_NobtPairwise, 
                contrast=c("DAP", "0", "3"), 
                res=Res_DAP_03)
ResSig_DAP_03 <- subset(Res_DAP_03, padj<=0.05)
Up03 <- list(as.data.frame(rownames(ResSig_DAP_03)[ResSig_DAP_03$log2FoldChange>=2]))
Down03 <- list(as.data.frame(rownames(ResSig_DAP_03)[ResSig_DAP_03$log2FoldChange<=2]))

# 0 v 6 DPA
Res_DAP_06 <- results(DDS_NobtPairwise, 
                      contrast=c("DAP", "0", "6"), 
                      alpha=0.05)
Res_DAP_06 <- lfcShrink(DDS_NobtPairwise, 
                       contrast=c("DAP", "0", "6"), 
                       res=Res_DAP_06)
ResSig_DAP_06 <- subset(Res_DAP_06, padj<=0.05)
Up06 <- list(as.data.frame(rownames(ResSig_DAP_06)[ResSig_DAP_06$log2FoldChange>=2]))
Down06 <- list(as.data.frame(rownames(ResSig_DAP_06)[ResSig_DAP_06$log2FoldChange<=2]))

# 0 v 11 DPA
Res_DAP_011 <- results(DDS_NobtPairwise, 
                      contrast=c("DAP", "0", "11"), 
                      alpha=0.05)
Res_DAP_011 <- lfcShrink(DDS_NobtPairwise, 
                       contrast=c("DAP", "0", "11"), 
                       res=Res_DAP_011)
ResSig_DAP_011 <- subset(Res_DAP_011, padj<=0.05)
Up011 <- list(as.data.frame(rownames(ResSig_DAP_011)[ResSig_DAP_011$log2FoldChange>=2]))
Down011 <- list(as.data.frame(rownames(ResSig_DAP_011)[ResSig_DAP_011$log2FoldChange<=2]))

# 3 v 6 DPA
Res_DAP_36 <- results(DDS_NobtPairwise, 
                      contrast=c("DAP", "3", "6"), 
                      alpha=0.05)
Res_DAP_36 <- lfcShrink(DDS_NobtPairwise, 
                       contrast=c("DAP", "3", "6"), 
                       res=Res_DAP_36)
ResSig_DAP_36 <- subset(Res_DAP_36, padj<=0.05)
Up36 <- list(as.data.frame(rownames(ResSig_DAP_36)[ResSig_DAP_36$log2FoldChange>=2]))
Down36 <- list(as.data.frame(rownames(ResSig_DAP_36)[ResSig_DAP_36$log2FoldChange<=2]))

# 3 v 11 DPA
Res_DAP_311 <- results(DDS_NobtPairwise, 
                      contrast=c("DAP", "3", "11"), 
                      alpha=0.05)
Res_DAP_311 <- lfcShrink(DDS_NobtPairwise, 
                       contrast=c("DAP", "3", "11"), 
                       res=Res_DAP_311)
ResSig_DAP_311 <- subset(Res_DAP_311, padj<=0.05)
Up311 <- list(as.data.frame(rownames(ResSig_DAP_311)[ResSig_DAP_311$log2FoldChange>=2]))
Down311 <- list(as.data.frame(rownames(ResSig_DAP_311)[ResSig_DAP_311$log2FoldChange<=2]))

# 6 v 11 DPA
Res_DAP_611 <- results(DDS_NobtPairwise,
                       contrast=c("DAP", "6", "11"),
                       alpha=0.05)
Res_DAP_611 <- lfcShrink(DDS_NobtPairwise, 
                       contrast=c("DAP", "6", "11"), 
                       res=Res_DAP_611)
ResSig_DAP_611 <- subset(Res_DAP_611, padj<=0.05)
Up611 <- list(as.data.frame(rownames(ResSig_DAP_611)[ResSig_DAP_611$log2FoldChange>=2]))
Down611 <- list(as.data.frame(rownames(ResSig_DAP_611)[ResSig_DAP_611$log2FoldChange<=2]))

List_Nobt <- Map(list,
                 Up03,Up06,Up011,Up36,Up311,Up611,
                 Down03,Down06,Down011,Down36,Down311,Down611)
Gene_List_Nobt <- as.data.frame(unique(unlist(List_Nobt)))


# Make a table of orthogenes and their presence/absence in a datas --------
for (GeneTable in c(Up03,Up06,Up011,Up36,Up311,Up611,
                    Down03,Down06,Down011,Down36,Down311,Down611)) {
  Common_Nobt<-as.data.frame(Gene_List_Nobt[,1] %in% as.data.frame(GeneTable)[,1])
  Common_Nobt[,1][Common_Nobt[,1]==TRUE]<- 1
  Gene_List_Nobt<-cbind(Gene_List_Nobt,Common_Nobt)
}
# Rename the columns for easier ID later
colnames(Gene_List_Nobt) <- c("ID","Up03","Up06","Up011","Up36","Up311","Up611",
                              'Down03',"Down06","Down011","Down36","Down311","Down611")
colnames(Gene_List_Nobt) <- c("ID","Up 1-2","Up 1-3","Up 1-Br","Up 2-3","Up 2-Br","Up 3-Br",
                              'Down 1-2',"Down 1-3","Down 1-Br","Down 2-3","Down 2-Br","Down 3-Br")

# Make the plot with Upset ------------------------------------------------
upset(Gene_List_Nobt, 
      sets = colnames(Gene_List_Nobt)[-1],
      mainbar.y.label = "Number of genes", 
      sets.x.label = "Total number of genes", 
      text.scale = 1.5,
      point.size = 3, 
      line.size = 1,
      nintersects = 100)


# Arabidopsis -------------------------------------------------------------
Expt_TAIR <- readRDS("DEGAnalysis/RNA-seq/Expt_TAIR.rds")
colData(Expt_TAIR)$DAP <- as.factor(colData(Expt_TAIR)$DAP)
DDS_TAIRPairwise <- DESeqDataSet(Expt_TAIR,
                                 design= ~ DAP)
DDS_TAIRPairwise$Accession <- sub("\\.\\d", "", DDS_TAIRPairwise$Accession) #converted to char FYI
DDS_TAIRPairwise <- collapseReplicates(DDS_TAIRPairwise,
                                       DDS_TAIRPairwise$Accession)
DDS_TAIRPairwise <- estimateSizeFactors(DDS_TAIRPairwise)
DDS_TAIRPairwise <- DESeq(DDS_TAIRPairwise)

# 3 v 6 DPA
TRes_DAP_36 <- results(DDS_TAIRPairwise, 
                      contrast=c("DAP", "3", "6"), 
                      alpha=0.05)
TRes_DAP_36 <- lfcShrink(DDS_TAIRPairwise, 
                        contrast=c("DAP", "3", "6"), 
                        res=TRes_DAP_36)
TResSig_DAP_36 <- subset(TRes_DAP_36, padj<=0.05)
TUp13 <- list(as.data.frame(rownames(TResSig_DAP_36)[TResSig_DAP_36$log2FoldChange>=2]))
TDown13 <- list(as.data.frame(rownames(TResSig_DAP_36)[TResSig_DAP_36$log2FoldChange<=2]))

# 3 v 9 DPA
TRes_DAP_39 <- results(DDS_TAIRPairwise, 
                      contrast=c("DAP", "3", "9"), 
                      alpha=0.05)
TRes_DAP_39 <- lfcShrink(DDS_TAIRPairwise, 
                        contrast=c("DAP", "3", "9"), 
                        res=TRes_DAP_39)
TRes_DAP_39 <- subset(TRes_DAP_39, padj<=0.05)
TUp39 <- list(as.data.frame(rownames(TRes_DAP_39)[TRes_DAP_39$log2FoldChange>=2]))
TDown39 <- list(as.data.frame(rownames(TRes_DAP_39)[TRes_DAP_39$log2FoldChange<=2]))

# 3 v 12 DPA
TRes_DAP_312 <- results(DDS_TAIRPairwise, 
                       contrast=c("DAP", "3", "12"), 
                       alpha=0.05)
TRes_DAP_312 <- lfcShrink(DDS_TAIRPairwise, 
                         contrast=c("DAP", "3", "12"), 
                         res=TRes_DAP_312)
TResSig_DAP_312 <- subset(TRes_DAP_312, padj<=0.05)
TUp312 <- list(as.data.frame(rownames(TResSig_DAP_312)[TResSig_DAP_312$log2FoldChange>=2]))
TDown312 <- list(as.data.frame(rownames(TResSig_DAP_312)[TResSig_DAP_312$log2FoldChange<=2]))

# 6 v 9 DPA
TRes_DAP_69 <- results(DDS_TAIRPairwise, 
                      contrast=c("DAP", "6", "9"), 
                      alpha=0.05)
TRes_DAP_69 <- lfcShrink(DDS_TAIRPairwise, 
                        contrast=c("DAP", "6", "9"), 
                        res=TRes_DAP_69)
TResSig_DAP_69 <- subset(TRes_DAP_69, padj<=0.05)
TUp69 <- list(as.data.frame(rownames(TResSig_DAP_69)[TResSig_DAP_69$log2FoldChange>=2]))
TDown69 <- list(as.data.frame(rownames(TResSig_DAP_69)[TResSig_DAP_69$log2FoldChange<=2]))

# 6 v 12 DPA
TRes_DAP_612 <- results(DDS_TAIRPairwise, 
                       contrast=c("DAP", "6", "12"), 
                       alpha=0.05)
TRes_DAP_612 <- lfcShrink(DDS_TAIRPairwise, 
                         contrast=c("DAP", "6", "12"), 
                         res=TRes_DAP_612)
TResSig_DAP_612 <- subset(TRes_DAP_612, padj<=0.05)
TUp612 <- list(as.data.frame(rownames(TResSig_DAP_612)[TResSig_DAP_612$log2FoldChange>=2]))
TDown612 <- list(as.data.frame(rownames(TResSig_DAP_612)[TResSig_DAP_612$log2FoldChange<=2]))

# 9 v 12 DPA
TRes_DAP_912 <- results(DDS_TAIRPairwise,
                       contrast=c("DAP", "9", "12"),
                       alpha=0.05)
TRes_DAP_912 <- lfcShrink(DDS_TAIRPairwise, 
                         contrast=c("DAP", "9", "12"), 
                         res=TRes_DAP_912)
TResSig_DAP_912 <- subset(TRes_DAP_912, padj<=0.05)
TUp912 <- list(as.data.frame(rownames(TResSig_DAP_912)[TResSig_DAP_912$log2FoldChange>=2]))
TDown912 <- list(as.data.frame(rownames(TResSig_DAP_912)[TResSig_DAP_912$log2FoldChange<=2]))

List_TAIR <- Map(list,
                 TUp36,TUp39,TUp312,TUp69,TUp612,TUp912,
                 TDown36,TDown39,TDown312,TDown69,TDown612,TDown912)
Gene_List_TAIR <- as.data.frame(unique(unlist(List_TAIR)))


# Make a table of orthogenes and their presence/absence in a datas --------
for (GeneTable in c(TUp36,TUp39,TUp312,TUp69,TUp612,TUp912,
                    TDown36,TDown39,TDown312,TDown69,TDown612,TDown912)) {
  Common_TAIR<-as.data.frame(Gene_List_TAIR[,1] %in% as.data.frame(GeneTable)[,1])
  Common_TAIR[,1][Common_TAIR[,1]==TRUE]<- 1
  Gene_List_TAIR<-cbind(Gene_List_TAIR,Common_TAIR)
}
# Rename the columns for easier ID later
colnames(Gene_List_TAIR) <- c("ID","TUp36","TUp39","TUp312","TUp69","TUp612","TUp912",
                              "TDown36","TDown39","TDown312","TDown69","TDown612","TDown912")
colnames(Gene_List_TAIR) <- c("ID","Up 1-2","Up 1-3","Up 1-Br","Up 2-3","Up 2-Br","Up 3-Br",
                              'Down 1-2',"Down 1-3","Down 1-Br","Down 2-3","Down 2-Br","Down 3-Br")


# Make the plot with Upset ------------------------------------------------
upset(Gene_List_TAIR, 
      sets = colnames(Gene_List_TAIR)[-1],
      mainbar.y.label = "Number of genes", 
      sets.x.label = "Total number of genes", 
      text.scale = 1.5,
      point.size = 3, 
      line.size = 1,
      nintersects = 100)


# Make Slyc Data ----------------------------------------------------------
Expt_Slyc <- readRDS("DEGAnalysis/RNA-seq/Expt_SlycSE.rds")
colData(Expt_Slyc)$DAP <- as.factor(colData(Expt_Slyc)$DAP)
DDS_SlycPairwise <- DESeqDataSet(Expt_Slyc,
                                 design= ~ DAP)
DDS_SlycPairwise <- DESeq(subset(DDS_SlycPairwise, select=as.numeric(DAP)<45))

# 1 v 3 DPA
SRes_DAP_13 <- results(DDS_SlycPairwise, 
                       contrast=c("DAP", "1", "3"), 
                       alpha=0.05)
SRes_DAP_13 <- lfcShrink(DDS_SlycPairwise, 
                         contrast=c("DAP", "1", "3"), 
                         res=SRes_DAP_13)
SResSig_DAP_13 <- subset(SRes_DAP_13, padj<=0.05)
SUp13 <- list(as.data.frame(rownames(SResSig_DAP_13)[SResSig_DAP_13$log2FoldChange>=2]))
SDown13 <- list(as.data.frame(rownames(SResSig_DAP_13)[SResSig_DAP_13$log2FoldChange<=2]))

# 0 v 15 DPA
SRes_DAP_115 <- results(DDS_SlycPairwise, 
                       contrast=c("DAP", "1", "15"), 
                       alpha=0.05)
SRes_DAP_115 <- lfcShrink(DDS_SlycPairwise, 
                         contrast=c("DAP", "1", "15"), 
                         res=SRes_DAP_115)
SRes_DAP_115 <- subset(SRes_DAP_115, padj<=0.05)
SUp115 <- list(as.data.frame(rownames(SRes_DAP_115)[SRes_DAP_115$log2FoldChange>=2]))
SDown115 <- list(as.data.frame(rownames(SRes_DAP_115)[SRes_DAP_115$log2FoldChange<=2]))

# 0 v 11 DPA
SRes_DAP_135 <- results(DDS_SlycPairwise, 
                        contrast=c("DAP", "1", "35"), 
                        alpha=0.05)
SRes_DAP_135 <- lfcShrink(DDS_SlycPairwise, 
                          contrast=c("DAP", "1", "35"), 
                          res=SRes_DAP_135)
SResSig_DAP_135 <- subset(SRes_DAP_135, padj<=0.05)
SUp135 <- list(as.data.frame(rownames(SResSig_DAP_135)[SResSig_DAP_135$log2FoldChange>=2]))
SDown135 <- list(as.data.frame(rownames(SResSig_DAP_135)[SResSig_DAP_135$log2FoldChange<=2]))

# 3 v 15 DPA
SRes_DAP_315 <- results(DDS_SlycPairwise, 
                       contrast=c("DAP", "3", "15"), 
                       alpha=0.05)
SRes_DAP_315 <- lfcShrink(DDS_SlycPairwise, 
                         contrast=c("DAP", "3", "15"), 
                         res=SRes_DAP_315)
SResSig_DAP_315 <- subset(SRes_DAP_315, padj<=0.05)
SUp315 <- list(as.data.frame(rownames(SResSig_DAP_315)[SResSig_DAP_315$log2FoldChange>=2]))
SDown315 <- list(as.data.frame(rownames(SResSig_DAP_315)[SResSig_DAP_315$log2FoldChange<=2]))

# 3 v 35 DPA
SRes_DAP_335 <- results(DDS_SlycPairwise, 
                        contrast=c("DAP", "3", "35"), 
                        alpha=0.05)
SRes_DAP_335 <- lfcShrink(DDS_SlycPairwise, 
                          contrast=c("DAP", "3", "35"), 
                          res=SRes_DAP_335)
SResSig_DAP_335 <- subset(SRes_DAP_335, padj<=0.05)
SUp335 <- list(as.data.frame(rownames(SResSig_DAP_335)[SResSig_DAP_335$log2FoldChange>=2]))
SDown335 <- list(as.data.frame(rownames(SResSig_DAP_335)[SResSig_DAP_335$log2FoldChange<=2]))

# 15 v 35 DPA
SRes_DAP_1535 <- results(DDS_SlycPairwise,
                        contrast=c("DAP", "15", "35"),
                        alpha=0.05)
SRes_DAP_1535 <- lfcShrink(DDS_SlycPairwise, 
                          contrast=c("DAP", "15", "35"), 
                          res=SRes_DAP_1535)
SResSig_DAP_1535 <- subset(SRes_DAP_1535, padj<=0.05)
SUp1535 <- list(as.data.frame(rownames(SResSig_DAP_1535)[SResSig_DAP_1535$log2FoldChange>=2]))
SDown1535 <- list(as.data.frame(rownames(SResSig_DAP_1535)[SResSig_DAP_1535$log2FoldChange<=2]))

List_Slyc <- Map(list,
                 SUp13,SUp115,SUp135,SUp315,SUp335,SUp1535,
                 SDown13,SDown115,SDown135,SDown315,SDown335,SDown1535)
Gene_List_Slyc <- as.data.frame(unique(unlist(List_Slyc)))


# Make a table of orthogenes and their presence/absence in a datas --------
for (GeneTable in c(SUp13,SUp115,SUp135,SUp315,SUp335,SUp1535,
                    SDown13,SDown115,SDown135,SDown315,SDown335,SDown1535)) {
  Common_Slyc<-as.data.frame(Gene_List_Slyc[,1] %in% as.data.frame(GeneTable)[,1])
  Common_Slyc[,1][Common_Slyc[,1]==TRUE]<- 1
  Gene_List_Slyc<-cbind(Gene_List_Slyc,Common_Slyc)
}
# Rename the columns for easier ID later
colnames(Gene_List_Slyc) <- c("ID","Up13","Up115","Up135","Up315","Up335","Up1535",
                              'Down13',"Down115","Down135","Down315","Down335","Down1535")
colnames(Gene_List_Slyc) <- c("ID","Up 1-2","Up 1-3","Up 1-Br","Up 2-3","Up 2-Br","Up 3-Br",
                              'Down 1-2',"Down 1-3","Down 1-Br","Down 2-3","Down 2-Br","Down 3-Br")


# Make the plot with Upset ------------------------------------------------
upset(Gene_List_Slyc, 
      sets = colnames(Gene_List_Slyc)[-1],
      mainbar.y.label = "Number of genes", 
      sets.x.label = "Total number of genes", 
      text.scale = 1.5,
      point.size = 3, 
      line.size = 1,
      nintersects = 100)


# Pimp --------------------------------------------------------------------
# Make Slyc Data ----------------------------------------------------------
Expt_Spimp <- readRDS("DEGAnalysis/RNA-seq/Expt_SpimpSE.rds")
colData(Expt_Spimp)$DAP <- as.factor(colData(Expt_Spimp)$DAP)
DDS_SpimpPairwise <- DESeqDataSet(Expt_Spimp,
                                 design= ~ DAP)
DDS_SpimpPairwise <- DESeq(subset(DDS_SpimpPairwise, select=as.numeric(DAP)<45))

# 1 v 3 DPA
SpRes_DAP_13 <- results(DDS_SpimpPairwise, 
                       contrast=c("DAP", "1", "3"), 
                       alpha=0.05)
SpRes_DAP_13 <- lfcShrink(DDS_SpimpPairwise, 
                         contrast=c("DAP", "1", "3"), 
                         res=SpRes_DAP_13)
SpResSig_DAP_13 <- subset(SpRes_DAP_13, padj<=0.05)
SpUp13 <- list(as.data.frame(rownames(SpResSig_DAP_13)[SpResSig_DAP_13$log2FoldChange>=2]))
SpDown13 <- list(as.data.frame(rownames(SpResSig_DAP_13)[SpResSig_DAP_13$log2FoldChange<=2]))

# 0 v 15 DPA
SpRes_DAP_115 <- results(DDS_SpimpPairwise, 
                        contrast=c("DAP", "1", "15"), 
                        alpha=0.05)
SpRes_DAP_115 <- lfcShrink(DDS_SpimpPairwise, 
                          contrast=c("DAP", "1", "15"), 
                          res=SpRes_DAP_115)
SpRes_DAP_115 <- subset(SpRes_DAP_115, padj<=0.05)
SpUp115 <- list(as.data.frame(rownames(SpRes_DAP_115)[SpRes_DAP_115$log2FoldChange>=2]))
SpDown115 <- list(as.data.frame(rownames(SpRes_DAP_115)[SpRes_DAP_115$log2FoldChange<=2]))

# 0 v 11 DPA
SpRes_DAP_135 <- results(DDS_SpimpPairwise, 
                        contrast=c("DAP", "1", "35"), 
                        alpha=0.05)
SpRes_DAP_135 <- lfcShrink(DDS_SpimpPairwise, 
                          contrast=c("DAP", "1", "35"), 
                          res=SpRes_DAP_135)
SpResSig_DAP_135 <- subset(SpRes_DAP_135, padj<=0.05)
SpUp135 <- list(as.data.frame(rownames(SpResSig_DAP_135)[SpResSig_DAP_135$log2FoldChange>=2]))
SpDown135 <- list(as.data.frame(rownames(SpResSig_DAP_135)[SpResSig_DAP_135$log2FoldChange<=2]))

# 3 v 15 DPA
SpRes_DAP_315 <- results(DDS_SpimpPairwise, 
                        contrast=c("DAP", "3", "15"), 
                        alpha=0.05)
SpRes_DAP_315 <- lfcShrink(DDS_SpimpPairwise, 
                          contrast=c("DAP", "3", "15"), 
                          res=SpRes_DAP_315)
SpResSig_DAP_315 <- subset(SpRes_DAP_315, padj<=0.05)
SpUp315 <- list(as.data.frame(rownames(SpResSig_DAP_315)[SpResSig_DAP_315$log2FoldChange>=2]))
SpDown315 <- list(as.data.frame(rownames(SpResSig_DAP_315)[SpResSig_DAP_315$log2FoldChange<=2]))

# 3 v 11 DPA
SpRes_DAP_335 <- results(DDS_SpimpPairwise, 
                        contrast=c("DAP", "3", "35"), 
                        alpha=0.05)
SpRes_DAP_335 <- lfcShrink(DDS_SpimpPairwise, 
                          contrast=c("DAP", "3", "35"), 
                          res=SpRes_DAP_335)
SpResSig_DAP_335 <- subset(SpRes_DAP_335, padj<=0.05)
SpUp335 <- list(as.data.frame(rownames(SpResSig_DAP_335)[SpResSig_DAP_335$log2FoldChange>=2]))
SpDown335 <- list(as.data.frame(rownames(SpResSig_DAP_335)[SpResSig_DAP_335$log2FoldChange<=2]))

# 15 v 11 DPA
SpRes_DAP_1535 <- results(DDS_SpimpPairwise,
                         contrast=c("DAP", "15", "35"),
                         alpha=0.05)
SpRes_DAP_1535 <- lfcShrink(DDS_SpimpPairwise, 
                           contrast=c("DAP", "15", "35"), 
                           res=SpRes_DAP_1535)
SpResSig_DAP_1535 <- subset(SpRes_DAP_1535, padj<=0.05)
SpUp1535 <- list(as.data.frame(rownames(SpResSig_DAP_1535)[SpResSig_DAP_1535$log2FoldChange>=2]))
SpDown1535 <- list(as.data.frame(rownames(SpResSig_DAP_1535)[SpResSig_DAP_1535$log2FoldChange<=2]))

List_Spimp <- Map(list,
                 SpUp13,SpUp115,SpUp135,SpUp315,SpUp335,SpUp1535,
                 SpDown13,SpDown115,SpDown135,SpDown315,SpDown335,SpDown1535)
Gene_List_Spimp <- as.data.frame(unique(unlist(List_Spimp)))


# Make a table of orthogenes and their presence/absence in a datas --------
for (GeneTable in c(SpUp13,SpUp115,SpUp135,SpUp315,SpUp335,SpUp1535,
                    SpDown13,SpDown115,SpDown135,SpDown315,SpDown335,SpDown1535)) {
  Common_Spimp<-as.data.frame(Gene_List_Spimp[,1] %in% as.data.frame(GeneTable)[,1])
  Common_Spimp[,1][Common_Spimp[,1]==TRUE]<- 1
  Gene_List_Spimp<-cbind(Gene_List_Spimp,Common_Spimp)
}
# Rename the columns for easier ID later
colnames(Gene_List_Spimp) <- c("ID","Up13","Up115","Up135","Up315","Up335","Up1535",
                              'Down13',"Down115","Down135","Down315","Down335","Down1535")
colnames(Gene_List_Spimp) <- c("ID","Up 1-2","Up 1-3","Up 1-Br","Up 2-3","Up 2-Br","Up 3-Br",
                              'Down 1-2',"Down 1-3","Down 1-Br","Down 2-3","Down 2-Br","Down 3-Br")

# Make the plot with Upset ------------------------------------------------
upset(Gene_List_Spimp, 
      sets = colnames(Gene_List_Spimp)[-1],
      mainbar.y.label = "Number of genes", 
      sets.x.label = "Total number of genes", 
      text.scale = 1.5,
      point.size = 3, 
      line.size = 1,
      nintersects = 100)



