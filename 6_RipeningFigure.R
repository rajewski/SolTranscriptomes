# Make an upset plot with the Nobt and Slyc Ripening orthogenes
library("DESeq2")
library("UpSetR")
library("readr")
source("X_Functions.R")
# Read in the data and filter ---------------------------------------------
DDS_NobtRipe_Ortho <- readRDS("DEGAnalysis/RNA-seq/DDS_NobtRipe_Ortho.rds")
DDS_Slyc1545_Ortho <- readRDS("DEGAnalysis/RNA-seq/DDS_Slyc1545_Ortho.rds")
DDS_Slyc1535_Ortho <- readRDS("DEGAnalysis/RNA-seq/DDS_Slyc1535_Ortho.rds")
# DDS_Slyc3545_Ortho <- readRDS("DEGAnalysis/RNA-seq/DDS_Slyc3545_Ortho.rds")

ResultsNobt <- results(DDS_NobtRipe_Ortho, alpha=0.05)
ResultsNobt <- subset(ResultsNobt, padj<=0.05)
UpNobt <- list(as.data.frame(rownames(ResultsNobt)[ResultsNobt$log2FoldChange>2]))
DownNobt <- list(as.data.frame(rownames(ResultsNobt)[ResultsNobt$log2FoldChange<=-2]))

ResultsSlyc1535 <- results(DDS_Slyc1535_Ortho, alpha=0.05)
ResultsSlyc1535 <- subset(ResultsSlyc1535, padj<=0.05)
UpSlyc1535 <- list(as.data.frame(rownames(ResultsSlyc1535)[ResultsSlyc1535$log2FoldChange>=2]))
DownSlyc1535 <- list(as.data.frame(rownames(ResultsSlyc1535)[ResultsSlyc1535$log2FoldChange<=-2]))

ResultsSlyc1545 <- results(DDS_Slyc1545_Ortho, alpha=0.05)
ResultsSlyc1545 <- subset(ResultsSlyc1545, padj<=0.05)
UpSlyc1545 <- list(as.data.frame(rownames(ResultsSlyc1545)[ResultsSlyc1545$log2FoldChange>=2]))
DownSlyc1545 <- list(as.data.frame(rownames(ResultsSlyc1545)[ResultsSlyc1545$log2FoldChange<=-2]))

# ResultsSlyc3545 <- results(DDS_Slyc3545_Ortho, alpha=0.05)
# ResultsSlyc3545 <- subset(ResultsSlyc3545, padj<=0.05)
# UpSlyc3545 <- list(as.data.frame(rownames(ResultsSlyc3545)[ResultsSlyc3545$log2FoldChange>=2]))
# DownSlyc3545 <- list(as.data.frame(rownames(ResultsSlyc3545)[ResultsSlyc3545$log2FoldChange<=-2]))

List_Ripening <- Map(list,
                     UpNobt,
                     DownNobt,
                     UpSlyc1535,
                     DownSlyc1535,
                     UpSlyc1545,
                     DownSlyc1545)
Gene_List_Ripening <- as.data.frame(unique(unlist(List_Ripening)))


# Make a table of orthogenes and their presence/absence in a datas --------
for (GeneTable in c(UpNobt, DownNobt, UpSlyc1535, DownSlyc1535, UpSlyc1545, DownSlyc1545)) {
  Common_Ripening<-as.data.frame(Gene_List_Ripening[,1] %in% as.data.frame(GeneTable)[,1])
  Common_Ripening[,1][Common_Ripening[,1]==TRUE]<- 1
  Gene_List_Ripening<-cbind(Gene_List_Ripening,Common_Ripening)
}
# Rename the columns for easier ID later
colnames(Gene_List_Ripening) <- c("ID", "Tobacco Up", "Tobacco Down", "Tomato 3-Breaker Up", "Tomato 3-Breaker Down", "Tomato 3-Ripe Up", "Tomato 3-Ripe Down")


# Make the plot with Upset ------------------------------------------------
upset(Gene_List_Ripening, 
      sets = colnames(Gene_List_Ripening)[-1],
      mainbar.y.label = "Number of genes", 
      sets.x.label = "Total number of genes", 
      text.scale = 1.5,
      point.size = 3, 
      line.size = 1,
      nintersects = 100)

# Create named factor of presence absence for orthogenes and some combinations
Gene_List_Intersect <- data.frame(ID=Gene_List_Ripening[,1])
# One-Set
Gene_List_Intersect$`Tobacco Up` <- Gene_List_Ripening[,2]*!apply(Gene_List_Ripening[,-c(1,2)],1,FUN = max)
Gene_List_Intersect$`Tobacco Down` <- Gene_List_Ripening[,3]*!apply(Gene_List_Ripening[,-c(1,3)],1,FUN = max)
Gene_List_Intersect$`Tomato 3-Ripe Up` <- Gene_List_Ripening[,6]*!apply(Gene_List_Ripening[,-c(1,6)],1,FUN = max)
Gene_List_Intersect$`Tomato 3-Ripe Down` <- Gene_List_Ripening[,7]*!apply(Gene_List_Ripening[,-c(1,7)],1,FUN = max)
Gene_List_Intersect$`Tomato 3-Breaker Up` <- Gene_List_Ripening[,4]*!apply(Gene_List_Ripening[,-c(1,4)],1,FUN = max)
Gene_List_Intersect$`Tomato 3-Breaker Down` <- Gene_List_Ripening[,5]*!apply(Gene_List_Ripening[,-c(1,5)],1,FUN = max)
# Two-Set
Gene_List_Intersect$`Tobacco Down v Tomato 3-Ripe Up` <- Gene_List_Ripening[,3]*Gene_List_Ripening[,6]*!apply(Gene_List_Ripening[,-c(1,3,6)],1,FUN = max)
Gene_List_Intersect$`Tobacco Down v Tomato 3-Ripe Down` <- Gene_List_Ripening[,3]*Gene_List_Ripening[,7]*!apply(Gene_List_Ripening[,-c(1,3,7)],1,FUN = max)
Gene_List_Intersect$`Tobacco Down v Tomato 3-Breaker Up` <- Gene_List_Ripening[,3]*Gene_List_Ripening[,4]*!apply(Gene_List_Ripening[,-c(1,3,4)],1,FUN = max)
Gene_List_Intersect$`Tobacco Down v Tomato 3-Breaker Down` <- Gene_List_Ripening[,3]*Gene_List_Ripening[,5]*!apply(Gene_List_Ripening[,-c(1,3,5)],1,FUN = max)
Gene_List_Intersect$`Tobacco Up v Tomato 3-Ripe Up` <- Gene_List_Ripening[,2]*Gene_List_Ripening[,6]*!apply(Gene_List_Ripening[,-c(1,2,6)],1,FUN = max)
Gene_List_Intersect$`Tobacco Up v Tomato 3-Ripe Down` <- Gene_List_Ripening[,2]*Gene_List_Ripening[,7]*!apply(Gene_List_Ripening[,-c(1,2,7)],1,FUN = max)
Gene_List_Intersect$`Tobacco Up v Tomato 3-Breaker Down` <- Gene_List_Ripening[,2]*Gene_List_Ripening[,5]*!apply(Gene_List_Ripening[,-c(1,2,5)],1,FUN = max)
Gene_List_Intersect$`Tobacco Up v Tomato 3-Breaker Up` <- Gene_List_Ripening[,2]*Gene_List_Ripening[,4]*!apply(Gene_List_Ripening[,-c(1,2,4)],1,FUN = max)
# Three Set
Gene_List_Intersect$`Tobacco Down v Tomato 3-Breaker/Ripe Up` <- Gene_List_Ripening[,3]*Gene_List_Ripening[,4]*Gene_List_Ripening[,6]*!apply(Gene_List_Ripening[,-c(1,3,4,6)],1,FUN = max)
Gene_List_Intersect$`Tobacco Down v Tomato 3-Breaker/Ripe Down` <- Gene_List_Ripening[,3]*Gene_List_Ripening[,5]*Gene_List_Ripening[,7]*!apply(Gene_List_Ripening[,-c(1,3,5,7)],1,FUN = max)
Gene_List_Intersect$`Tobacco Up v Tomato 3-Breaker/Ripe Up` <- Gene_List_Ripening[,2]*Gene_List_Ripening[,4]*Gene_List_Ripening[,6]*!apply(Gene_List_Ripening[,-c(1,2,4,6)],1,FUN = max)
Gene_List_Intersect$`Tobacco Up v Tomato 3-Breaker/Ripe Down` <- Gene_List_Ripening[,2]*Gene_List_Ripening[,5]*Gene_List_Ripening[,7]*!apply(Gene_List_Ripening[,-c(1,2,5,7)],1,FUN = max)

# Make a list of the GO Tables
GO_Tables <- list()
for (i in 2:dim(Gene_List_Intersect)[2]) {
  tmpList <- as.factor(Gene_List_Intersect[,i])
  names(tmpList) <- Gene_List_Intersect$ID
  GO_Tables[[i-1]] <- GOEnrich(gene2go = "DEGAnalysis/Pfam/Ortho.gene2go.tsv",
                             GOIs=tmpList)
}
names(GO_Tables) <- names(Gene_List_Intersect)[-1]
capture.output(GO_Tables, file="DEGAnalysis/RNA-seq/RipeningGOTables.txt")

# Find some key genes
impt_genes <- c("CNR","NOR","RIN","MC","J","FUL1","FUL2","MBP10","MBP20","NR","TAGL1", "ACS2", "ACS4", "ACO1", "ACO2", "ACO3", "ACO4", "ACO5", "ACO6", "ACO7", "TAG1", "AGL11", "MBP3", "PSY1", "PL1", "Cel2", "PGA2A", "EXP1", "GAD1", "GAD2", "GAD3", "CHS-1", "CHS-2", "GRAS38", "CTOMT1", "TOMLOXC" )
names(impt_genes) <- c("Solyc02g077920.4.1","Solyc10g006880.3.1","Solyc05g012020.4.1","Solyc05g056620.2.1","Solyc11g010570.2.1","Solyc06g069430.3.1","Solyc03g114830.3.1","Solyc02g065730.2.1","Solyc02g089210.4.1","Solyc09g075440.4.1","Solyc07g055920.4.1", "Solyc01g095080.3.1", "Solyc05g050010.3.1", "Solyc07g049530.3.1", "Solyc12g005940.2.1", "Solyc07g049550.3.1", "Solyc02g081190.4.1", "Solyc07g026650.3.1", "Solyc02g036350.3.1", "Solyc06g060070.3.1", "Solyc02g071730.4.1", "Solyc11g028020.3.1", "Solyc06g064840.4.1" ,"Solyc03g031860.3.1", "Solyc03g111690.4.1", "Solyc09g010210.3.1", "Solyc10g080210.2.1", "Solyc06g051800.3.1", "Solyc03g098240.3.1", "Solyc11g011920.2.1", "Solyc01g005000.3.1", "Solyc09g091510.3.1", "Solyc05g053550.3.1", "Solyc07g052960.3.1", "Solyc10g005060.4.1", "Solyc01g006540.4.1")
orthomapping <- read.table("Orthofinder/OrthoFinder/Results_May17/Orthogroups/Orthogroups.tsv", stringsAsFactors = F, sep="\t", header=T)[,c(1,4)]
orthomapping <- orthomapping[orthomapping$Solanum %in% names(impt_genes),]
orthomapping$Abbr <- impt_genes[names(impt_genes) %in% orthomapping$Solanum]
# Find presence of orthogroups of important genes in sets
Gene_List_Ripening[Gene_List_Ripening$ID %in% orthomapping$Orthogroup,]
