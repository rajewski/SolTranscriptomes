# Make an upset plot with the Nobt and Slyc Ripening orthogenes
library("DESeq2")
library("UpSetR")
library("readr")

#Read in the data and filter
DDS_NobtRipe_Ortho <- readRDS("DEGAnalysis/RNA-seq/DDS_NobtRipe_Ortho.rds")
DDS_SlycRipe_Ortho <- readRDS("DEGAnalysis/RNA-seq/DDS_SlycRipe_Ortho.rds")

ResultsNobt <- results(DDS_NobtRipe_Ortho, alpha=0.05)
ResultsNobt <- subset(ResultsNobt, padj<=0.05)
UpNobt <- list(as.data.frame(rownames(ResultsNobt)[ResultsNobt$log2FoldChange>2]))
DownNobt <- list(as.data.frame(rownames(ResultsNobt)[ResultsNobt$log2FoldChange<=-2]))

ResultsSlyc <- results(DDS_SlycRipe_Ortho, alpha=0.05)
ResultsSlyc <- subset(ResultsSlyc, padj<=0.05)
UpSlyc <- list(as.data.frame(rownames(ResultsSlyc)[ResultsSlyc$log2FoldChange>=2]))
DownSlyc <- list(as.data.frame(rownames(ResultsSlyc)[ResultsSlyc$log2FoldChange<=-2]))

List_Ripening <- Map(list, UpNobt, DownNobt, UpSlyc, DownSlyc)
Gene_List_Ripening <- as.data.frame(unique(unlist(List_Ripening)))

# Make a table oforthogenes and their presence/absence in a dataset
for (GeneTable in c(UpNobt, DownNobt, UpSlyc, DownSlyc)) {
  Common_Ripening<-as.data.frame(Gene_List_Ripening[,1] %in% as.data.frame(GeneTable)[,1])
  Common_Ripening[,1][Common_Ripening[,1]==TRUE]<- 1
  Gene_List_Ripening<-cbind(Gene_List_Ripening,Common_Ripening)
}
# Rename the columns for easier ID later
colnames(Gene_List_Ripening) <- c("ID", "Tobacco Up", "Tobacco Down", "Tomato Up", "Tomato Down")

# make the plot with Upset
upset(Gene_List_Ripening, 
      sets = colnames(Gene_List_Ripening)[-1],
      mainbar.y.label = "Number of genes", 
      sets.x.label = "Total number of genes", 
      text.scale = 1.5, 
      group.by = "degree", 
      order.by="freq", 
      point.size = 3, 
      line.size = 1,
      nintersects = 15)

