library("GenomicAlignments")
library("Rsamtools")
library("GenomicFeatures")
library("splines")
library("DESeq2")
library("ggplot2")
library("magrittr")
library("DEGreport")
library("dplyr")
library("tibble")

DESeqCluster <- function(dds=dds,
                         numGenes=c("1000", "3000", "all"),
                         FDRthreshold=0.01,
                         timeVar="DAP",
                         CaseCtlVar="Genotype"){
  #this section is HEAVILY borrowed from:
  #https://hbctraining.github.io/DGE_workshop/lessons/08_DGE_LRT.html
  numGenes=match.arg(numGenes)
  message("Normalizing counts")
  rld <- rlog(dds)
  message("Done.")
  message("Now for some diagnostic plots:")
  plot( assay(rld)[ , 1:2], col=rgb(0,0,0,.2), pch=16, cex=0.3, )
  plot <- plotPCA(rld, intgroup = c(CaseCtlVar, timeVar))
  print(plot)
  message("Filtering DEGs to FDR threshold less than ", FDRthreshold)
  res_LRT <- results(dds)
  sig_res <- res_LRT %>%
    data.frame() %>%
    rownames_to_column(var="gene") %>% 
    as_tibble() %>% 
    filter(padj < FDRthreshold)
  message(nrow(sig_res), " genes remaining")
  # Subset results for faster clustering
  if (numGenes=="all") {
    numGenes <- nrow(sig_res)
    message("Using all ", numGenes, " DEGs. Maybe consider fewer to speed clustering?")
  } else {
    message("Using ", numGenes, " for clustering")
  }
  clustering_sig_genes <- sig_res %>%
    arrange(padj) %>%
    head(n=as.numeric(numGenes))
  # Obtain rlog values for those significant genes
  cluster_rlog <- rld[clustering_sig_genes$gene, ]
  colData(cluster_rlog)[,timeVar] <- as.factor(colData(cluster_rlog)[,timeVar])
  if (length(unique(colData(cluster_rlog)[,CaseCtlVar]))>1) {
    message("Two entries detected for ", CaseCtlVar,". Plots will be colored by ", CaseCtlVar)
    clusters <- degPatterns(assay(cluster_rlog), metadata = colData(cluster_rlog), time = timeVar, col=CaseCtlVar, reduce = T)
  } else {
    message("Only one entry detected for ", CaseCtlVar, ". Ignoring ", CaseCtlVar, " for plotting.")
    clusters <- degPatterns(assay(cluster_rlog), metadata = colData(cluster_rlog), time = timeVar, col=NULL, reduce = T)
  }
  return(clusters)
}

#set a command line argument in order to run these analyses individually
args <- commandArgs(trailingOnly = TRUE)
if (args =="TAIR"){
  # TAIR 
  TAIR10dds <- readRDS("DEGAnalysis/RNA-seq/TAIR10dds.rds")
  TAIR10cluster <- DESeqCluster(TAIR10dds, numGenes = "all")
  saveRDS(TAIR10cluster, "DEGAnalysis/RNA-seq/TAIR10Allcluster.rds")
  
  X <- split(TAIR10cluster$df, TAIR10cluster$df$cluster)
  for (i in 1:length(X)) {
    write.table(row.names(X[[i]]), 
                file=paste0("DEGAnalysis/RNA-seq/TAIR10_Cluster_", i, ".txt"),
                row.names = FALSE,
                quote = FALSE,
                col.names = FALSE)
  }
}

if (args=="SlycSRA") {
  # Slyc SRA
  SlycSRAdds <- readRDS("DEGAnalysis/RNA-seq/SlycSRAdds.rds")
  SlycSRAcluster <- DESeqCluster(SlycSRAdds, numGenes = "all")
  saveRDS(SlycSRAcluster, "DEGAnalysis/RNA-seq/SlycSRAAllcluster.rds")
  
  X <- split(SlycSRAcluster$df, SlycSRAcluster$df$cluster)
  for (i in 1:length(X)) {
    write.table(row.names(X[[i]]), 
                file=paste0("DEGAnalysis/RNA-seq/SlycSRA_All_Cluster_", i, ".txt"),
                row.names = FALSE,
                quote = FALSE,
                col.names = FALSE)
  }
}

if (args=="SlycIH") {
  # Slyc IH
  SlycIHdds <- readRDS("DEGAnalysis/RNA-seq/SlycIHdds.rds")
  SlycIHcluster <- DESeqCluster(SlycIHdds, numGenes = "all")
  saveRDS(SlycIHcluster, "DEGAnalysis/RNA-seq/SlycIHAllcluster.rds")
  
  X <- split(SlycIHcluster$df, SlycIHcluster$df$cluster)
  for (i in 1:length(X)) {
    write.table(row.names(X[[i]]), 
                file=paste0("DEGAnalysis/RNA-seq/SlycIH_All_Cluster_", i, ".txt"),
                row.names = FALSE,
                quote = FALSE,
                col.names = FALSE)
  }
}

if (args=="SlycIH3Stage") {
  # Slyc IH
  SlycIHdds <- readRDS("DEGAnalysis/RNA-seq/SlycIHdds.rds")
  SlycIHcluster <- DESeqCluster(SlycIHdds, numGenes = "all")
  saveRDS(SlycIHcluster, "DEGAnalysis/RNA-seq/SlycIHAllcluster.rds")
  
  X <- split(SlycIHcluster$df, SlycIHcluster$df$cluster)
  for (i in 1:length(X)) {
    write.table(row.names(X[[i]]), 
                file=paste0("DEGAnalysis/RNA-seq/SlycIH_All_Cluster_", i, ".txt"),
                row.names = FALSE,
                quote = FALSE,
                col.names = FALSE)
  }
}

if (args=="Nobt") {
  # Nobt
  Nobtdds <- readRDS("DEGAnalysis/RNA-seq/Nobtdds.rds")
  Nobtcluster <- DESeqCluster(Nobtdds, numGenes = "all")
  saveRDS(Nobtcluster, "DEGAnalysis/RNA-seq/NobtAllcluster.rds")
  
  X <- split(Nobtcluster$df, Nobtcluster$df$cluster)
  for (i in 1:length(X)) {
    write.table(row.names(X[[i]]), 
                file=paste0("DEGAnalysis/RNA-seq/Nobt_All_Cluster_", i, ".txt"),
                row.names = FALSE,
                quote = FALSE,
                col.names = FALSE)
  }
}