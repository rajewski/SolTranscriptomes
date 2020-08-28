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
source("X_Functions.R")

#set a command line argument in order to run these analyses individually
args <- commandArgs(trailingOnly = TRUE)
if (args =="TAIR"){
  # TAIR 
  TAIR10dds <- readRDS("DEGAnalysis/RNA-seq/DDS_TAIR.rds")
  TAIR10cluster <- DESeqCluster(TAIR10dds, numGenes = "all")
  saveRDS(TAIR10cluster, "DEGAnalysis/RNA-seq/Cluster_TAIR.rds")
  }


if (args=="SlycSRA") {
  # Slyc SRA
  SlycSRAdds <- readRDS("DEGAnalysis/RNA-seq/DDS_SlycSRA.rds")
  SlycSRAcluster <- DESeqCluster(SlycSRAdds, numGenes = "all")
  saveRDS(SlycSRAcluster, "DEGAnalysis/RNA-seq/Cluster_SlycSRA.rds")
}

if (args=="Slyc") {
  # Slyc IH
  SlycIHdds <- readRDS("DEGAnalysis/RNA-seq/DDS_Slyc.rds")
  SlycIHcluster <- DESeqCluster(SlycIHdds, numGenes = "all")
  saveRDS(SlycIHcluster, "DEGAnalysis/RNA-seq/Cluster_Slyc.rds")
}

if (args=="Nobt") {
  # Nobt
  Nobtdds <- readRDS("DEGAnalysis/RNA-seq/DDS_Nobt.rds")
  Nobtcluster <- DESeqCluster(Nobtdds, numGenes = "all")
  saveRDS(Nobtcluster, "DEGAnalysis/RNA-seq/Cluster_Nobt.rds")
  }