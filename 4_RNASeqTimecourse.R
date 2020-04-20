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

# New Functions -----------------------------------------------------------

# Define two new functions to 1) fit a spline model of DEGs and 2) cluster them

DESeqSpline <- function(se=se, 
                        timeVar="DAP",
                        CaseCtlVar="Genotype",
                        SetDF="") {
  #calculate spline df as the number of times -1
  if (is.na(as.numeric(SetDF))){
    dfSpline <- (length(unique(colData(se)[,timeVar]))-1)
  } else if (as.numeric(SetDF)>0) {
    dfSpline <- SetDF
  } else {
    stop("Either leave SetDF blank or provide a positive integer for spline degrees of freedom")
  }
  message("Fitting spline regression with ", dfSpline, " degrees of freedom")
  design <- ns(colData(se)[,timeVar], df=dfSpline)
  colnames(design) <- paste0("spline", seq(1:dim(design)[2]))
  colData(se) <- cbind(colData(se), design)
  if (length(unique(colData(se)[,CaseCtlVar]))>1) {
    message("Two entries detected for ", CaseCtlVar,". Incorporating this into the model")
    dds <- DESeqDataSet(se, design = as.formula(paste0("~", CaseCtlVar, "+" ,paste(paste0(CaseCtlVar,":",colnames(design)), collapse = "+"))))
  } else {
    message("Only one entry detected for ", CaseCtlVar,". Ignoring this variable for model building")
    dds <- DESeqDataSet(se, design = as.formula(paste0("~",paste(colnames(design), collapse = "+"))))
  }
  dds <- estimateSizeFactors(dds)
  if (length(unique(colData(se)[,CaseCtlVar]))>1) {
    dds <- DESeq(dds,test="LRT", reduced = as.formula(paste0("~",paste(colnames(design), collapse = "+"))))
  } else {
    dds <- DESeq(dds,test="LRT", reduced = ~ 1)
  }
  return(dds)
}


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

# Prep Inputs -------------------------------------------------------------
# Borrorwed heavily from https://www.bioconductor.org/help/course-materials/2015/LearnBioconductorFeb2015/B02.1.1_RNASeqLab.html#construct

# Load a gene list by exon for counting (or make and save one)
tryCatch(Slycgenes <- readRDS("DEGAnalysis/RNA-seq/Slycgenes.rds"), error=function(e){
  Slyctxdb <- makeTxDbFromGFF("SlycDNA/ITAG4.0_gene_models.gff", organism="Solanum lycopersicum")
  Slycgenes <- exonsBy(Slyctxdb, by="tx", use.names=TRUE)
  saveRDS(Slycgenes, "DEGAnalysis/RNA-seq/Slycgenes.rds")
})
tryCatch(TAIR10genes <- readRDS("DEGAnalysis/RNA-seq/TAIR10genes.rds"), error=function(e){
  TAIR10txdb <- makeTxDbFromGFF("ExternalData/TAIR10/TAIR10.gff3", organism="Arabidopsis thaliana")
  TAIR10genes <- exonsBy(TAIR10txdb, by="tx", use.names=TRUE)
  saveRDS(TAIR10genes, "DEGAnalysis/RNA-seq/TAIR10genes.rds")
})
tryCatch(Nobtgenes <- readRDS("DEGAnalysis/RNA-seq/Nobtgenes.rds"), error=function(e){
  Nobttxdb <- makeTxDbFromGFF("NobtDNA/NIOBT_r1.0.update.gff", organism="Nicotiana obtusifolia")
  Nobtgenes <- exonsBy(Nobttxdb, by="tx", use.names=TRUE)
  saveRDS(Nobtgenes, "DEGAnalysis/RNA-seq/Nobtgenes.rds")
})

# Read in the Sample list
metadata <- read.table("DEGAnalysis//RNA-seq/metadata.tsv", header=T, sep="")
metadata$Path <- as.character(metadata$Path)
# Get list of BAMs for each expt
NobtBamFiles <- BamFileList(metadata$Path[metadata$Species=="Tobacco"], yieldSize=2000000)
TAIR10BamFiles <- BamFileList(metadata$Path[metadata$Species=="Arabidopsis"], yieldSize=2000000)
SlycSRABamFiles <- BamFileList(metadata$Path[metadata$Species=="Tomato" & metadata$PE==0], yieldSize=2000000)
SlycIHBamFiles <- BamFileList(metadata$Path[metadata$Species=="Tomato" & metadata$PE==1], yieldSize=50000)
SpimpBamFiles <- BamFileList(metadata$Path[metadata$Species=="Pimpinellifolium"], yieldSize=2000000)
seqinfo(TAIR10BamFiles[1]) #check that it worked

# Count Reads -------------------------------------------------------------
#To my knowledge the SRA experiments were not strand-specific
tryCatch(TAIR10Expt <- readRDS("DEGAnalysis/RNA-seq/TAIR10Expt.rds"), error=function(e){
  TAIR10Expt <- summarizeOverlaps(features=TAIR10genes,
                                  reads=TAIR10BamFiles,
                                  mode="Union",
                                  singleEnd=TRUE,
                                  ignore.strand=TRUE)
  colData(TAIR10Expt) <- DataFrame(metadata[metadata$Species=="Arabidopsis",])
  saveRDS(TAIR10Expt, "DEGAnalysis/RNA-seq/TAIR10Expt.rds")
})
#Option to subset to stages 1-3
TAIR10Expt_3stage <- subset(TAIR10Expt, select=DAP<12)

tryCatch(SlycSRAExpt <- readRDS("DEGAnalysis/RNA-seq/SlycSRAExpt.rds"), error=function(e){
  SlycSRAExpt <- summarizeOverlaps(features=Slycgenes,
                                   reads=SlycSRABamFiles,
                                   mode="Union",
                                   singleEnd=TRUE,
                                   ignore.strand=TRUE)
  colData(SlycSRAExpt) <- DataFrame(metadata[metadata$Species=="Tomato" & metadata$PE==0,])
  saveRDS(SlycSRAExpt, "DEGAnalysis/RNA-seq/SlycSRAExpt.rds")
})
tryCatch(SlycIHExpt <- readRDS("DEGAnalysis/RNA-seq/SlycIHExpt.rds"), error=function(e){
  SlycIHPart<-list()
  for (i in 1:length(metadata$Accession[metadata$Species=="Tomato" & metadata$PE==1])) {
    SlycIHPart[[i]] <- summarizeOverlaps(feature=Slycgenes,
                                         reads=SlycIHBamFiles[i],
                                         mode="Union",
                                         singleEnd=FALSE,
                                         ignore.strand=FALSE)
  }
  SlycIHExpt <- do.call(cbind,SlycIHPart)
  colData(SlycIHExpt) <- DataFrame(metadata[metadata$Species=="Tomato" & metadata$PE==1,])
  saveRDS(SlycIHExpt, "DEGAnalysis/RNA-seq/SlycIHExpt.rds")
})
#Option to subset to stages 1-3
SlycIHExpt_3stage <- subset(SlycIHExpt, select=DAP<35)

tryCatch(SpimpExpt <- readRDS("DEGAnalysis/RNA-seq/SpimpExpt.rds"), error=function(e){
  SpimpPart<-list()
  for (i in 1:length(metadata$Accession[metadata$Species=="Pimpinellifolium"])) {
    SpimpPart[[i]] <- summarizeOverlaps(feature=Slycgenes,
                                         reads=SpimpBamFiles[i],
                                         mode="Union",
                                         singleEnd=FALSE,
                                         ignore.strand=FALSE)
  }
  SpimpExpt <- do.call(cbind,SpimpPart)
  colData(SpimpExpt) <- DataFrame(metadata[metadata$Species=="Pimpinellifolium",])
  saveRDS(SpimpExpt, "DEGAnalysis/RNA-seq/SpimpExpt.rds")
})


tryCatch(NobtExpt <- readRDS("DEGAnalysis/RNA-seq/NobtExpt.rds"), error=function(e){
  NobtExpt <- summarizeOverlaps(feature=Nobtgenes,
                                reads=NobtBamFiles,
                                mode="Union",
                                singleEnd=FALSE,
                                ignore.strand=FALSE,
                                fragments=TRUE,
                                BPPARAM=SerialParam())
  colData(NobtExpt) <- DataFrame(metadata[metadata$Species=="Tobacco",])
  saveRDS(NobtExpt, "DEGAnalysis/RNA-seq/NobtExpt.rds")
})


# Design and DE Testing ----------------------------------------------------
# With length(DAP)-1 degrees of freedom 
tryCatch(TAIR10dds <- readRDS("DEGAnalysis/RNA-seq/TAIR10dds.rds"), error=function(e){
  TAIR10dds <- DESeqSpline(TAIR10Expt)
  saveRDS(TAIR10dds, "DEGAnalysis/RNA-seq/TAIR10dds.rds")
})
tryCatch(TAIR10dds_3stage <- readRDS("DEGAnalysis/RNA-seq/TAIR10dds_3stage.rds"), error=function(e){
  TAIR10dds_3stage <- DESeqSpline(TAIR10Expt_3stage)
  saveRDS(TAIR10dds_3stage, "DEGAnalysis/RNA-seq/TAIR10dds_3stage.rds")
})
tryCatch(SlycSRAdds <- readRDS("DEGAnalysis/RNA-seq/SlycSRAdds.rds"), error=function(e){
  SlycSRAdds <- DESeqSpline(SlycSRAExpt)
  saveRDS(SlycSRAdds, "DEGAnalysis/RNA-seq/SlycSRAdds.rds")
})
tryCatch(SlycIHdds <- readRDS("DEGAnalysis/RNA-seq/SlycIHdds.rds"), error=function(e){
  SlycIHdds <- DESeqSpline(SlycIHExpt)
  saveRDS(SlycIHdds, "DEGAnalysis/RNA-seq/SlycIHdds.rds")
})
tryCatch(SlycIHdds_3stage <- readRDS("DEGAnalysis/RNA-seq/SlycIHdds_3stage.rds"), error=function(e){
  SlycIHdds_3stage <- DESeqSpline(SlycIHExpt_3stage)
  saveRDS(SlycIHdds_3stage, "DEGAnalysis/RNA-seq/SlycIHdds_3stage.rds")
})
tryCatch(Nobtdds <- readRDS("DEGAnalysis/RNA-seq/Nobtdds.rds"), error=function(e){
  Nobtdds <- DESeqSpline(NobtExpt)
  saveRDS(Nobtdds, "DEGAnalysis/RNA-seq/Nobtdds.rds")
})




# Play around with individual genes ---------------------------------------
Exampledds <- SlycIHdds_3stage #assign one dds as the example to streamline code
ExampleRes <- results(Exampledds) #get results
ExampleResSig <- subset(ExampleRes, padj < 0.05) #subset by FDR
head(ExampleResSig[order(ExampleResSig$padj ), ]) #see best fitting genes for spline model
#Examine an individual Gene
topGene <- rownames(ExampleRes)[which.min(ExampleRes$padj)]
colData(Exampledds)$DAP <- as.factor(colData(Exampledds)$DAP)
plotCounts(Exampledds, gene=topGene, intgroup=c("Genotype", "DAP"), normalized = T) #plot best fitting gene

# Get a set of FUL genes for each species. Only use one of these
FULgenes<-c(FUL.1="AT5G60910.1",
            FUL.2="AT5G60910.2",
            AGL79="AT3G30260.1")
FULgenes<-c(SlFUL1="Solyc06g069430.3.1",
            SlFUL2="Solyc03g114830.3.1",
            SlMBP10="Solyc02g065730.2.1",
            SlMBP20="Solyc02g089210.4.1" )
FULgenes<-c(NoFUL1="NIOBTv3_g28929-D2.t1",
            NoFUL2="NIOBTv3_g39464.t1",
            NoMBP10="NIOBTv3_g07845.t1",
            NoMBP20="NIOBT_gMBP20.t1" )
for (i in 1:length(FULgenes)) {
  pdf(file=paste0("DEGAnalysis/RNA-seq/Plot_IH_3Stage_", names(FULgenes[i]), ".pdf"), # _SRA v _IH on Slyc
      width=6,
      height=4)
  plotCounts(Exampledds,
             gene=FULgenes[i],
             intgroup=c("DAP"), # Remove Genotype for nonSlycSRA
             main=paste0(names(FULgenes[i]), 
                         " p=",
                         formatC(ExampleRes$padj[rownames(ExampleRes)==FULgenes[i]], format = "e", digits = 1)),
             xlab="Days Post Anthesis",
             normalized=T,
             replace=T)
  dev.off()
}
 

# Clustering --------------------------------------------------------------
# For all genes, this is best run noninteractively with the 4_NoninteractiveClustering.sh script
tryCatch(TAIR10Allcluster <- readRDS("DEGAnalysis/RNA-seq/TAIR10Allcluster.rds"), error=function(e){
  TAIR10Allcluster <- DESeqCluster(TAIR10dds, numGenes = "all")
  saveRDS(TAIR10Allcluster, "DEGAnalysis/RNA-seq/TAIR10Allcluster.rds")
})
tryCatch(TAIR103Stagecluster <- readRDS("DEGAnalysis/RNA-seq/TAIR103Stagecluster.rds"), error=function(e){
  TAIR10Allcluster <- DESeqCluster(TAIR10dds_3stage, numGenes = "all")
  saveRDS(TAIR103Stagecluster, "DEGAnalysis/RNA-seq/TAIR103Stagecluster.rds")
})
tryCatch(SlycSRAAllcluster <- readRDS("DEGAnalysis/RNA-seq/SlycSRAAllcluster.rds"), error=function(e){
  SlycSRAAllcluster <- DESeqCluster(SlycSRAdds, numGenes = "all")
  saveRDS(SlycSRAAllcluster, "DEGAnalysis/RNA-seq/SlycSRAAllcluster.rds")
})
tryCatch(SlycIHAllcluster <- readRDS("DEGAnalysis/RNA-seq/SlycIHAllcluster.rds"), error=function(e){
  SlycIHAllcluster <- DESeqCluster(SlycIHdds, numGenes = "all")
  saveRDS(SlycIHAllcluster, "DEGAnalysis/RNA-seq/SlycIHAllcluster.rds")
})

tryCatch(SlycIH3Stagecluster <- readRDS("DEGAnalysis/RNA-seq/SlycIH3Stagecluster.rds"), error=function(e){
  SlycIH3Stagecluster <- DESeqCluster(SlycIHdds_3stage, numGenes = "all")
  saveRDS(SlycIH3Stagecluster, "DEGAnalysis/RNA-seq/SlycIH3Stagecluster.rds")
})

tryCatch(NobtAllcluster <- readRDS("DEGAnalysis/RNA-seq/NobtAllcluster.rds"), error=function(e){
  NobtAllcluster <- DESeqCluster(Nobtdds, numGenes = "all")
  saveRDS(NobtAllcluster, "DEGAnalysis/RNA-seq/NobtAllcluster.rds")
})

# Plot Cluster Profiles ---------------------------------------------------
ClusterforPlotting <- TAIR103Stagecluster
PlotCluster <-degPlotCluster(ClusterforPlotting$normalized,
                             time="DAP",
                             boxes=T,
                             points=F,
                             #color="Genotype",
                             #lines=F
                             )
PlotCluster + theme_minimal() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave(filename = "DEGAnalysis/RNA-seq/Plot_TAIR_3Stage_ClusterProfiles.pdf",
       width=11,
       height=7)

# NoFUL2 and NoMBP10 are in cluster 2
# In SlycIH FUL1 is cluster 3
# in SlycSRA FUL1 and FUL2 are cluster 1
# FUL and AGL79 are not DE

# Save Cluster genes to a file --------------------------------------------
X <- split(NobtAllcluster$df, NobtAllcluster$df$cluster)
for (i in 1:length(X)) {
  write.table(row.names(X[[i]]), 
              file=paste0("DEGAnalysis/RNA-seq/Nobt_Cluster_", max(X[[i]]$cluster), ".txt"),
              row.names = FALSE,
              quote = FALSE,
              col.names = FALSE)
}
