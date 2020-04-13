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
                        CaseCtlVar="Genotype") {
  #calculate spline df as the number of times -1
  dfSpline <- (length(unique(colData(se)[,timeVar]))-1)
  message("Fitting spline regresion with ", dfSpline, " degrees of freedom")
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
metadata <- read.table("DEGAnalysis//RNA-seq/SampleList.txt", header=T, sep="\t")
metadata$Path <- NA # Add Paths to BAM file
metadata$Path[metadata$Species=="Arabidopsis"] <- paste0("DEGAnalysis/STAR/TAIR10/",metadata$Accession[metadata$Species=="Arabidopsis"],".Aligned.sortedByCoord.out.bam")
metadata$Path[metadata$Species=="Tomato"] <- paste0("DEGAnalysis/STAR/Slyc/",metadata$Accession[metadata$Species=="Tomato"],".Aligned.sortedByCoord.out.bam")
metadata$Path[metadata$Species=="Tobacco"] <- paste0("DEGAnalysis/STAR/Nobt/",metadata$Accession[metadata$Species=="Tobacco"],".Aligned.sortedByCoord.out.bam")
# Get list of BAMs for each expt
NobtBamFiles <- BamFileList(metadata$Path[metadata$Species=="Tobacco"], yieldSize=2000000)
TAIR10BamFiles <- BamFileList(metadata$Path[metadata$Species=="Arabidopsis"], yieldSize=2000000)
SlycSRABamFiles <- BamFileList(metadata$Path[metadata$Species=="Tomato" & metadata$PE==0], yieldSize=2000000)
SlycIHBamFiles <- BamFileList(metadata$Path[metadata$Species=="Tomato" & metadata$PE==1], yieldSize=50000) #lower yield size?
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
# I wrote a function to do what is currently my default analysis uniformly on every expt

# Load (or make and save) DESeq datasets for each experiment
tryCatch(TAIR10dds <- readRDS("DEGAnalysis/RNA-seq/TAIR10dds.rds"), error=function(e){
  TAIR10dds <- DESeqSpline(TAIR10Expt)
  saveRDS(TAIR10dds, "DEGAnalysis/RNA-seq/TAIR10dds.rds")
})
tryCatch(SlycSRAdds <- readRDS("DEGAnalysis/RNA-seq/SlycSRAdds.rds"), error=function(e){
  SlycSRAdds <- DESeqSpline(SlycSRAExpt)
  saveRDS(SlycSRAdds, "DEGAnalysis/RNA-seq/SlycSRAdds.rds")
})
tryCatch(SlycIHdds <- readRDS("DEGAnalysis/RNA-seq/SlycIHdds.rds"), error=function(e){
  SlycIHdds <- DESeqSpline(SlycIHExpt)
  saveRDS(SlycIHdds, "DEGAnalysis/RNA-seq/SlycIHdds.rds")
})
tryCatch(Nobtdds <- readRDS("DEGAnalysis/RNA-seq/Nobtdds.rds"), error=function(e){
  Nobtdds <- DESeqSpline(NobtExpt)
  saveRDS(Nobtdds, "DEGAnalysis/RNA-seq/Nobtdds.rds")
})

# Play around with individual genes ---------------------------------------
Exampledds <- SlycIHdds #assign one dds as the example to streamline code
ExampleRes <- results(Exampledds) #get results
ExampleResSig <- subset(ExampleRes, padj < 0.05) #subset by FDR
head(ExampleResSig[order(ExampleResSig$padj ), ]) #see best fitting genes for spline model
#Examine an individual Gene
topGene <- rownames(ExampleRes)[which.min(ExampleRes$padj)]
colData(Exampledds)$DAP <- as.factor(colData(Exampledds)$DAP)
plotCounts(Exampledds, gene=topGene, intgroup="DAP", normalized = T) #plot best fitting gene

#FUL AT5G60910 
#FUL1 Solyc06g069430.3 NIOBTv3_g28929-D2
#FUL1 Solyc03g114830.3 NIOBTv3_g39464
#MBP10 Solyc02g065730.2 NIOBTv307845
#MBP20 Solyc02g089210.4 NIOBT_gMBP20
plotCounts(Exampledds, gene="NIOBTv3_g39464", intgroup="DAP",normalized=T) #FRUITFULL

# Clustering --------------------------------------------------------------

tryCatch(TAIR10cluster <- readRDS("DEGAnalysis/RNA-seq/TAIR10cluster.rds"), error=function(e){
  TAIR10cluster <- DESeqCluster(TAIR10dds, numGenes = "3000")
  saveRDS(TAIR10cluster, "DEGAnalysis/RNA-seq/TAIR10cluster.rds")
})
tryCatch(SlycSRAcluster <- readRDS("DEGAnalysis/RNA-seq/SlycSRAcluster.rds"), error=function(e){
  SlycSRAcluster <- DESeqCluster(SlycSRAdds, numGenes = "3000")
  saveRDS(SlycSRAcluster, "DEGAnalysis/RNA-seq/SlycSRAcluster.rds")
})
tryCatch(SlycIHcluster <- readRDS("DEGAnalysis/RNA-seq/SlycIHcluster.rds"), error=function(e){
  SlycIHcluster <- DESeqCluster(SlycIHdds, numGenes = "3000")
  saveRDS(SlycIHcluster, "DEGAnalysis/RNA-seq/SlycIHcluster.rds")
})

tryCatch(Nobtcluster <- readRDS("DEGAnalysis/RNA-seq/Nobtcluster.rds"), error=function(e){
  Nobtcluster <- DESeqCluster(Nobtdds, numGenes = "3000")
  saveRDS(Nobtcluster, "DEGAnalysis/RNA-seq/Nobtcluster.rds")
})


# Save Cluster genes to a file --------------------------------------------
X <- split(Nobtcluster$df, Nobtcluster$df$cluster)
for (i in 1:length(X)) {
  write.table(row.names(X[[i]]), 
              file=paste0("DEGAnalysis/RNA-seq/Nobt_Cluster_", i, ".txt"),
              row.names = FALSE,
              quote = FALSE,
              col.names = FALSE)
}
