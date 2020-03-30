library("GenomicAlignments")
library("Rsamtools")
library("GenomicFeatures")

# Prep Inputs -------------------------------------------------------------
#Borrorwed heavily from https://www.bioconductor.org/help/course-materials/2015/LearnBioconductorFeb2015/B02.1.1_RNASeqLab.html#construct
#Read in the Sample list
metadata <- read.table("DEGAnalysis/SampleList.txt", header=T, sep="\t")
# Add Path to BAM file
metadata$Path <- NA
metadata$Path[metadata$Species=="Arabidopsis"] <- paste0("DEGAnalysis/STAR/TAIR10/",metadata$Accession[metadata$Species=="Arabidopsis"],".Aligned.sortedByCoord.out.bam")
metadata$Path[metadata$Species=="Tomato"] <- paste0("DEGAnalysis/STAR/Slyc/",metadata$Accession[metadata$Species=="Tomato"],".Aligned.sortedByCoord.out.bam")
metadata$Path[metadata$Species=="Tobacco"] <- paste0("DEGAnalysis/STAR/Nobt/",metadata$Accession[metadata$Species=="Tobacco"],".Aligned.sortedByCoord.out.bam")

#Read in BAM files
NobtBamFiles <- BamFileList(metadata$Path[metadata$Species=="Tobacco"], yieldSize=2000000)
TAIR10BamFiles <- BamFileList(metadata$Path[metadata$Species=="Arabidopsis"], yieldSize=2000000)
SlycSRABamFiles <- BamFileList(metadata$Path[metadata$Species=="Tomato" & metadata$PE==0], yieldSize=2000000)
SlycIHBamFiles <- BamFileList(metadata$Path[metadata$Species=="Tomato" & metadata$PE==1], yieldSize=50000) #lower yield size?
seqinfo(TAIR10BamFiles[1]) #check that it worked

#Import (or create and save) Transcript databases
tryCatch(Slyctxdb <- loadDb("DEGAnalysis/SlycTxDb.sqlite"),error=function(e){
  Slyctxdb <- makeTxDbFromGFF("SlycDNA/ITAG4.0_gene_models.gff", organism="Solanum lycopersicum")
  saveDb(Slyctxdb, "DEGAnalysis/SlycTxDb.sqlite")
})
tryCatch(TAIR10txdb <- loadDb("DEGAnalysis/TAIR10TxDb.sqlite"), error=function(e){
  TAIR10txdb <- makeTxDbFromGFF("ExternalData/TAIR10/TAIR10.gff3", organism="Arabidopsis thaliana")
  saveDb(TAIR10txdb, "DEGAnalysis/TAIR10TxDb.sqlite")
})
tryCatch(Nobttxdb <- loadDb("DEGAnalysis/NobtTxDb.sqlite"),error=function(e){
  Nobttxdb <- makeTxDbFromGFF("NobtDNA/NIOBT_r1.0.update.gff", organism="Nicotiana obtusifolia")
  saveDb(Nobttxdb, "DEGAnalysis/NobtTxDb.sqlite")
})

# Load a gene list by exon for counting (or make and save one)
tryCatch(Slycgenes <- readRDS("DEGAnalysis/Slycgenes.rds"), error=function(e){
  Slycgenes <- exonsBy(Slyctxdb, by="gene")
  saveRDS(Slycgenes, "DEGAnalysis/Slycgenes.rds")
})
tryCatch(TAIR10genes <- readRDS("DEGAnalysis/TAIR10genes.rds"), error=function(e){
  TAIR10genes <- exonsBy(TAIR10txdb, by="gene")
  saveRDS(TAIR10genes, "DEGAnalysis/TAIR10genes.rds")
})
tryCatch(Nobtgenes <- readRDS("DEGAnalysis/Nobtgenes.rds"), error=function(e){
  Nobtgenes <- exonsBy(Nobttxdb, by="gene")
  saveRDS(Nobtgenes, "DEGAnalysis/Nobtgenes.rds")
})


# Count Reads -------------------------------------------------------------
#To my knowledge the SRA experiments were not strand-specific
tryCatch(TAIR10Expt <- readRDS("DEGAnalysis/TAIR10Expt.rds"), error=function(e){
  TAIR10Expt <- summarizeOverlaps(features=TAIR10genes,
                                  reads=TAIR10BamFiles,
                                  mode="Union",
                                  singleEnd=TRUE,
                                  ignore.strand=TRUE)
  colData(TAIR10Expt) <- DataFrame(metadata[metadata$Species=="Arabidopsis",])
  saveRDS(TAIR10Expt, "DEGAnalysis/TAIR10Expt.rds")
})
tryCatch(SlycSRAExpt <- readRDS("DEGAnalysis/SlycSRAExpt.rds"), error=function(e){
  SlycSRAExpt <- summarizeOverlaps(features=Slycgenes,
                                   reads=SlycSRABamFiles,
                                   mode="Union",
                                   singleEnd=TRUE,
                                   ignore.strand=TRUE)
  colData(SlycSRAExpt) <- DataFrame(metadata[metadata$Species=="Tomato" & metadata$PE==0,])
  saveRDS(SlycSRAExpt, "DEGAnalysis/SlycSRAExpt.rds")
})
tryCatch(SlycIHExpt <- readRDS("DEGAnalysis/SlycIHExpt.rds"), error=function(e){
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
  saveRDS(SlycIHExpt, "DEGAnalysis/SlycIHExpt.rds")
})
tryCatch(NobtExpt <- readRDS("DEGAnalysis/NobtExpt.rds"), error=function(e){
  NobtExpt <- summarizeOverlaps(feature=Nobtgenes,
                                reads=NobtBamFiles,
                                mode="Union",
                                singleEnd=FALSE,
                                ignore.strand=FALSE,
                                fragments=TRUE,
                                BPPARAM=SerialParam())
  colData(NobtExpt) <- DataFrame(metadata[metadata$Species=="Tobacco",])
  saveRDS(NobtExpt, "DEGAnalysis/NobtExpt.rds")
})

# Design and DE Testing ----------------------------------------------------
# I wrote a function to do what is currently my default analysis uniformly on every expt
library("splines")
library("DESeq2")
library("ggplot2")

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

# Load (or make and save) DESeq datasets for each experiment
tryCatch(TAIR10dds <- readRDS("DEGAnalysis/TAIR10dds.rds"), error=function(e){
  TAIR10dds <- DESeqSpline(TAIR10Expt)
  saveRDS(TAIR10dds, "DEGAnalysis/TAIR10dds.rds")
})
tryCatch(SlycSRAdds <- readRDS("DEGAnalysis/SlycSRAdds.rds"), error=function(e){
  SlycSRAdds <- DESeqSpline(SlycSRAExpt)
  saveRDS(SlycSRAdds, "DEGAnalysis/SlycSRAdds.rds")
})
tryCatch(SlycIHdds <- readRDS("DEGAnalysis/SlycIHdds.rds"), error=function(e){
  SlycIHdds <- DESeqSpline(SlycIHExpt)
  saveRDS(SlycIHdds, "DEGAnalysis/SlycIHdds.rds")
})
tryCatch(Nobtdds <- readRDS("DEGAnalysis/Nobtdds.rds"), error=function(e){
  Nobtdds <- DESeqSpline(NobtExpt)
  saveRDS(Nobtdds, "DEGAnalysis/Nobtdds.rds")
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
#FUL1 Solyc06g069430.3
#FUL1 Solyc03g114830.3
#MBP10 Solyc02g065730.2
#MBP20 Solyc02g089210.4
plotCounts(Exampledds, gene="Solyc06g069430.3", intgroup="DAP",normalized=T) #FRUITFULL

# Clustering --------------------------------------------------------------
#this section is sorta experimental and is HEAVILY borrowed fromL
#https://hbctraining.github.io/DGE_workshop/lessons/08_DGE_LRT.html
library("magrittr")
library("DEGreport")
library("dplyr")
library("tibble")

DESeqCluster <- function(dds=dds,
                         numGenes=c("1000", "3000", "all"),
                         FDRthreshold=0.01,
                         timeVar="DAP",
                         CaseCtlVar="Genotype"){
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

TAIR10cluster <- DESeqCluster(TAIR10dds, numGenes = "3000")
saveRDS(TAIR10cluster, "DEGAnalysis/TAIR10cluster.rds")

SlycSRAcluster <- DESeqCluster(SlycSRAdds, numGenes = "3000")
saveRDS(SlycSRAcluster, "DEGAnalysis/SlycSRAcluster.rds")


