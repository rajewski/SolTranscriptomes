# Prep Inputs -------------------------------------------------------------
library("GenomicAlignments")
library("Rsamtools")
library("GenomicFeatures")

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
DESeqSpline <- function(se, 
                        dfSpline=3, 
                        timeVar="DAP", 
                        CaseCtl=FALSE, 
                        CaseCtlVar="Genotype",
                        DiagnosticPlots=TRUE) {
  require("splines", "DESeq2", "ggplot2")
  design <- ns(colData(se)[,timeVar], df=dfSpline)
  colnames(design) <- paste0("spline", seq(1:dim(design)[2]))
  colData(se) <- cbind(colData(se), design)
  if (CaseCtl) {
    dds <- DESeqDataSet(se, design = as.formula(paste0("~", CaseCtlVar, "+" ,paste(paste0(CaseCtlVar,":",colnames(design)), collapse = "+"))))
  } else {
    dds <- DESeqDataSet(se, design = as.formula(paste0("~",paste(colnames(design), collapse = "+"))))
  }
  print("Normalizing counts...")
  rld <- rlog(dds)
  print("Done. Now for some diagnostic plots.")
  dds <- estimateSizeFactors(dds)
  plot( assay(rld)[ , 1:2], col=rgb(0,0,0,.2), pch=16, cex=0.3, )
  #This isnt printing because it's a ggplot object, I think
  plotPCA(rld, intgroup = c(if(CaseCtl){CaseCtlVar}, timeVar))
  if (CaseCtl) {
    dds <- DESeq(dds,test="LRT", reduced = as.formula(paste0("~",paste(colnames(design), collapse = "+"))))
  } else {
    dds <- DESeq(dds,test="LRT", reduced = ~ 1)
  }
  print("Getting results...")
  res <- results(dds)
  return(res)
}

#consider analyzing all tomato datasets together
TAIR10res <- DESeqSpline(TAIR10Expt)
SlycSRAres <- DESeqSpline(SlycIHExpt, dfSpline = 2, CaseCtl=FALSE, DiagnosticPlots = TRUE)
SlycSRAdds <- DESeqDataSet(SlycSRAExpt, design = ~ Genotype + Timepoint)
SlycIHdds <- DESeqDataSet(SlycIHExpt, design = ~ Timepoint)
Nobtdds <- DESeqDataSet(NobtExpt, design = ~ Timepoint)


TAIR10res <- results(TAIR10dds)
summary(TAIR10res)

#Apply a FDR cutoff and sort to show the most extreme LFC
TAIR10resSig <- subset(TAIR10res, padj < 0.05)
head(TAIR10resSig[ order( TAIR10resSig$padj ), ])
head(TAIR10resSig[ order( -TAIR10resSig$log2FoldChange ), ])

#Look at other contrasts not reported by default
results(TAIR10dds, contrast=c("DAP", 3, 6))

#Examine an individual Gene
topGene <- rownames(TAIR10res)[which.min(TAIR10res$padj)]
plotCounts(TAIR10dds, gene=topGene, intgroup="DAP", normalized = T)
plotCounts(TAIR10dds, gene="AT5G60910", intgroup=c("DAP"),normalized=T) #FRUITFULL

# Clustering --------------------------------------------------------------
#this section is sorta experimental and is HEAVILY borrowed fromL
#https://hbctraining.github.io/DGE_workshop/lessons/08_DGE_LRT.html
library("magrittr")
BiocManager::install("DEGreport")
library("DEGreport")
library("tibble")
sig_res_TAIR10 <- TAIR10res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < 0.001)
# Subset results for faster cluster finding (for classroom demo purposes) originally at n=1000
clustering_sig_genesTAIR10 <- sig_res_TAIR10 %>%
  arrange(padj) %>%
  head(n=6000)
# Obtain rlog values for those significant genes
cluster_rlog_TAIR10 <- rld[clustering_sig_genesTAIR10$gene, ]
colData(cluster_rlog_TAIR10)$DAP <- as.factor(colData(cluster_rlog_TAIR10)$DAP)
clusters <- degPatterns(assay(cluster_rlog_TAIR10), metadata = colData(cluster_rlog_TAIR10), time = "DAP", col=NULL)
clusters <- degPatterns(assay(cluster_rlog_TAIR10), metadata = colData(cluster_rlog_TAIR10), time = "DAP", col=NULL, reduce = T)
