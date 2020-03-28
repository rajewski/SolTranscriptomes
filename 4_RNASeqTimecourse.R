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
#BiocManager::install("Rsamtools")
library("Rsamtools")
NobtBamFiles <- BamFileList(metadata$Path[metadata$Species=="Tobacco"], yieldSize=2000000)
TAIR10BamFiles <- BamFileList(metadata$Path[metadata$Species=="Arabidopsis"], yieldSize=2000000)
SlycSRABamFiles <- BamFileList(metadata$Path[metadata$Species=="Tomato" & metadata$PE==0], yieldSize=2000000)
SlycIHBamFiles <- BamFileList(metadata$Path[metadata$Species=="Tomato" & metadata$PE==1], yieldSize=50000) #lower yield size?
seqinfo(TAIR10BamFiles[1]) #check that it worked

#BiocManager::install("GenomicFeatures")
library("GenomicFeatures")
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

#Create exon list for each species for counting
(Slycgenes <- exonsBy(Slyctxdb, by="gene"))
(TAIR10genes <- exonsBy(TAIR10txdb, by="gene"))
(Nobtgenes <- exonsBy(Nobttxdb, by="gene"))


# Count Reads -------------------------------------------------------------
library("GenomicAlignments")
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
# I think this is the best way to process the in house (IH) PE strand-specific libraries.
# I could also use preprocess.reads=invertStrand from
# https://support.bioconductor.org/p/65844/ but this seems unusual
# This dataset is too large to read in at once, maybe not with the smaller yield size, but regardless this is the easier way I can find to read it in.
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

# Experimental Design Set up ----------------------------------------------------
#spline regression will be better suited for this
library("splines")
full_design <- model.matrix(formula(~ ns(colData(TAIR10Expt)$DAP, df = 3)))
library("DESeq2")
#consider analyzing all tomato datasets together
TAIR10dds <- DESeqDataSet(TAIR10Expt, design = ~ Timepoint)
TAIR10dds <- estimateSizeFactors(TAIR10dds)
SlycSRAdds <- DESeqDataSet(SlycSRAExpt, design = ~ Genotype + Timepoint)
SlycSRAdds <- estimateSizeFactors(SLYCSRAdds)
SlycIHdds <- DESeqDataSet(SlycIHExpt, design = ~ Timepoint)
SlycIHdds <- estimateSizeFactors(SlycIHdds)
Nobtdds <- DESeqDataSet(NobtExpt, design = ~ Timepoint)
Nobtdds <- estimateSizeFactors(Nobtdds)


# Exploratory Analyses ----------------------------------------------------
#This should be run for each Expt, so I am making the code once with dummy variables
#that can map back to each expt.
dds <- TAIR10dds
#plot a comparison of log2 and regularized log transformations
#the rlog is for exploratory analyses NOT for differential expression
rld <- rlog(dds)
par( mfrow = c( 1, 2 ) )
dds <- estimateSizeFactors(dds)
plot( log2( 1 + counts(dds, normalized=TRUE)[ , 1:2] ),
      col=rgb(0,0,0,.2), pch=16, cex=0.3 )
plot( assay(rld)[ , 1:2],
      col=rgb(0,0,0,.2), pch=16, cex=0.3 )

#Check of the euclidean distances between samples/treatments make sense
sampleDists <- dist( t( assay(rld) ) )
library("gplots")
library("RColorBrewer")
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$Genotype, rld$Timepoint, sep="-" )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
hc <- hclust(sampleDists)
heatmap.2( sampleDistMatrix, Rowv=as.dendrogram(hc),
           symm=TRUE, trace="none", col=colors,
           margins=c(2,10), labCol=FALSE )

#Check again but with Poisson distances
library("PoiClaClu")
poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( dds$Genotype, dds$Timepoint, sep="-" )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
hc <- hclust(poisd$dd)
heatmap.2( samplePoisDistMatrix, Rowv=as.dendrogram(hc),
           symm=TRUE, trace="none", col=colors,
           margins=c(2,10), labCol=FALSE )

#Make a PCA
plotPCA(rld, intgroup = c("Genotype", "Timepoint"))

#Make an MDS Plot with rlog counts
library("ggplot2")
mds <- data.frame(cmdscale(sampleDistMatrix))
mds <- data.frame(cbind(mds, colData(rld)))
qplot(X1,X2,color=Timepoint,shape=Genotype,data=mds)

#Make an MDS plot with poisson counts
mds <- data.frame(cmdscale(samplePoisDistMatrix))
mds <- data.frame(cbind(mds, colData(dds)))
qplot(X1,X2,color=Timepoint,shape=Genotype,data=mds)


# Differential Expression -------------------------------------------------

#do the actual DE testing
TAIR10dds <- DESeq(TAIR10dds)

#Summarize the results
TAIR10res <- results(TAIR10dds)
summary(TAIR10res)

#Apply a FDR cutoff and sort to show the most extreme LFC
TAIR10resSig <- subset(TAIR10res, padj < 0.1)
head(TAIR10resSig[ order( TAIR10resSig$log2FoldChange ), ])
head(TAIR10resSig[ order( -TAIR10resSig$log2FoldChange ), ])

#Look at other contrasts not reported by default
results(TAIR10dds, contrast=c("Timepoint", "DAP3", "DAP6"))

#Examine an individual Gene
topGene <- rownames(TAIR10res)[which.min(TAIR10res$padj)]
plotCounts(dds, gene=topGene, intgroup=c("Timepoint"))
plotCounts(dds, gene="AT5G60910", intgroup=c("Timepoint")) #FRUITFULL

#Look at dispersion. Notsure how useful this is
plotDispEsts(TAIR10dds)

