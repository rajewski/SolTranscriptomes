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
library("tidyr")
source("X_Functions.R")

# Prep Inputs -------------------------------------------------------------
# Borrowed heavily from https://www.bioconductor.org/help/course-materials/2015/LearnBioconductorFeb2015/B02.1.1_RNASeqLab.html#construct

# Load a gene list by exon for counting (or make and save one)
Slycgenes <- tryCatch(readRDS("DEGAnalysis/RNA-seq/Slycgenes.rds"),
                      error=function(e){
                        Slyctxdb <- makeTxDbFromGFF("SlycDNA/ITAG4.0_gene_models.gff",
                                                    organism="Solanum lycopersicum")
                        Slycgenes <- exonsBy(Slyctxdb, by="tx", use.names=TRUE)
                        saveRDS(Slycgenes, "DEGAnalysis/RNA-seq/Slycgenes.rds")
                        return(Slycgenes)})
TAIR10genes <- tryCatch(readRDS("DEGAnalysis/RNA-seq/TAIR10genes.rds"),
                        error=function(e){
                          TAIR10txdb <- makeTxDbFromGFF("ExternalData/TAIR10/TAIR10.gff3",
                                                        organism="Arabidopsis thaliana")
                          TAIR10genes <- exonsBy(TAIR10txdb, by="tx", use.names=TRUE)
                          saveRDS(TAIR10genes, "DEGAnalysis/RNA-seq/TAIR10genes.rds")
                          return(TAIR10genes)})
Nobtgenes <- tryCatch(readRDS("DEGAnalysis/RNA-seq/Nobtgenes.rds"),
                      error=function(e){
                        Nobttxdb <- makeTxDbFromGFF("NobtDNA/NIOBT_r1.0.update.gff",
                                                    organism="Nicotiana obtusifolia")
                        Nobtgenes <- exonsBy(Nobttxdb, by="tx", use.names=TRUE)
                        saveRDS(Nobtgenes, "DEGAnalysis/RNA-seq/Nobtgenes.rds")
                        return(Nobtgenes)})

# Read in the Sample list
metadata <- read.table("DEGAnalysis/RNA-seq/metadata.tsv", header=T, sep="")
metadata$Path <- as.character(metadata$Path)

# Count Reads -------------------------------------------------------------
Expt_TAIR <- tryCatch(readRDS("DEGAnalysis/RNA-seq/Expt_TAIR.rds"),
                      error=function(e){
                        TAIR10BamFiles <- BamFileList(metadata$Path[metadata$Species=="Arabidopsis"],
                                         yieldSize=2000000)
                        Expt_TAIR <- summarizeOverlaps(features=TAIR10genes,
                                                       reads=TAIR10BamFiles,
                                                       mode="Union",
                                                       singleEnd=TRUE,
                                                       ignore.strand=TRUE)
                        colData(Expt_TAIR) <- DataFrame(metadata[metadata$Species=="Arabidopsis",])
                        saveRDS(Expt_TAIR, "DEGAnalysis/RNA-seq/Expt_TAIR.rds")
                        return(Expt_TAIR)})

Expt_Nobt <- tryCatch(readRDS("DEGAnalysis/RNA-seq/Expt_Nobt.rds"),
                      error=function(e){
                        NobtBamFiles <- BamFileList(metadata$Path[metadata$Species=="Tobacco"],
                                                    yieldSize=2000000)
                        Expt_Nobt <- summarizeOverlaps(feature=Nobtgenes,
                                                       reads=NobtBamFiles,
                                                       mode="Union",
                                                       singleEnd=FALSE,
                                                       ignore.strand=FALSE,
                                                       fragments=TRUE,
                                                       BPPARAM=SerialParam())
                        colData(Expt_Nobt) <- DataFrame(metadata[metadata$Species=="Tobacco",])
                        saveRDS(Expt_Nobt, "DEGAnalysis/RNA-seq/Expt_Nobt.rds")
                        return(Expt_Nobt)})

Expt_NobtSE <- tryCatch(readRDS("DEGAnalysis/RNA-seq/Expt_NobtSE.rds"),
                        error=function(e){
                          NobtSEBamFiles <- BamFileList(metadata$Path[metadata$Species=="Tobacco" & metadata$PE==0],
                                                        yieldSize = 2000000)
                          Expt_NobtSE <- summarizeOverlaps(features=Nobtgenes,
                                                           reads=NobtSEBamFiles,
                                                           mode="Union",
                                                           singleEnd=TRUE,
                                                           ignore.strand=FALSE) 
                          colData(Expt_NobtSE) <- DataFrame(metadata[metadata$Species=="Tobacco" & metadata$PE==0,])
                          saveRDS(Expt_NobtSE, "DEGAnalysis/RNA-seq/Expt_NobtSE.rds")
                          return(Expt_NobtSE)})

# Combine the two Nobt datasets
Expt_Nobt_All <- tryCatch(readRDS("DEGAnalysis/RNA-seq/Expt_Nobt_All.rds"),
                          error=function(e){
                            Expt_Nobt <- readRDS("DEGAnalysis/RNA-seq/Expt_Nobt.rds")
                            Expt_NobtSE <- readRDS("DEGAnalysis/RNA-seq/Expt_NobtSE.rds")
                            Expt_Nobt_All <- do.call(cbind, list(Expt_Nobt, Expt_NobtSE))
                            saveRDS(Expt_Nobt_All, "DEGAnalysis/RNA-seq/Expt_Nobt_All.rds")
                            return(Expt_Nobt_All)})

Expt_SlycSRA <- tryCatch(readRDS("DEGAnalysis/RNA-seq/Expt_SlycSRA.rds"),
                         error=function(e){
                           SlycSRABamFiles <- BamFileList(metadata$Path[metadata$Species=="Tomato" & metadata$PE==0],
                                                          yieldSize=2000000)
                           Expt_SlycSRA <- summarizeOverlaps(features=Slycgenes,
                                                             reads=SlycSRABamFiles,
                                                             mode="Union",
                                                             singleEnd=TRUE,
                                                             ignore.strand=TRUE)
                           colData(Expt_SlycSRA) <- DataFrame(metadata[metadata$Species=="Tomato" & metadata$PE==0,])
                           saveRDS(Expt_SlycSRA, "DEGAnalysis/RNA-seq/Expt_SlycSRA.rds")
                           return(Expt_SlycSRA)})

Expt_Slyc <- tryCatch(readRDS("DEGAnalysis/RNA-seq/Expt_Slyc.rds"),
                      error=function(e){
                        SlycIHBamFiles <- BamFileList(metadata$Path[metadata$Species=="Tomato" & metadata$PE==1],
                                         yieldSize=50000)
                        SlycIHPart<-list()
                        for (i in 1:length(metadata$Accession[metadata$Species=="Tomato" & metadata$PE==1])) {
                          SlycIHPart[[i]] <- summarizeOverlaps(feature=Slycgenes,
                                                               reads=SlycIHBamFiles[i],
                                                               mode="Union",
                                                               singleEnd=FALSE,
                                                               ignore.strand=FALSE)}
                        Expt_Slyc <- do.call(cbind,SlycIHPart)
                        colData(Expt_Slyc) <- DataFrame(metadata[metadata$Species=="Tomato" & metadata$PE==1,])
                        saveRDS(Expt_Slyc, "DEGAnalysis/RNA-seq/Expt_Slyc.rds")
                        return(Expt_Slyc)})

Expt_Spimp <- tryCatch(readRDS("DEGAnalysis/RNA-seq/Expt_Spimp.rds"),
         error=function(e){
           SpimpBamFiles <- BamFileList(metadata$Path[metadata$Species=="Pimpinellifolium"],
                                        yieldSize=2000000)
           SpimpPart<-list()
           for (i in 1:length(metadata$Accession[metadata$Species=="Pimpinellifolium"])) {
             SpimpPart[[i]] <- summarizeOverlaps(feature=Slycgenes,
                                                 reads=SpimpBamFiles[i],
                                                 mode="Union",
                                                 singleEnd=FALSE,
                                                 ignore.strand=FALSE)}
           Expt_Spimp <- do.call(cbind,SpimpPart)
           colData(Expt_Spimp) <- DataFrame(metadata[metadata$Species=="Pimpinellifolium",])
           saveRDS(Expt_Spimp, "DEGAnalysis/RNA-seq/Expt_Spimp.rds")
           return(Expt_Spimp)})

# Combine Spimp and Slyc to look for genes that are different between them
Expt_Solanum <- do.call(cbind, list(Expt_Slyc, Expt_Spimp))

# Orthogroups -------------------------------------------------------------
# This section is an attempt to do a cross-species comparison using the orthogroups assigned to the 4 species by Orthofinder.
Expt_All_Ortho <- tryCatch(readRDS("DEGAnalysis/RNA-seq/Expt_All_Ortho.rds"),
                           error=function(e){
                             Expt_Nobt_Ortho <- ConvertGenes2Orthos(OrthogroupMappingFile = "Orthofinder/OrthoFinder/Results_Apr23/Orthogroups/Orthogroups.tsv",
                                                                    GeneWiseExpt = Expt_Nobt,
                                                                    SingleCopyOrthoOnly = TRUE)
                             Expt_Slyc_Ortho <- ConvertGenes2Orthos(OrthogroupMappingFile = "Orthofinder/OrthoFinder/Results_Apr23/Orthogroups/Orthogroups.tsv",
                                                                    GeneWiseExpt = subset(Expt_Slyc, select=DAP<35),
                                                                    SingleCopyOrthoOnly = TRUE)
                             Expt_Spimp_Ortho <- ConvertGenes2Orthos(OrthogroupMappingFile = "Orthofinder/OrthoFinder/Results_Apr23/Orthogroups/Orthogroups.tsv",
                                                                     GeneWiseExpt = subset(Expt_Spimp, select=DAP<35),
                                                                     SingleCopyOrthoOnly = TRUE)
                             Expt_TAIR_Ortho <- ConvertGenes2Orthos(OrthogroupMappingFile = "Orthofinder/OrthoFinder/Results_Apr23/Orthogroups/Orthogroups.tsv",
                                                                    GeneWiseExpt = subset(Expt_TAIR, select=DAP<12),
                                                                    SingleCopyOrthoOnly = TRUE)
                             do.call(cbind, list(Expt_Nobt_Ortho,
                                                 Expt_Slyc_Ortho,
                                                 Expt_Spimp_Ortho,
                                                 Expt_TAIR_Ortho))
                             #Add Stage variable to normalize DAP across species
                             Expt_All_Ortho$Stage <- c(1,1,1,2,2,3,2,3,3,
                                                       1,1,1,2,2,2,3,3,3,
                                                       1,1,1,2,2,2,3,3,3,
                                                       1,1,1,2,2,2,3,3,3)
                             saveRDS(Expt_All_Ortho, file="DEGAnalysis/RNA-seq/Expt_All_Ortho.rds")
                             return(Expt_All_Ortho)})

Expt_NobtRipe_Ortho <- tryCatch(readRDS("DEGAnalysis/RNA-seq/Expt_NobtRipe_Ortho.rds"),
                                error=function(e){
                                  Expt_NobtRipe_Ortho <- ConvertGenes2Orthos(OrthogroupMappingFile = "Orthofinder/OrthoFinder/Results_May17/Orthogroups/Orthogroups.tsv",
                                                                         GeneWiseExpt = subset(Expt_Nobt_All, select=DAP>3),
                                                                         SingleCopyOrthoOnly = TRUE)
                                  Expt_NobtRipe_Ortho$Stage <- c(3,3,3,3.5,3.5,3.5,3.5,3.5,3.5)
                                  saveRDS(Expt_NobtRipe_Ortho, "DEGAnalysis/RNA-seq/Expt_NobtRipe_Ortho.rds")
                                  return(Expt_NobtRipe_Ortho)})
Expt_SlycRipe_Ortho <- tryCatch(readRDS("DEGAnalysis/RNA-seq/Expt_SlycRipe_Ortho.rds"),
                                error=function(e){
                                  Expt_SlycRipe_Ortho <- ConvertGenes2Orthos(OrthogroupMappingFile = "Orthofinder/OrthoFinder/Results_May17/Orthogroups/Orthogroups.tsv",
                                                                             GeneWiseExpt = subset(Expt_Slyc, select=DAP>15),
                                                                             SingleCopyOrthoOnly = TRUE)
                                  Expt_SlycRipe_Ortho$Stage <- c(3,3,3,3.5,3.5,3.5)
                                  saveRDS(Expt_SlycRipe_Ortho, "DEGAnalysis/RNA-seq/Expt_SlycRipe_Ortho.rds")
                                  return(Expt_SlycRipe_Ortho)})
Expt_Ripe_Ortho <- tryCatch(readRDS("DEGAnalysis/RNA-seq/Expt_Ripe_Ortho.rds"),
                            error=function(e){
                              Expt_Nobt_Ortho <- readRDS("DEGAnalysis/RNA-seq/Expt_NobtRipe_Ortho.rds")
                              Expt_Slyc_Ortho <- readRDS("DEGAnalysis/RNA-seq/Expt_SlycRipe_Ortho.rds")
                              Expt_Ripe_Ortho <- do.call(cbind, list(Expt_Nobt_Ortho,
                                                                     Expt_Slyc_Ortho))
                              saveRDS(Expt_Ripe_Ortho, "DEGAnalysis/RNA-seq/Expt_Ripe_Ortho.rds")
                              return(Expt_Ripe_Ortho)})

# Design and DE Testing ----------------------------------------------------
DDS_TAIR <- tryCatch(readRDS("DEGAnalysis/RNA-seq/DDS_TAIR.rds"),
                     error=function(e){
                       DDS_TAIR <- DESeqSpline(Expt_TAIR)
                       saveRDS(DDS_TAIR, "DEGAnalysis/RNA-seq/DDS_TAIR.rds")
                       return(DDS_TAIR)})

DDS_Nobt <- tryCatch(readRDS("DEGAnalysis/RNA-seq/DDS_Nobt.rds"),
                     error=function(e){
                       DDS_Nobt <- DESeqSpline(Expt_Nobt)
                       saveRDS(DDS_Nobt, "DEGAnalysis/RNA-seq/DDS_Nobt.rds")
                       return(DDS_Nobt)})

DDS_Nobt_All <- tryCatch(read.RDS("DEGAnalysis/RNA-seq/DDS_Nobt_All.rds"),
                         error=function(e){
                           DDS_Nobt_All <- DESeqSpline(Expt_Nobt_All)
                           saveRDS(DDS_Nobt_All, "DEGAnalysis/RNA-seq/DDS_Nobt_All.rds")
                           return(DDS_Nobt_All)})

DDS_SlycSRA <- tryCatch(readRDS("DEGAnalysis/RNA-seq/DDS_SlycSRA.rds"),
                        error=function(e){
                          DDS_SlycSRA <- DESeqSpline(Expt_SlycSRA)
                          saveRDS(DDS_SlycSRA, "DEGAnalysis/RNA-seq/DDS_SlycSRA.rds")
                          return(DDS_SlycSRA)})

DDS_Slyc <- tryCatch(readRDS("DEGAnalysis/RNA-seq/DDS_Slyc.rds"),
                     error=function(e){
                       DDS_Slyc <- DESeqSpline(Expt_Slyc)
                       saveRDS(DDS_Slyc, "DEGAnalysis/RNA-seq/DDS_Slyc.rds")
                       return(DDS_Slyc)})

DDS_Spimp <- tryCatch(readRDS("DEGAnalysis/RNA-seq/DDS_Spimp.rds"),
                      error=function(e){
                        DDS_Spimp <- DESeqSpline(Expt_Spimp)
                        saveRDS(DDS_Spimp, "DEGAnalysis/RNA-seq/DDS_Spimp.rds")
                        return(DDS_Spimp)})

DDS_Solanum <- tryCatch(readRDS("DEGAnalysis/RNA-seq/DDS_Solanum.rds"),
                        error=function(e){
                          DDS_Solanum <- DESeqSpline(se=Expt_Solanum,
                                                     CaseCtlVar = "Species")
                          saveRDS(DDS_Solanum, "DEGAnalysis/RNA-seq/DDS_Solanum.rds")
                          return(DDS_Solanum)})

DDS_Solanum_3DF <- tryCatch(readRDS("DEGAnalysis/RNA-seq/DDS_Solanum_3DF.rds"),
                            error=function(e){
                              DDS_Solanum <- DESeqSpline(se=Expt_Solanum,
                                                         CaseCtlVar = "Species",
                                                         SetDF = 3)
                              saveRDS(DDS_Solanum_3DF, "DEGAnalysis/RNA-seq/DDS_Solanum_3DF.rds")
                              return(DDS_Solanum_3DF)})

DDS_Solanum_3DF_Noise <- tryCatch(readRDS("DEGAnalysis/RNA-seq/DDS_Solanum_3DF_Noise.rds"),
                                  error=function(e){
                                    DESeqSpline(Expt_Solanum,
                                                vsNoise = TRUE,
                                                SetDF = 3)
                                    saveRDS(DDS_Solanum_3DF_Noise, "DEGAnalysis/RNA-seq/DDS_Solanum_3DF_Noise.rds")
                                    return(DDS_Solanum_3DF_Noise)})

# Test for DEGs with diff patterns by species
DDS_AllOrtho_DEGBySpecies <- tryCatch(readRDS("DEGAnalysis/RNA-seq/DDS_AllOrtho_DEGBySpecies.rds"),
                                      error=function(e){
                                        DDS_AllOrtho_DEGBySpecies <- DESeqSpline(Expt_All_Ortho,
                                                                                 CaseCtlVar = "Species",
                                                                                 timeVar = "Stage")
                                        # columns have duplicate names, which would be changed and mess up metadata mapping
                                        colnames(DDS_AllOrtho_DEGBySpecies) <- NULL 
                                        saveRDS(DDS_AllOrtho_DEGBySpecies, "DEGAnalysis/RNA-seq/DDS_AllOrtho_DEGBySpecies.rds")
                                        return(DDS_AllOrtho_DEGBySpecies)})

# Test for DEGs with diff patterns by fruit type
DDS_AllOrtho_DEGByFruit <- tryCatch(readRDS("DEGAnalysis/RNA-seq/DDS_AllOrtho_DEGByFruit.rds"),
                                    error=function(e){
                                      DDS_AllOrtho_DEGByFruit <- DESeqSpline(Expt_All_Ortho,
                                                                             CaseCtlVar = "Fruit",
                                                                             timeVar = "Stage")
                                      # columns have duplicate names, which would be changed and mess up metadata mapping
                                      colnames(DDS_AllOrtho_DEGByFruit) <- NULL 
                                      saveRDS(DDS_AllOrtho_DEGByFruit, "DEGAnalysis/RNA-seq/DDS_AllOrtho_DEGByFruit.rds")
                                      return(DDS_AllOrtho_DEGByFruit)})

# Find the genes with a common pattern across species
DDS_AllOrtho_Noise <- tryCatch(readRDS("DEGAnalysis/RNA-seq/DDS_AllOrtho_Noise.rds"),
                               error=function(e){
                                 DESeqSpline(Expt_All_Ortho,
                                             vsNoise = TRUE,
                                             timeVar="Stage")
                                 colnames(DDS_AllOrtho_Noise) <- NULL 
                                 saveRDS(DDS_AllOrtho_Noise, "DEGAnalysis/RNA-seq/DDS_AllOrtho_Noise.rds")
                                 return(DDS_AllOrtho_Noise)})

# Test for DEGs in Dry Fruited only
DDS_DryOrtho <- tryCatch(readRDS("DEGAnalysis/RNA-seq/DDS_DryOrtho.rds"),
                         error=function(e){
                           DDS_DryOrtho <- DESeqSpline(Expt_Dry_Ortho,
                                                       CaseCtlVar = "Species",
                                                       timeVar = "Stage")
                           # columns have duplicate names, which would be changed and mess up metadata mapping
                           colnames(DDS_DryOrtho) <- NULL 
                           saveRDS(DDS_DryOrtho, "DEGAnalysis/RNA-seq/DDS_DryOrtho.rds")
                           return(DDS_DryOrtho)})

# Test for DE orthogenes at Ripening for Nobt only
DDS_NobtRipe_Ortho <- tryCatch(readRDS("DEGAnalysis/RNA-seq/DDS_NobtRipe_Ortho.rds"),
                               error=function(e){
                                 DDS_NobtRipe_Ortho <- DESeqDataSetFromMatrix(countData = assays(Expt_NobtRipe_Ortho)$counts,
                                                        colData = colData(Expt_NobtRipe_Ortho),
                                                        design = ~ Stage)
                                 DDS_NobtRipe_Ortho <- estimateSizeFactors(DDS_NobtRipe_Ortho)
                                 DDS_NobtRipe_Ortho <- DESeq(DDS_NobtRipe_Ortho,
                                                             test="LRT",
                                                             reduced = ~1)
                                 saveRDS(DDS_NobtRipe_Ortho, "DEGAnalysis/RNA-seq/DDS_NobtRipe_Ortho.rds")
                                 return(DDS_NobtRipe_Ortho)})
# Test for DE orthogenes at Ripening for Slyc Only
DDS_SlycRipe_Ortho <- tryCatch(readRDS("DEGAnalysis/RNA-seq/DDS_SlycRipe_Ortho.rds"),
                               error=function(e){
                                 DDS_SlycRipe_Ortho <- DESeqDataSetFromMatrix(countData = assays(Expt_SlycRipe_Ortho)$counts,
                                                        colData = colData(Expt_SlycRipe_Ortho),
                                                        design = ~ Stage)
                                 DDS_SlycRipe_Ortho <- estimateSizeFactors(DDS_SlycRipe_Ortho)
                                 DDS_SlycRipe_Ortho <- DESeq(DDS_SlycRipe_Ortho,
                                                             test="LRT",
                                                             reduced = ~1)
                                 saveRDS(DDS_SlycRipe_Ortho, "DEGAnalysis/RNA-seq/DDS_SlycRipe_Ortho.rds")
                                 return(DDS_SlycRipe_Ortho)})

# Play around with individual genes ---------------------------------------
Exampledds <- DDS_AllOrtho_DEGByFruit #assign one dds as the example to streamline code
ExampleRes <- results(Exampledds) #get results
ExampleResSig <- subset(ExampleRes, padj < 0.05) #subset by FDR
head(ExampleResSig[order(ExampleResSig$padj ), ]) #see best fitting genes for spline model
#Examine an individual Gene
topGene <- rownames(ExampleRes)[which.min(ExampleRes$padj)]
colData(Exampledds)$DAP <- as.factor(colData(Exampledds)$DAP)
#colData(Exampledds)$Stage <- as.factor(colData(Exampledds)$Stage)
plotCounts(Exampledds, gene=topGene, intgroup=c("Species", "DAP"), normalized = T) #plot best fitting gene

# Get a set of FUL genes for each species. Only use one of these
FULgenes<-c(FUL.1="AT5G60910.1",
            FUL.2="AT5G60910.2",
            AGL79="AT3G30260.1")
FULgenes<-c(SlFUL1="Solyc06g069430.3.1",
            SlFUL2="Solyc03g114830.3.1",
            SlMBP10="Solyc02g065730.2.1",
            SlMBP20="Solyc02g089210.4.1" )
FULgenes<-c(SpFUL1="Solyc06g069430.3.1",
            SpFUL2="Solyc03g114830.3.1",
            SpMBP10="Solyc02g065730.2.1",
            SpMBP20="Solyc02g089210.4.1" )
FULgenes<-c(NoFUL1="NIOBTv3_g28929-D2.t1",
            NoFUL2="NIOBTv3_g39464.t1",
            NoMBP10="NIOBTv3_g07845.t1",
            NoMBP20="NIOBT_gMBP20.t1" )
FULgenes<-c(euFULI="OG0003276",
            euFULII="OG0008754")
# Plot FUL Genes
for (i in 1:length(FULgenes)) {
  pdf(file=paste0("DEGAnalysis/RNA-seq/Plots/Slyc_", names(FULgenes[i]), ".pdf"), # _SRA v _IH on Slyc
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
Cluster_TAIR <- tryCatch(readRDS("DEGAnalysis/RNA-seq/Cluster_TAIR.rds"),
                         error=function(e){
                           Cluster_TAIR <- DESeqCluster(DDS_TAIR, numGenes = "all")
                           saveRDS(Cluster_TAIR, "DEGAnalysis/RNA-seq/Cluster_TAIR.rds")
                           return(Cluster_TAIR)})

Cluster_Nobt <- tryCatch(readRDS("DEGAnalysis/RNA-seq/Cluster_Nobt.rds"),
                         error=function(e){
                           Cluster_Nobt <- DESeqCluster(DDS_Nobt, numGenes = "all")
                           saveRDS(Cluster_Nobt, "DEGAnalysis/RNA-seq/Cluster_Nobt.rds")
                           return(Cluster_Nobt)})

Cluster_SlycSRA <- tryCatch(readRDS("DEGAnalysis/RNA-seq/Cluster_SlycSRA.rds"),
                            error=function(e){
                              Cluster_SlycSRA <- DESeqCluster(DDS_SlycSRA, numGenes = "all")
                              saveRDS(Cluster_SlycSRA, "DEGAnalysis/RNA-seq/Cluster_SlycSRA.rds")
                              return(Cluster_SlycSRA)})

Cluster_Slyc <- tryCatch(readRDS("DEGAnalysis/RNA-seq/Cluster_Slyc.rds"),
                         error=function(e){
                           Cluster_Slyc <- DESeqCluster(DDS_Slyc, numGenes = "all")
                           saveRDS(Cluster_Slyc, "DEGAnalysis/RNA-seq/Cluster_Slyc.rds")
                           return(Cluster_Slyc)})

Cluster_Spimp <- tryCatch(readRDS("DEGAnalysis/RNA-seq/Cluster_Spimp.rds"),
                          error=function(e){
                            Cluster_Spimp <- DESeqCluster(DDS_Spimp, numGenes = "all")
                            saveRDS(Cluster_Spimp, "DEGAnalysis/RNA-seq/Cluster_Spimp.rds")
                            return(Cluster_Spimp)})

Cluster_Solanum <- tryCatch(readRDS("DEGAnalysis/RNA-seq/Cluster_Solanum.rds"),
                            error=function(e){
                              Cluster_Solanum <- DESeqCluster(DDS_Solanum,  
                                                              numGenes = "all",
                                                              CaseCtlVar = "Species")
                              saveRDS(Cluster_Solanum, "DEGAnalysis/RNA-seq/Cluster_Solanum.rds")
                              return(Cluster_Solanum)})

Cluster_Solanum_3DF <- tryCatch(readRDS("DEGAnalysis/RNA-seq/Cluster_Solanum_3DF.rds"),
                                error=function(e){
                                  Cluster_Solanum_3DF <- DESeqCluster(DDS_Solanum_3DF,
                                                                      numGenes = "all",
                                                                      CaseCtlVar = "Species")
                                  saveRDS(Cluster_Solanum_3DF, "DEGAnalysis/RNA-seq/Cluster_Solanum_3DF.rds")
                                  return(Cluster_Solanum_3DF)})

Cluster_Solanum_3DF_Noise <- tryCatch(readRDS("DEGAnalysis/RNA-seq/Cluster_Solanum_3DF_Noise.rds"),
                                      error=function(e){
                                        Cluster_Solanum_3DF_Noise <- DESeqCluster(DDS_Solanum_3DF_Noise,
                                                                                  numGenes = "all")
                                        saveRDS(Cluster_Solanum_3DF_Noise, "DEGAnalysis/RNA-seq/Cluster_Solanum_3DF_Noise.rds")
                                        return(Cluster_Solanum_3DF_Noise)})

Cluster_AllOrtho_DEGBySpecies <- tryCatch(readRDS("DEGAnalysis/RNA-seq/Cluster_AllOrtho_DEGBySpecies.rds"),
                                          error=function(e){
                                            Cluster_AllOrtho_DEGBySpecies <- DESeqCluster(DDS_AllOrtho_DEGBySpecies,
                                                                                          numGenes = "all",
                                                                                          CaseCtlVar = "Species",
                                                                                          timeVar = "Stage")
                                            saveRDS(Cluster_AllOrtho_DEGBySpecies, "DEGAnalysis/RNA-seq/Cluster_AllOrtho_DEGBySpecies.rds")
                                            return(Cluster_AllOrtho_DEGBySpecies)})

Cluster_AllOrtho_DEGByFruit <- tryCatch(readRDS("DEGAnalysis/RNA-seq/Cluster_AllOrtho_DEGByFruit.rds"),
                                        error=function(e){
                                          Cluster_AllOrtho_DEGByFruit <- DESeqCluster(DDS_AllOrtho_DEGByFruit,
                                                                                      numGenes = "all",
                                                                                      CaseCtlVar = "Fruit",
                                                                                      timeVar = "Stage")
                                          saveRDS(Cluster_AllOrtho_DEGByFruit, "DEGAnalysis/RNA-seq/Cluster_AllOrtho_DEGByFruit.rds")
                                          return(Cluster_AllOrtho_DEGByFruit)})

Cluster_AllOrtho_Noise <- tryCatch(readRDS("DEGAnalysis/RNA-seq/Cluster_AllOrtho_Noise.rds"), 
                                   error=function(e){
                                     DESeqCluster(DDS_AllOrtho_Noise,
                                                  numGenes = "all",
                                                  timeVar = "Stage")
                                     saveRDS(Cluster_AllOrtho_Noise, "DEGAnalysis/RNA-seq/Cluster_AllOrtho_Noise.rds")
                                     return(Cluster_AllOrtho_Noise)})

Cluster_DryOrtho <- tryCatch(readRDS("DEGAnalysis/RNA-seq/Cluster_DryOrtho.rds"),
                             error=function(e){
                               Cluster_DryOrtho <- DESeqCluster(DDS_DryOrtho,
                                                                numGenes = "all",
                                                                CaseCtlVar = "Species",
                                                                timeVar = "Stage")
                               saveRDS(Cluster_DryOrtho, "DEGAnalysis/RNA-seq/Cluster_DryOrtho.rds")
                               return(Cluster_DryOrtho)})


# Plot Cluster Profiles ---------------------------------------------------
ClusterforPlotting <- Cluster_AllOrtho_Noise
PlotCluster <-degPlotCluster(ClusterforPlotting$normalized[ClusterforPlotting$normalized$cluster==2,],
                             time="Stage",
                             boxes=T,
                             points=F,
                             #color="Species",
                             lines=F
                             )
PlotCluster + theme_minimal() +
  theme(#legend.position = c(1, 0.0),
        #legend.justification = c(1, 0),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave(filename = "DEGAnalysis/RNA-seq/Plots/ClusterProfiles_AllOrtho_Noise.pdf",
       width=11,
       height=7)

# NoFUL2 and NoMBP10 are in cluster 2
# In SlycIH FUL1 is cluster 3
# in SlycSRA FUL1 and FUL2 are cluster 1
# FUL and AGL79 are not DE

# Save Cluster genes to a file --------------------------------------------
X <- split(Cluster_AllOrtho_Noise$df, Cluster_AllOrtho_Noise$df$cluster)
for (i in 1:length(X)) {
  write.table(row.names(X[[i]]), 
              file=paste0("DEGAnalysis/RNA-seq/Lists/AllOrtho_Noise_Cluster_", max(X[[i]]$cluster), ".txt"),
              row.names = FALSE,
              quote = FALSE,
              col.names = FALSE)
}

