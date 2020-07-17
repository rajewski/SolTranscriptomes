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
metadata <- read.table("DEGAnalysis/RNA-seq/metadata.tsv", header=T, sep="")
metadata$Path <- as.character(metadata$Path)
# Get list of BAMs for each expt
NobtBamFiles <- BamFileList(metadata$Path[metadata$Species=="Tobacco"], yieldSize=2000000)
NobtSEBamFiles <- BamFileList(metadata$Path[metadata$Species=="Tobacco" & metadata$PE==0], yieldSize = 2000000)
TAIR10BamFiles <- BamFileList(metadata$Path[metadata$Species=="Arabidopsis"], yieldSize=2000000)
SlycSRABamFiles <- BamFileList(metadata$Path[metadata$Species=="Tomato" & metadata$PE==0], yieldSize=2000000)
SlycIHBamFiles <- BamFileList(metadata$Path[metadata$Species=="Tomato" & metadata$PE==1], yieldSize=50000)
SpimpBamFiles <- BamFileList(metadata$Path[metadata$Species=="Pimpinellifolium"], yieldSize=2000000)
seqinfo(TAIR10BamFiles[1]) #check that it worked

# Count Reads -------------------------------------------------------------
#To my knowledge the SRA experiments were not strand-specific
tryCatch(Expt_TAIR <- readRDS("DEGAnalysis/RNA-seq/Expt_TAIR.rds"),
         error=function(e){
           Expt_TAIR <- summarizeOverlaps(features=TAIR10genes,
                                  reads=TAIR10BamFiles,
                                  mode="Union",
                                  singleEnd=TRUE,
                                  ignore.strand=TRUE)
  colData(Expt_TAIR) <- DataFrame(metadata[metadata$Species=="Arabidopsis",])
  saveRDS(Expt_TAIR, "DEGAnalysis/RNA-seq/Expt_TAIR.rds")
})
#Option to subset to stages 1-3
Expt_TAIR_3Stage <- subset(Expt_TAIR, select=DAP<12)

tryCatch(Expt_Nobt <- readRDS("DEGAnalysis/RNA-seq/Expt_Nobt.rds"),
         error=function(e){
           Expt_Nobt <- summarizeOverlaps(feature=Nobtgenes,
                                reads=NobtBamFiles,
                                mode="Union",
                                singleEnd=FALSE,
                                ignore.strand=FALSE,
                                fragments=TRUE,
                                BPPARAM=SerialParam())
  colData(Expt_Nobt) <- DataFrame(metadata[metadata$Species=="Tobacco",])
  saveRDS(Expt_Nobt, "DEGAnalysis/RNA-seq/Expt_Nobt.rds")
})

tryCatch(Expt_NobtSE <- readRDS("DEGAnalysis/RNA-seq/Expt_NobtSE.rds"),
         error=function(e){
           Expt_NobtSE <- summarizeOverlaps(features=Nobtgenes,
                                            reads=NobtSEBamFiles,
                                            mode="Union",
                                            singleEnd=TRUE,
                                            ignore.strand=FALSE) 
           colData(Expt_NobtSE) <- DataFrame(metadata[metadata$Species=="Tobacco" & metadata$PE==0,])
           saveRDS(Expt_NobtSE, "DEGAnalysis/RNA-seq/Expt_NobtSE.rds")
         })
#Combine the two Nobt datasets
Expt_Nobt_All <- do.call(cbind, list(Expt_Nobt, Expt_NobtSE))
saveRDS(Expt_Nobt_All, "DEGAnalysis/RNA-seq/Expt_Nobt_All.rds")
tryCatch(Expt_SlycSRA <- readRDS("DEGAnalysis/RNA-seq/Expt_SlycSRA.rds"),
         error=function(e){
           Expt_SlycSRA <- summarizeOverlaps(features=Slycgenes,
                                   reads=SlycSRABamFiles,
                                   mode="Union",
                                   singleEnd=TRUE,
                                   ignore.strand=TRUE)
  colData(Expt_SlycSRA) <- DataFrame(metadata[metadata$Species=="Tomato" & metadata$PE==0,])
  saveRDS(Expt_SlycSRA, "DEGAnalysis/RNA-seq/Expt_SlycSRA.rds")
})
tryCatch(Expt_Slyc <- readRDS("DEGAnalysis/RNA-seq/Expt_Slyc.rds"),
         error=function(e){
  SlycIHPart<-list()
  for (i in 1:length(metadata$Accession[metadata$Species=="Tomato" & metadata$PE==1])) {
    SlycIHPart[[i]] <- summarizeOverlaps(feature=Slycgenes,
                                         reads=SlycIHBamFiles[i],
                                         mode="Union",
                                         singleEnd=FALSE,
                                         ignore.strand=FALSE)
  }
  Expt_Slyc <- do.call(cbind,SlycIHPart)
  colData(Expt_Slyc) <- DataFrame(metadata[metadata$Species=="Tomato" & metadata$PE==1,])
  saveRDS(Expt_Slyc, "DEGAnalysis/RNA-seq/Expt_Slyc.rds")
})
#Option to subset to stages 1-3
Expt_Slyc_3Stage <- subset(Expt_Slyc, select=DAP<35)
tryCatch(Expt_Spimp <- readRDS("DEGAnalysis/RNA-seq/Expt_Spimp.rds"),
         error=function(e){
  SpimpPart<-list()
  for (i in 1:length(metadata$Accession[metadata$Species=="Pimpinellifolium"])) {
    SpimpPart[[i]] <- summarizeOverlaps(feature=Slycgenes,
                                         reads=SpimpBamFiles[i],
                                         mode="Union",
                                         singleEnd=FALSE,
                                         ignore.strand=FALSE)
  }
  Expt_Spimp <- do.call(cbind,SpimpPart)
  colData(Expt_Spimp) <- DataFrame(metadata[metadata$Species=="Pimpinellifolium",])
  saveRDS(Expt_Spimp, "DEGAnalysis/RNA-seq/Expt_Spimp.rds")
})
#Option to subset to stages 1-3
Expt_Spimp_3Stage <- subset(Expt_Spimp, select=DAP<35)
#Combine Spimp and Slyc to look for genes that are different between them
Expt_Solanum <- do.call(cbind, list(Expt_Slyc, Expt_Spimp))
#Expt_Solanum_3Stage <- do.call(cbind, list(Expt_Slyc_3Stage, Expt_Spimp_3Stage))

# Orthogroups -------------------------------------------------------------
# This section is an attempt to do a cross-species comparison using the orthogroups assigned to the 4 species by Orthofinder.
Expt_Nobt_Ortho <- ConvertGenes2Orthos(OrthogroupMappingFile = "Orthofinder/OrthoFinder/Results_Apr23/Orthogroups/Orthogroups.tsv",
                                     GeneWiseExpt = Expt_Nobt,
                                     SingleCopyOrthoOnly = TRUE)
Expt_Slyc_Ortho <- ConvertGenes2Orthos(OrthogroupMappingFile = "Orthofinder/OrthoFinder/Results_Apr23/Orthogroups/Orthogroups.tsv",
                                     GeneWiseExpt = Expt_Slyc_3Stage,
                                     SingleCopyOrthoOnly = TRUE)
Expt_Spimp_Ortho <- ConvertGenes2Orthos(OrthogroupMappingFile = "Orthofinder/OrthoFinder/Results_Apr23/Orthogroups/Orthogroups.tsv",
                                      GeneWiseExpt = Expt_Spimp_3Stage,
                                      SingleCopyOrthoOnly = TRUE)
Expt_TAIR_Ortho <- ConvertGenes2Orthos(OrthogroupMappingFile = "Orthofinder/OrthoFinder/Results_Apr23/Orthogroups/Orthogroups.tsv",
                                     GeneWiseExpt = Expt_TAIR_3Stage,
                                     SingleCopyOrthoOnly = TRUE)
Expt_All_Ortho <- do.call(cbind, list(Expt_Nobt_Ortho,
                                      Expt_Slyc_Ortho,
                                      Expt_Spimp_Ortho,
                                      Expt_TAIR_Ortho))
Expt_Dry_Ortho <- do.call(cbind, list(Expt_Nobt_Ortho,
                                      Expt_TAIR_Ortho))
#Add Stage variable to normalize DAP across species
# I could do this in the metadata, but then I would have to recount everything and that would take forever
Expt_All_Ortho$Stage <- c(1,1,1,2,2,3,2,3,3,
                          1,1,1,2,2,2,3,3,3,
                          1,1,1,2,2,2,3,3,3,
                          1,1,1,2,2,2,3,3,3)
saveRDS(Expt_All_Ortho, file="DEGAnalysis/RNA-seq/Expt_All_Ortho.rds")
Expt_Dry_Ortho$Stage <- c(1,1,1,2,2,3,2,3,3,
                          1,1,1,2,2,2,3,3,3)

# Design and DE Testing ----------------------------------------------------
tryCatch(DDS_TAIR <- readRDS("DEGAnalysis/RNA-seq/DDS_TAIR.rds"),
         error=function(e){
           DDS_TAIR <- DESeqSpline(Expt_TAIR)
  saveRDS(DDS_TAIR, "DEGAnalysis/RNA-seq/DDS_TAIR.rds")
})
tryCatch(DDS_TAIR_3Stage <- readRDS("DEGAnalysis/RNA-seq/DDS_TAIR_3Stage.rds"),
         error=function(e){
           DDS_TAIR_3Stage <- DESeqSpline(Expt_TAIR_3Stage)
  saveRDS(DDS_TAIR_3Stage, "DEGAnalysis/RNA-seq/DDS_TAIR_3Stage.rds")
})
tryCatch(DDS_Nobt <- readRDS("DEGAnalysis/RNA-seq/DDS_Nobt.rds"),
         error=function(e){
           DDS_Nobt <- DESeqSpline(Expt_Nobt)
  saveRDS(DDS_Nobt, "DEGAnalysis/RNA-seq/DDS_Nobt.rds")
})

tryCatch(DDS_Nobt_All <- read.RDS("DEGAnalysis/RNA-seq/DDS_Nobt_All.rds"),
         error=function(e){
           DDS_Nobt_All <- DESeqSpline(Expt_Nobt_All)
           saveRDS(DDS_Nobt_All, "DEGAnalysis/RNA-seq/DDS_Nobt_All.rds")
         })
tryCatch(DDS_SlycSRA <- readRDS("DEGAnalysis/RNA-seq/DDS_SlycSRA.rds"),
         error=function(e){
           DDS_SlycSRA <- DESeqSpline(Expt_SlycSRA)
           saveRDS(DDS_SlycSRA, "DEGAnalysis/RNA-seq/DDS_SlycSRA.rds")
})
tryCatch(DDS_Slyc <- readRDS("DEGAnalysis/RNA-seq/DDS_Slyc.rds"),
         error=function(e){
           DDS_Slyc <- DESeqSpline(Expt_Slyc)
  saveRDS(DDS_Slyc, "DEGAnalysis/RNA-seq/DDS_Slyc.rds")
})
tryCatch(DDS_Slyc_3Stage <- readRDS("DEGAnalysis/RNA-seq/DDS_Slyc_3Stage.rds"),
         error=function(e){
           DDS_Slyc_3Stage <- DESeqSpline(Expt_Slyc_3Stage)
  saveRDS(DDS_Slyc_3Stage, "DEGAnalysis/RNA-seq/DDS_Slyc_3Stage.rds")
})
tryCatch(DDS_Spimp <- readRDS("DEGAnalysis/RNA-seq/DDS_Spimp.rds"),
         error=function(e){
           DDS_Spimp <- DESeqSpline(Expt_Spimp)
  saveRDS(DDS_Spimp, "DEGAnalysis/RNA-seq/DDS_Spimp.rds")
})
tryCatch(DDS_Spimp_3Stage <- readRDS("DEGAnalysis/RNA-seq/DDS_Spimp_3Stage.rds"),
         error=function(e){
           DDS_Spimp_3Stage <- DESeqSpline(Expt_Spimp_3Stage)
           saveRDS(DDS_Spimp_3Stage, "DEGAnalysis/RNA-seq/DDS_Spimp_3Stage.rds")
})
tryCatch(DDS_Solanum <- readRDS("DEGAnalysis/RNA-seq/DDS_Solanum.rds"),
         error=function(e){
           DDS_Solanum <- DESeqSpline(se=Expt_Solanum,
                                     CaseCtlVar = "Species")
           saveRDS(DDS_Solanum, "DEGAnalysis/RNA-seq/DDS_Solanum.rds")
         })
DDS_Solanum_3DF <- DESeqSpline(se=Expt_Solanum,
                              CaseCtlVar = "Species",
                              SetDF = 3)
saveRDS(DDS_Solanum_3DF, "DEGAnalysis/RNA-seq/DDS_Solanum_3DF.rds")
tryCatch(DDS_Solanum_3Stage <- readRDS("DEGAnalysis/RNA-seq/DDS_Solanum_3Stage.rds"),
         error=function(e){
           DDS_Solanum_3Stage <- DESeqSpline(se=Expt_Solanum_3Stage,
                                           CaseCtlVar = "Species")
           saveRDS(DDS_Solanum_3Stage, "DEGAnalysis/RNA-seq/DDS_Solanum_3Stage.rds")
         })
# Test for DEGs with diff patterns by species
tryCatch(DDS_AllOrtho_DEGBySpecies <- readRDS("DEGAnalysis/RNA-seq/DDS_AllOrtho_DEGBySpecies.rds"),
         error=function(e){
           DDS_AllOrtho_DEGBySpecies <- DESeqSpline(Expt_All_Ortho,
                                             CaseCtlVar = "Species",
                                             timeVar = "Stage")
           # columns have duplicate names, which would be changed and mess up metadata mapping
           colnames(DDS_AllOrtho_DEGBySpecies) <- NULL 
           saveRDS(DDS_AllOrtho_DEGBySpecies, "DEGAnalysis/RNA-seq/DDS_AllOrtho_DEGBySpecies.rds")
         })
# Test for DEGs with diff patterns by fruit type
tryCatch(DDS_AllOrtho_DEGByFruit <- readRDS("DEGAnalysis/RNA-seq/DDS_AllOrtho_DEGByFruit.rds"),
         error=function(e){
           DDS_AllOrtho_DEGByFruit <- DESeqSpline(AllSpeciesOrthoExpt,
                                                  CaseCtlVar = "Fruit",
                                                  timeVar = "Stage")
           # columns have duplicate names, which would be changed and mess up metadata mapping
           colnames(DDS_AllOrtho_DEGByFruit) <- NULL 
           saveRDS(DDS_AllOrtho_DEGByFruit, "DEGAnalysis/RNA-seq/DDS_AllOrtho_DEGByFruit.rds")
         })
# Test for DEGs in Dry Fruited only
tryCatch(DDS_DryOrtho <- readRDS("DEGAnalysis/RNA-seq/DDS_DryOrtho.rds"),
         error=function(e){
           DDS_DryOrtho <- DESeqSpline(Expt_Dry_Ortho,
                                       CaseCtlVar = "Species",
                                       timeVar = "Stage")
           # columns have duplicate names, which would be changed and mess up metadata mapping
           colnames(DDS_DryOrtho) <- NULL 
           saveRDS(DDS_DryOrtho, "DEGAnalysis/RNA-seq/DDS_DryOrtho.rds")
         })

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
tryCatch(Cluster_TAIR <- readRDS("DEGAnalysis/RNA-seq/Cluster_TAIR.rds"),
         error=function(e){
           Cluster_TAIR <- DESeqCluster(DDS_TAIR, numGenes = "all")
  saveRDS(Cluster_TAIR, "DEGAnalysis/RNA-seq/Cluster_TAIR.rds")
})
tryCatch(Cluster_TAIR_3Stage <- readRDS("DEGAnalysis/RNA-seq/Cluster_TAIR_3Stage.rds"),
         error=function(e){
           Cluster_TAIR_3Stage <- DESeqCluster(DDS_TAIR_3Stage, numGenes = "all")
  saveRDS(Cluster_TAIR_3Stage, "DEGAnalysis/RNA-seq/Cluster_TAIR_3Stage.rds")
})
tryCatch(Cluster_Nobt <- readRDS("DEGAnalysis/RNA-seq/Cluster_Nobt.rds"),
         error=function(e){
           Cluster_Nobt <- DESeqCluster(DDS_Nobt, numGenes = "all")
           saveRDS(Cluster_Nobt, "DEGAnalysis/RNA-seq/Cluster_Nobt.rds")
         })
tryCatch(Cluster_SlycSRA <- readRDS("DEGAnalysis/RNA-seq/Cluster_SlycSRA.rds"),
         error=function(e){
           Cluster_SlycSRA <- DESeqCluster(DDS_SlycSRA, numGenes = "all")
  saveRDS(Cluster_SlycSRA, "DEGAnalysis/RNA-seq/Cluster_SlycSRA.rds")
})
tryCatch(Cluster_Slyc <- readRDS("DEGAnalysis/RNA-seq/Cluster_Slyc.rds"),
         error=function(e){
           Cluster_Slyc <- DESeqCluster(DDS_Slyc, numGenes = "all")
  saveRDS(Cluster_Slyc, "DEGAnalysis/RNA-seq/Cluster_Slyc.rds")
})
tryCatch(Cluster_Slyc_3Stage <- readRDS("DEGAnalysis/RNA-seq/Cluster_Slyc_3Stage.rds"),
         error=function(e){
           Cluster_Slyc_3Stage <- DESeqCluster(DDS_Slyc_3Stage, numGenes = "all")
  saveRDS(Cluster_Slyc_3Stage, "DEGAnalysis/RNA-seq/Cluster_Slyc_3Stage.rds")
})
tryCatch(Cluster_Spimp <- readRDS("DEGAnalysis/RNA-seq/Cluster_Spimp.rds"),
         error=function(e){
           Cluster_Spimp <- DESeqCluster(DDS_Spimp, numGenes = "all")
  saveRDS(Cluster_Spimp, "DEGAnalysis/RNA-seq/Cluster_Spimp.rds")
})
tryCatch(Cluster_Spimp_3Stage <- readRDS("DEGAnalysis/RNA-seq/Cluster_Spimp_3Stage.rds"),
         error=function(e){
           Cluster_Spimp_3Stage <- DESeqCluster(DDS_Spimp_3Stage, numGenes = "all")
  saveRDS(Cluster_Spimp_3Stage, "DEGAnalysis/RNA-seq/Cluster_Spimp_3Stage.rds")
})
tryCatch(Cluster_Solanum <- readRDS("DEGAnalysis/RNA-seq/Cluster_Solanum.rds"),
         error=function(e){
           Cluster_Solanum <- DESeqCluster(DDS_Solanum,
                                          numGenes = "all",
                                          CaseCtlVar = "Species")
           saveRDS(Cluster_Solanum, "DEGAnalysis/RNA-seq/Cluster_Solanum.rds")
         })
tryCatch(Cluster_Solanum_3Stage <- readRDS("DEGAnalysis/RNA-seq/Cluster_Solanum_3Stage.rds"),
         error=function(e){
           Cluster_Solanum_3Stage <- DESeqCluster(DDS_Solanum_3Stage,
                                          numGenes = "all",
                                          CaseCtlVar = "Species")
           saveRDS(Cluster_Solanum_3Stage, "DEGAnalysis/RNA-seq/Cluster_Solanum_3Stage.rds")
         })
tryCatch(Cluster_Solanum_3DF <- readRDS("DEGAnalysis/RNA-seq/Cluster_Solanum_3DF.rds"),
         error=function(e){
           Cluster_Solanum_3DF <- DESeqCluster(DDS_Solanum_3DF,
                                                  numGenes = "all",
                                                  CaseCtlVar = "Species")
           saveRDS(Cluster_Solanum_3DF, "DEGAnalysis/RNA-seq/Cluster_Solanum_3DF.rds")
         })
tryCatch(Cluster_AllOrtho_DEGBySpecies <- readRDS("DEGAnalysis/RNA-seq/Cluster_AllOrtho_DEGBySpecies.rds"),
         error=function(e){
           Cluster_AllOrtho_DEGBySpecies <- DESeqCluster(DDS_AllOrtho_DEGBySpecies,
                                                  numGenes = "all",
                                                  CaseCtlVar = "Species",
                                                  timeVar = "Stage")
           saveRDS(Cluster_AllOrtho_DEGBySpecies, "DEGAnalysis/RNA-seq/Cluster_AllOrtho_DEGBySpecies.rds")
         })
tryCatch(Cluster_AllOrtho_DEGByFruit <- readRDS("DEGAnalysis/RNA-seq/Cluster_AllOrtho_DEGByFruit.rds"),
         error=function(e){
           Cluster_AllOrtho_DEGByFruit <- DESeqCluster(DDS_AllOrtho_DEGByFruit,
                                                       numGenes = "all",
                                                       CaseCtlVar = "Fruit",
                                                       timeVar = "Stage")
           saveRDS(Cluster_AllOrtho_DEGByFruit, "DEGAnalysis/RNA-seq/Cluster_AllOrtho_DEGByFruit.rds")
         })

tryCatch(Cluster_DryOrtho <- readRDS("DEGAnalysis/RNA-seq/Cluster_DryOrtho.rds"),
         error=function(e){
           Cluster_DryOrtho <- DESeqCluster(DDS_DryOrtho,
                                            numGenes = "all",
                                            CaseCtlVar = "Species",
                                            timeVar = "Stage")
           saveRDS(Cluster_DryOrtho, "DEGAnalysis/RNA-seq/Cluster_DryOrtho.rds")
         })


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


# NonDE Genes -------------------------------------------------------------

# I want to look at the genes that show the same expression pattern in some interspecific comparison. I'll start with SP and AC together using the vsNoise setting
DDS_Solanum_3DF_Noise <- DESeqSpline(Expt_Solanum,
                                     vsNoise = TRUE,
                                     SetDF = 3)
tryCatch(Cluster_Solanum_3DF_Noise <- readRDS("DEGAnalysis/RNA-seq/Cluster_Solanum_3DF_Noise.rds"),
         error=function(e){
           Cluster_Solanum_3DF_Noise <- DESeqCluster(DDS_Solanum_3DF_Noise,
                                                     numGenes = "all")
           saveRDS(Cluster_Solanum_3DF_Noise, "DEGAnalysis/RNA-seq/Cluster_Solanum_3DF_Noise.rds")
         })


# Do it with the pre-ripening data set
DDS_AllOrtho_Noise <- DESeqSpline(Expt_All_Ortho,
                                  vsNoise = TRUE)
colnames(DDS_AllOrtho_Noise) <- NULL 
saveRDS(DDS_AllOrtho_Noise, "DEGAnalysis/RNA-seq/DDS_AllOrtho_Noise.rds")
Cluster_AllOrtho_Noise <- DESeqCluster(DDS_AllOrtho_Noise,
                                       numGenes = "all",
                                       timeVar = "Stage")
saveRDS(Cluster_AllOrtho_Noise, "DEGAnalysis/RNA-seq/Cluster_AllOrtho_Noise.rds")
