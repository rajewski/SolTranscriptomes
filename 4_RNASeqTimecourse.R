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
library("BiocParallel")
source("X_Functions.R")

# Prep Inputs -------------------------------------------------------------
# Borrowed heavily from https://www.bioconductor.org/help/course-materials/2015/LearnBioconductorFeb2015/B02.1.1_RNASeqLab.html#construct

# Load a gene list by exon for counting (or make and save one)
Melongenes <- tryCatch(readRDS("DEGAnalysis/RNA-seq/Melongenes.rds"),
                       error=function(e){
                         Melontxdb <- makeTxDbFromGFF("ExternalData/C_melo/CM3.5.1_gene.gff",
                                                      organism="Cucumis melo")
                         Melongenes <- exonsBy(Melontxdb, by="tx", use.names=TRUE)
                         saveRDS(Melongenes, "DEGAnalysis/RNA-seq/Melongenes.rds")
                         return(Melongenes)
                       })
Nobtgenes <- tryCatch(readRDS("DEGAnalysis/RNA-seq/Nobtgenes.rds"),
                      error=function(e){
                        Nobttxdb <- makeTxDbFromGFF("NobtDNA/NIOBT_r1.0.update.gff",
                                                    organism="Nicotiana obtusifolia")
                        Nobtgenes <- exonsBy(Nobttxdb, by="tx", use.names=TRUE)
                        saveRDS(Nobtgenes, "DEGAnalysis/RNA-seq/Nobtgenes.rds")
                        return(Nobtgenes)})

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
# Read in the Sample list
metadata <- read.table("DEGAnalysis/RNA-seq/metadata.tsv", header=T, sep="")
metadata$Path <- as.character(metadata$Path)

# Count Reads -------------------------------------------------------------
Expt_Nobt <- tryCatch(readRDS("DEGAnalysis/RNA-seq/Expt_Nobt.rds"), #PE stage 1-3 data
                      error=function(e){
                        NobtBamFiles <- BamFileList(metadata$Path[metadata$Species=="Tobacco" & metadata$PE==1], yieldSize=2000000)
                        Expt_Nobt <- summarizeOverlaps(feature=Nobtgenes,
                                                       reads=NobtBamFiles,
                                                       mode="Union",
                                                       singleEnd=FALSE,
                                                       ignore.strand=FALSE,
                                                       BPPARAM=SerialParam()) 
                        colData(Expt_Nobt) <- DataFrame(metadata[metadata$Species=="Tobacco" & metadata$PE==1,])
                        saveRDS(Expt_Nobt, "DEGAnalysis/RNA-seq/Expt_Nobt.rds")
                        return(Expt_Nobt)})

Expt_NobtSE <- tryCatch(readRDS("DEGAnalysis/RNA-seq/Expt_NobtSE.rds"), #All data as SE
                        error=function(e){
                          NobtSEBamFiles <- BamFileList(metadata$Path[metadata$Species=="Tobacco" & metadata$PE==0], yieldSize = 2000000)
                          Expt_NobtSE <- summarizeOverlaps(features=Nobtgenes,
                                                           reads=NobtSEBamFiles,
                                                           mode="Union",
                                                           singleEnd=TRUE,
                                                           ignore.strand=FALSE) 
                          colData(Expt_NobtSE) <- DataFrame(metadata[metadata$Species=="Tobacco" & metadata$PE==0,])
                          saveRDS(Expt_NobtSE, "DEGAnalysis/RNA-seq/Expt_NobtSE.rds")
                          return(Expt_NobtSE)})

Expt_Melon <- tryCatch(readRDS("DEGAnalysis/RNA-seq/Expt_Melon.rds"),
                      error=function(e){
                        MelonBamFiles <- BamFileList(metadata$Path[metadata$Species=="Melon"],
                                                      yieldSize=2000000)
                        Expt_Melon <- summarizeOverlaps(features=Melongenes,
                                                       reads=MelonBamFiles,
                                                       mode="Union",
                                                       singleEnd=TRUE,
                                                       ignore.strand=TRUE,
                                                       BPPARAM=SerialParam())
                        colData(Expt_Melon) <- DataFrame(metadata[metadata$Species=="Melon",])
                        saveRDS(Expt_Melon, "DEGAnalysis/RNA-seq/Expt_Melon.rds")
                        return(Expt_Melon)})

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

Expt_SlycSE <- tryCatch(readRDS("DEGAnalysis/RNA-seq/Expt_SlycSE.rds"),
                        error=function(e){
                          SlycSEBamFiles <- BamFileList(metadata$Path[grep(x = metadata$Path,pattern="*Slyc_SE*")],
                                                        yieldSize=50000)
                          Expt_SlycSE <- summarizeOverlaps(feature=Slycgenes,
                                                           reads=SlycSEBamFiles,
                                                           mode="Union",
                                                           singleEnd=TRUE,
                                                           ignore.strand=FALSE)
                          colData(Expt_SlycSE) <- DataFrame(metadata[grep(x = metadata$Path,pattern="*Slyc_SE*"),])
                          saveRDS(Expt_SlycSE, "DEGAnalysis/RNA-seq/Expt_SlycSE.rds")
                          return(Expt_SlycSE)})

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

Expt_SpimpSE <- tryCatch(readRDS("DEGAnalysis/RNA-seq/Expt_SpimpSE.rds"),
                         error=function(e){
                           SpimpSEBamFiles <- BamFileList(metadata$Path[grep(x = metadata$Path,pattern="*Spimp_SE*")],
                                                          yieldSize=50000)
                           Expt_SpimpSE <- summarizeOverlaps(feature=Slycgenes,
                                                             reads=SpimpSEBamFiles,
                                                             mode="Union",
                                                             singleEnd=TRUE,
                                                             ignore.strand=FALSE,
                                                             BPPARAM=SerialParam())
                           colData(Expt_SpimpSE) <- DataFrame(metadata[grep(x = metadata$Path,pattern="*Spimp_SE*"),])
                           saveRDS(Expt_SpimpSE, "DEGAnalysis/RNA-seq/Expt_SpimpSE.rds")
                           return(Expt_SpimpSE)}) #Alex Redo these 8/5/20

Expt_TAIR <- tryCatch(readRDS("DEGAnalysis/RNA-seq/Expt_TAIR.rds"),
                      error=function(e){
                        TAIR10BamFiles <- BamFileList(metadata$Path[metadata$Species=="Arabidopsis"],
                                         yieldSize=2000000)
                        Expt_TAIR <- summarizeOverlaps(features=TAIR10genes,
                                                       reads=TAIR10BamFiles,
                                                       mode="Union",
                                                       singleEnd=TRUE,
                                                       ignore.strand=TRUE,
                                                       BPPARAM=SerialParam())
                        colData(Expt_TAIR) <- DataFrame(metadata[metadata$Species=="Arabidopsis",])
                        saveRDS(Expt_TAIR, "DEGAnalysis/RNA-seq/Expt_TAIR.rds")
                        return(Expt_TAIR)})

# Combine Spimp and Slyc to look for genes that are different between them
Expt_Solanum <- tryCatch(readRDS("DEGAnalysis/RNA-seq/Expt_Solanum.rds"),
                        error=function(e){
                        do.call(cbind, list(Expt_Slyc, Expt_Spimp))})

# Orthogroups -------------------------------------------------------------
# This section does a cross-species comparison using the orthogroups assigned to the 4 species by Orthofinder.
Expt_All_Ortho <- tryCatch(readRDS("DEGAnalysis/RNA-seq/Expt_All_Ortho.rds"),
                  error=function(e){
                   mapping <- "Orthofinder/OrthoFinder/Results_May17/Orthogroups/Orthogroups.tsv"
                   Expt_Nobt_Ortho <- ConvertGenes2Orthos(OrthogroupMappingFile = mapping,
                                         GeneWiseExpt = Expt_NobtSE,
                                         SingleCopyOrthoOnly = TRUE)
                   Expt_Slyc_Ortho <- ConvertGenes2Orthos(OrthogroupMappingFile = mapping,
                                         GeneWiseExpt = subset(Expt_SlycSE, select=DAP<45),
                                         SingleCopyOrthoOnly = TRUE)
                   Expt_Spimp_Ortho <- ConvertGenes2Orthos(OrthogroupMappingFile = mapping,
                                         GeneWiseExpt = subset(Expt_SpimpSE, select=DAP<45),
                                         SingleCopyOrthoOnly = TRUE)
                   Expt_TAIR_Ortho <- ConvertGenes2Orthos(OrthogroupMappingFile = mapping,
                                         GeneWiseExpt = Expt_TAIR,
                                         SingleCopyOrthoOnly = TRUE)
                   Expt_All_Ortho <- do.call(cbind, list(Expt_Nobt_Ortho,
                                       Expt_Slyc_Ortho,
                                       Expt_Spimp_Ortho,
                                       Expt_TAIR_Ortho))
                   #Add Stage variable to normalize DAP across species
                   Expt_All_Ortho$Stage <- c(3.5,3.5,3.5,3.5,3.5,3.5,1,1,1,2,2,2,3,3,3,
                                             1,1,1,2,2,2,3,3,3,3.5,3.5,3.5,
                                             1,1,1,2,2,2,3,3,3,3.5,3.5,3.5,
                                             1,1,1,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3.5,3.5,3.5)
                   saveRDS(Expt_All_Ortho, file="DEGAnalysis/RNA-seq/Expt_All_Ortho.rds")
                   return(Expt_All_Ortho)})

Expt_All_Ortho_Unripe <- tryCatch(readRDS("DEGAnalysis/RNA-seq/Expt_All_Ortho_Unripe.rds"),
                           error=function(e){
                             mapping <- "Orthofinder/OrthoFinder/Results_May17/Orthogroups/Orthogroups.tsv"
                             Expt_Nobt_Ortho <- ConvertGenes2Orthos(OrthogroupMappingFile = mapping,
                                                  GeneWiseExpt = subset(Expt_NobtSE, select=DAP<11),
                                                  SingleCopyOrthoOnly = TRUE)
                             Expt_Slyc_Ortho <- ConvertGenes2Orthos(OrthogroupMappingFile = mapping,
                                                  GeneWiseExpt = subset(Expt_SlycSE, select=DAP<35),
                                                  SingleCopyOrthoOnly = TRUE)
                             Expt_Spimp_Ortho <- ConvertGenes2Orthos(OrthogroupMappingFile = mapping,
                                                  GeneWiseExpt = subset(Expt_SpimpSE, select=DAP<35),
                                                  SingleCopyOrthoOnly = TRUE)
                             Expt_TAIR_Ortho <- ConvertGenes2Orthos(OrthogroupMappingFile = mapping,
                                                  GeneWiseExpt = subset(Expt_TAIR, select=DAP<12),
                                                  SingleCopyOrthoOnly = TRUE)
                             Expt_All_Ortho_Unripe <- do.call(cbind, list(Expt_Nobt_Ortho,
                                                  Expt_Slyc_Ortho,
                                                  Expt_Spimp_Ortho,
                                                  Expt_TAIR_Ortho))
                             #Add Stage variable to normalize DAP across species
                             Expt_All_Ortho_Unripe$Stage <- c(1,1,1,2,2,2,3,3,3,
                                                              1,1,1,2,2,2,3,3,3,
                                                              1,1,1,2,2,2,3,3,3,
                                                              1,1,1,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3)
                             saveRDS(Expt_All_Ortho_Unripe, file="DEGAnalysis/RNA-seq/Expt_All_Ortho_Unripe.rds")
                             return(Expt_All_Ortho_Unripe)})

Expt_Five_Ortho_Unripe <- tryCatch(readRDS("DEGAnalysis/RNA-seq/Expt_Five_Ortho_Unripe.rds"),
                                  error=function(e){
                                    mapping <- "Orthofinder/OrthoFinder/Results_Aug31/Orthogroups/Orthogroups.tsv"
                    Expt_Nobt_Ortho <- ConvertGenes2Orthos(OrthogroupMappingFile = mapping,
                                                           GeneWiseExpt = subset(Expt_NobtSE, select=DAP<11),
                                                           SingleCopyOrthoOnly = TRUE)
                    Expt_Melon_Ortho <- ConvertGenes2Orthos(OrthogroupMappingFile = mapping,
                                                            GeneWiseExpt = subset(Expt_Melon, select=DAP<40),
                                                            SingleCopyOrthoOnly = TRUE)
                    Expt_Slyc_Ortho <- ConvertGenes2Orthos(OrthogroupMappingFile = mapping,
                                                           GeneWiseExpt = subset(Expt_SlycSE, select=DAP<35),
                                                           SingleCopyOrthoOnly = TRUE)
                    Expt_Spimp_Ortho <- ConvertGenes2Orthos(OrthogroupMappingFile = mapping,
                                                            GeneWiseExpt = subset(Expt_SpimpSE, select=DAP<35),
                                                            SingleCopyOrthoOnly = TRUE)
                    Expt_TAIR_Ortho <- ConvertGenes2Orthos(OrthogroupMappingFile = mapping,
                                                           GeneWiseExpt = subset(Expt_TAIR, select=DAP<12),
                                                           SingleCopyOrthoOnly = TRUE)
                    Expt_Five_Ortho_Unripe <- do.call(cbind, list(Expt_Nobt_Ortho,
                                                                 Expt_Melon_Ortho,
                                                                 Expt_Slyc_Ortho,
                                                                 Expt_Spimp_Ortho,
                                                                 Expt_TAIR_Ortho))
                    #Add Stage variable to normalize DAP across species
                    Expt_Five_Ortho_Unripe$Stage <- c(1,1,1,2,2,3,3,3,2,
                                                     1,1,1,2,2,2,3,3,3,
                                                     1,1,1,2,2,2,3,3,3,
                                                     1,1,1,2,2,2,3,3,3,
                                                     1,1,1,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3)
                    saveRDS(Expt_Five_Ortho_Unripe, file="DEGAnalysis/RNA-seq/Expt_Five_Ortho_Unripe.rds")
                    return(Expt_Five_Ortho_Unripe)})

Expt_SolOrtho <- tryCatch(readRDS("DEGAnalysis/RNA-seq/Expt_SolOrtho.rds"),
                             error=function(e){
                               mapping <- "Orthofinder/OrthoFinder/Results_May17/Orthogroups/Orthogroups.tsv"
                               Expt_Nobt_Ortho <- ConvertGenes2Orthos(OrthogroupMappingFile = mapping,
                                                      GeneWiseExpt = Expt_NobtSE,
                                                      SingleCopyOrthoOnly = TRUE,
                                                      Arabidopsis = FALSE)
                               Expt_Slyc_Ortho <- ConvertGenes2Orthos(OrthogroupMappingFile = mapping,
                                                      GeneWiseExpt = subset(Expt_SlycSE, select=DAP<45),
                                                      SingleCopyOrthoOnly = TRUE,
                                                      Arabidopsis = FALSE)
                               Expt_Spimp_Ortho <- ConvertGenes2Orthos(OrthogroupMappingFile = mapping,
                                                      GeneWiseExpt = subset(Expt_SpimpSE, select=DAP<45),
                                                      SingleCopyOrthoOnly = TRUE,
                                                      Arabidopsis = FALSE)
                               Expt_SolOrtho <- do.call(cbind, list(Expt_Nobt_Ortho,
                                                                     Expt_Slyc_Ortho,
                                                                     Expt_Spimp_Ortho))
                               #Add Stage variable to normalize DAP across species
                               Expt_SolOrtho$Stage <- c(3.5,3.5,3.5,3.5,3.5,3.5,1,1,1,2,2,2,3,3,3,
                                                         1,1,1,2,2,2,3,3,3,3.5,3.5,3.5,
                                                         1,1,1,2,2,2,3,3,3,3.5,3.5,3.5)
                               saveRDS(Expt_SolOrtho, file="DEGAnalysis/RNA-seq/Expt_SolOrtho.rds")
                               return(Expt_SolOrtho)})

Expt_SolOrtho_Unripe_Unripe <- tryCatch(readRDS("DEGAnalysis/RNA-seq/Expt_SolOrtho_Unripe.rds"),
                          error=function(e){
                            mapping <- "Orthofinder/OrthoFinder/Results_May17/Orthogroups/Orthogroups.tsv"
                            Expt_Nobt_Ortho <- ConvertGenes2Orthos(OrthogroupMappingFile = mapping,
                                                   GeneWiseExpt = subset(Expt_NobtSE, select=DAP<11),
                                                   SingleCopyOrthoOnly = TRUE,
                                                   Arabidopsis = FALSE)
                            Expt_Slyc_Ortho <- ConvertGenes2Orthos(OrthogroupMappingFile = mapping,
                                                   GeneWiseExpt = subset(Expt_SlycSE, select=DAP<35),
                                                   SingleCopyOrthoOnly = TRUE,
                                                   Arabidopsis = FALSE)
                            Expt_Spimp_Ortho <- ConvertGenes2Orthos(OrthogroupMappingFile = mapping,
                                                   GeneWiseExpt = subset(Expt_SpimpSE, select=DAP<35),
                                                   SingleCopyOrthoOnly = TRUE,
                                                   Arabidopsis = FALSE)
                            Expt_SolOrtho_Unripe <- do.call(cbind, list(Expt_Nobt_Ortho,
                                                                 Expt_Slyc_Ortho,
                                                                 Expt_Spimp_Ortho))
                            #Add Stage variable to normalize DAP across species
                            Expt_SolOrtho_Unripe$Stage <- c(1,1,1,2,2,2,3,3,3,
                                                     1,1,1,2,2,2,3,3,3,
                                                     1,1,1,2,2,2,3,3,3)
                            saveRDS(Expt_SolOrtho_Unripe, file="DEGAnalysis/RNA-seq/Expt_SolOrtho_Unripe.rds")
                            return(Expt_SolOrtho_Unripe)})

Expt_NobtRipe_Ortho <- tryCatch(readRDS("DEGAnalysis/RNA-seq/Expt_NobtRipe_Ortho.rds"),
                                error=function(e){
                        Expt_NobtRipe_Ortho <- ConvertGenes2Orthos(OrthogroupMappingFile = "Orthofinder/OrthoFinder/Results_May17/Orthogroups/Orthogroups_NoArabidopsis.tsv",
                                                GeneWiseExpt = subset(Expt_NobtSE, select=DAP>3),
                                                SingleCopyOrthoOnly = TRUE,
                                                Arabidopsis = FALSE)
                        Expt_NobtRipe_Ortho$Stage <- c(3,3,3,3.5,3.5,3.5,3.5,3.5,3.5)
                        saveRDS(Expt_NobtRipe_Ortho, "DEGAnalysis/RNA-seq/Expt_NobtRipe_Ortho.rds")
                        return(Expt_NobtRipe_Ortho)})

Expt_Slyc1545_Ortho <- tryCatch(readRDS("DEGAnalysis/RNA-seq/Expt_Slyc1545_Ortho.rds"),
                                error=function(e){
                         tmp <- subset(Expt_SlycSE, select=DAP>11)
                         tmp <- subset(tmp, select=DAP!=35)
                         Expt_Slyc1545_Ortho <- ConvertGenes2Orthos(OrthogroupMappingFile = "Orthofinder/OrthoFinder/Results_May17/Orthogroups/Orthogroups.tsv",
                                                GeneWiseExpt = tmp,
                                                SingleCopyOrthoOnly = TRUE,
                                                Arabidopsis = FALSE)
                         Expt_Slyc1545_Ortho$Stage <- c(3,3,3,4,4,4)
                         saveRDS(Expt_Slyc1545_Ortho, "DEGAnalysis/RNA-seq/Expt_Slyc1545_Ortho.rds")
                         rm(tmp)
                         return(Expt_Slyc1545_Ortho)})
Expt_Slyc1535_Ortho <- tryCatch(readRDS("DEGAnalysis/RNA-seq/Expt_Slyc1535_Ortho.rds"),
                                error=function(e){
                          tmp <- subset(Expt_SlycSE, select=DAP>11)
                          tmp <- subset(tmp, select=DAP<45)
                          Expt_Slyc1535_Ortho <- ConvertGenes2Orthos(OrthogroupMappingFile = "Orthofinder/OrthoFinder/Results_May17/Orthogroups/Orthogroups.tsv",
                                                GeneWiseExpt = tmp,
                                                SingleCopyOrthoOnly = TRUE,
                                                Arabidopsis = FALSE)
                          Expt_Slyc1535_Ortho$Stage <- c(3,3,3,3.5,3.5,3.5)
                          saveRDS(Expt_Slyc1535_Ortho, "DEGAnalysis/RNA-seq/Expt_Slyc1535_Ortho.rds")
                          rm(tmp)
                          return(Expt_Slyc1535_Ortho)})
Expt_Slyc3545_Ortho <- tryCatch(readRDS("DEGAnalysis/RNA-seq/Expt_Slyc3545_Ortho.rds"),
                                error=function(e){
                                  tmp <- subset(Expt_SlycSE, select=DAP>15)
                                  Expt_Slyc3545_Ortho <- ConvertGenes2Orthos(OrthogroupMappingFile = "Orthofinder/OrthoFinder/Results_May17/Orthogroups/Orthogroups.tsv",
                                                GeneWiseExpt = tmp,
                                                SingleCopyOrthoOnly = TRUE,
                                                Arabidopsis = FALSE)
                                  Expt_Slyc3545_Ortho$Stage <- c(3,3,3,3.5,3.5,3.5)
                                  saveRDS(Expt_Slyc3545_Ortho, "DEGAnalysis/RNA-seq/Expt_Slyc3545_Ortho.rds")
                                  rm(tmp)
                                  return(Expt_Slyc3545_Ortho)})

# Design and DE Testing ----------------------------------------------------
DDS_Melon <- tryCatch(readRDS("DEGAnalysis/RNA-seq/DDS_Melon.rds"),
                      error=function(e){
                        DDS_Melon <- DESeqSpline(Expt_Melon)
                        saveRDS(DDS_Melon, "DEGAnalysis/RNA-seq/DDS_Melon.rds")
                        return(DDS_Melon)})

DDS_Nobt <- tryCatch(readRDS("DEGAnalysis/RNA-seq/DDS_Nobt.rds"),
                     error=function(e){
                       DDS_Nobt <- DESeqSpline(Expt_Nobt)
                       saveRDS(DDS_Nobt, "DEGAnalysis/RNA-seq/DDS_Nobt.rds")
                       return(DDS_Nobt)})

DDS_NobtSE <- tryCatch(readRDS("DEGAnalysis/RNA-seq/DDS_NobtSE.rds"),
                       error=function(e){
                         DDS_NobtSE <- DESeqSpline(Expt_NobtSE,
                                                   CollapseTechRep = TRUE)
                         saveRDS(DDS_NobtSE, "DEGAnalysis/RNA-seq/DDS_NobtSE.rds")
                         return(DDS_NobtSE)})

DDS_Slyc <- tryCatch(readRDS("DEGAnalysis/RNA-seq/DDS_Slyc.rds"),
                     error=function(e){
                       DDS_Slyc <- DESeqSpline(Expt_Slyc)
                       saveRDS(DDS_Slyc, "DEGAnalysis/RNA-seq/DDS_Slyc.rds")
                       return(DDS_Slyc)})

DDS_SlycSRA <- tryCatch(readRDS("DEGAnalysis/RNA-seq/DDS_SlycSRA.rds"),
                        error=function(e){
                          DDS_SlycSRA <- DESeqSpline(Expt_SlycSRA)
                          saveRDS(DDS_SlycSRA, "DEGAnalysis/RNA-seq/DDS_SlycSRA.rds")
                          return(DDS_SlycSRA)})

DDS_Spimp <- tryCatch(readRDS("DEGAnalysis/RNA-seq/DDS_Spimp.rds"),
                      error=function(e){
                        DDS_Spimp <- DESeqSpline(Expt_Spimp)
                        saveRDS(DDS_Spimp, "DEGAnalysis/RNA-seq/DDS_Spimp.rds")
                        return(DDS_Spimp)})

DDS_TAIR <- tryCatch(readRDS("DEGAnalysis/RNA-seq/DDS_TAIR.rds"),
                     error=function(e){
                       DDS_TAIR <- DESeqSpline(Expt_TAIR,
                                               CollapseTechRep = TRUE)
                       saveRDS(DDS_TAIR, "DEGAnalysis/RNA-seq/DDS_TAIR.rds")
                       return(DDS_TAIR)}) 

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
                                                timeVar = "Stage",
                                                CollapseTechRep = TRUE)
                # columns have duplicate names, which would be changed and mess up metadata mapping
                colnames(DDS_AllOrtho_DEGBySpecies) <- NULL 
                saveRDS(DDS_AllOrtho_DEGBySpecies, "DEGAnalysis/RNA-seq/DDS_AllOrtho_DEGBySpecies.rds")
                return(DDS_AllOrtho_DEGBySpecies)})

# Test for DEGs with diff patterns by fruit type
DDS_AllOrtho_DEGByFruit <- tryCatch(readRDS("DEGAnalysis/RNA-seq/DDS_AllOrtho_DEGByFruit.rds"),
                error=function(e){
                  DDS_AllOrtho_DEGByFruit <- DESeqSpline(Expt_All_Ortho,
                                                         CaseCtlVar = "Fruit",
                                                         timeVar = "Stage",
                                                         CollapseTechRep = TRUE)
                  # columns have duplicate names, which would be changed and mess up metadata mapping
                  colnames(DDS_AllOrtho_DEGByFruit) <- NULL 
                  saveRDS(DDS_AllOrtho_DEGByFruit, "DEGAnalysis/RNA-seq/DDS_AllOrtho_DEGByFruit.rds")
                  return(DDS_AllOrtho_DEGByFruit)})

# Find the genes with a common pattern across species
DDS_AllOrtho_Noise <- tryCatch(readRDS("DEGAnalysis/RNA-seq/DDS_AllOrtho_Noise.rds"),
                error=function(e){
                  DDS_AllOrtho_Noise <- DESeqSpline(Expt_All_Ortho,
                                                    vsNoise = TRUE,
                                                    timeVar="Stage",
                                                    CollapseTechRep = TRUE)
                  colnames(DDS_AllOrtho_Noise) <- NULL 
                  saveRDS(DDS_AllOrtho_Noise, "DEGAnalysis/RNA-seq/DDS_AllOrtho_Noise.rds")
                  return(DDS_AllOrtho_Noise)})

#Just the solanaceous species
DDS_SolOrtho_DEGBySpecies <- tryCatch(readRDS("DEGAnalysis/RNA-seq/DDS_SolOrtho_DEGBySpecies.rds"),
                error=function(e){
                 DDS_SolOrtho_DEGBySpecies <- DESeqSpline(Expt_SolOrtho,
                                                          CaseCtlVar = "Species",
                                                          timeVar = "Stage",
                                                          CollapseTechRep = TRUE)
                 # columns have duplicate names, which would be changed and mess up metadata mapping
                 colnames(DDS_SolOrtho_DEGBySpecies) <- NULL 
                 saveRDS(DDS_SolOrtho_DEGBySpecies, "DEGAnalysis/RNA-seq/DDS_SolOrtho_DEGBySpecies.rds")
                 return(DDS_SolOrtho_DEGBySpecies)})
DDS_SolOrtho_DEGByFruit <- tryCatch(readRDS("DEGAnalysis/RNA-seq/DDS_SolOrtho_DEGByFruit.rds"),
               error=function(e){
                 DDS_SolOrtho_DEGByFruit <- DESeqSpline(Expt_SolOrtho,
                                                        CaseCtlVar = "Fruit",
                                                        timeVar = "Stage",
                                                        CollapseTechRep = TRUE)
                 # columns have duplicate names, which would be changed and mess up metadata mapping
                 colnames(DDS_SolOrtho_DEGByFruit) <- NULL 
                 saveRDS(DDS_SolOrtho_DEGByFruit, "DEGAnalysis/RNA-seq/DDS_SolOrtho_DEGByFruit.rds")
                 return(DDS_SolOrtho_DEGByFruit)})
DDS_SolOrtho_Noise <- tryCatch(readRDS("DEGAnalysis/RNA-seq/DDS_SolOrtho_Noise.rds"),
              error=function(e){
                DDS_SolOrtho_Noise <- DESeqSpline(Expt_SolOrtho,
                                                  vsNoise = TRUE,
                                                  timeVar="Stage",
                                                  CollapseTechRep = TRUE)
                colnames(DDS_SolOrtho_Noise) <- NULL 
                saveRDS(DDS_SolOrtho_Noise, "DEGAnalysis/RNA-seq/DDS_SolOrtho_Noise.rds")
                return(DDS_SolOrtho_Noise)})

# Unripe Data -------
# Test for DEGs with diff patterns by species
DDS_AllOrtho_Unripe_DEGBySpecies <- tryCatch(readRDS("DEGAnalysis/RNA-seq/DDS_AllOrtho_Unripe_DEGBySpecies.rds"),
              error=function(e){
              DDS_AllOrtho_Unripe_DEGBySpecies <- DESeqSpline(Expt_All_Ortho_Unripe,
                                                   CaseCtlVar = "Species",
                                                   timeVar = "Stage",
                                                   CollapseTechRep = TRUE)
              # columns have duplicate names, which would be changed and mess up metadata mapping
              colnames(DDS_AllOrtho_Unripe_DEGBySpecies) <- NULL 
              saveRDS(DDS_AllOrtho_Unripe_DEGBySpecies, "DEGAnalysis/RNA-seq/DDS_AllOrtho_Unripe_DEGBySpecies.rds")
              return(DDS_AllOrtho_Unripe_DEGBySpecies)})

# Test for DEGs with diff patterns by fruit type
DDS_AllOrtho_Unripe_DEGByFruit <- tryCatch(readRDS("DEGAnalysis/RNA-seq/DDS_AllOrtho_Unripe_DEGByFruit.rds"),
              error=function(e){
              DDS_AllOrtho_Unripe_DEGByFruit <- DESeqSpline(Expt_All_Ortho_Unripe,
                                    CaseCtlVar = "Fruit",
                                    timeVar = "Stage",
                                    CollapseTechRep = TRUE)
              # columns have duplicate names, which would be changed and mess up metadata mapping
              colnames(DDS_AllOrtho_Unripe_DEGByFruit) <- NULL 
              saveRDS(DDS_AllOrtho_Unripe_DEGByFruit, "DEGAnalysis/RNA-seq/DDS_AllOrtho_Unripe_DEGByFruit.rds")
              return(DDS_AllOrtho_Unripe_DEGByFruit)})

# Find the genes with a common pattern across species
DDS_AllOrtho_Unripe_Noise <- tryCatch(readRDS("DEGAnalysis/RNA-seq/DDS_AllOrtho_Unripe_Noise.rds"),
                               error=function(e){
                                 DDS_AllOrtho_Unripe_Noise <- DESeqSpline(Expt_All_Ortho_Unripe,
                                                                   vsNoise = TRUE,
                                                                   timeVar="Stage",
                                                                   CollapseTechRep = TRUE)
                                 colnames(DDS_AllOrtho_Unripe_Noise) <- NULL 
                                 saveRDS(DDS_AllOrtho_Unripe_Noise, "DEGAnalysis/RNA-seq/DDS_AllOrtho_Unripe_Noise.rds")
                                 return(DDS_AllOrtho_Unripe_Noise)})

# Test for DEGs with diff patterns by species in five taxon set
DDS_FiveOrtho_Unripe_DEGBySpecies <- tryCatch(readRDS("DEGAnalysis/RNA-seq/DDS_FiveOrtho_Unripe_DEGBySpecies.rds"),
                                              error=function(e){
                                                DDS_FiveOrtho_Unripe_DEGBySpecies <- DESeqSpline(Expt_Five_Ortho_Unripe,
                                                                                                 CaseCtlVar = "Species",
                                                                                                 timeVar = "Stage",
                                                                                                 CollapseTechRep = TRUE)
                                                # columns have duplicate names, which would be changed and mess up metadata mapping
                                                colnames(DDS_FiveOrtho_Unripe_DEGBySpecies) <- NULL 
                                                saveRDS(DDS_FiveOrtho_Unripe_DEGBySpecies, "DEGAnalysis/RNA-seq/DDS_FiveOrtho_Unripe_DEGBySpecies.rds")
                                                return(DDS_FiveOrtho_Unripe_DEGBySpecies)})

# Test for DEGs with diff patterns by fruit type
DDS_FiveOrtho_Unripe_DEGByFruit <- tryCatch(readRDS("DEGAnalysis/RNA-seq/DDS_FiveOrtho_Unripe_DEGByFruit.rds"),
                                            error=function(e){
                                              DDS_FiveOrtho_Unripe_DEGByFruit <- DESeqSpline(Expt_Five_Ortho_Unripe,
                                                                                             CaseCtlVar = "Fruit",
                                                                                             timeVar = "Stage",
                                                                                             CollapseTechRep = TRUE)
                                              # columns have duplicate names, which would be changed and mess up metadata mapping
                                              colnames(DDS_FiveOrtho_Unripe_DEGByFruit) <- NULL 
                                              saveRDS(DDS_FiveOrtho_Unripe_DEGByFruit, "DEGAnalysis/RNA-seq/DDS_FiveOrtho_Unripe_DEGByFruit.rds")
                                              return(DDS_FiveOrtho_Unripe_DEGByFruit)})

# Find the genes with a common pattern across species
DDS_FiveOrtho_Unripe_Noise <- tryCatch(readRDS("DEGAnalysis/RNA-seq/DDS_FiveOrtho_Unripe_Noise.rds"),
                                       error=function(e){
                                         DDS_FiveOrtho_Unripe_Noise <- DESeqSpline(Expt_Five_Ortho_Unripe,
                                                                                   vsNoise = TRUE,
                                                                                   timeVar="Stage",
                                                                                   CollapseTechRep = TRUE)
                                         colnames(DDS_FiveOrtho_Unripe_Noise) <- NULL 
                                         saveRDS(DDS_FiveOrtho_Unripe_Noise, "DEGAnalysis/RNA-seq/DDS_FiveOrtho_Unripe_Noise.rds")
                                         return(DDS_FiveOrtho_Unripe_Noise)})

#Sol species only
DDS_SolOrtho_Unripe_DEGBySpecies <- tryCatch(readRDS("DEGAnalysis/RNA-seq/DDS_SolOrtho_Unripe_DEGBySpecies.rds"),
                 error=function(e){
                   DDS_SolOrtho_Unripe_DEGBySpecies <- DESeqSpline(Expt_SolOrtho_Unripe,
                       CaseCtlVar = "Species",
                       timeVar = "Stage",
                       CollapseTechRep = TRUE)
                   # columns have duplicate names, which would be changed and mess up metadata mapping
                   colnames(DDS_SolOrtho_Unripe_DEGBySpecies) <- NULL 
                   saveRDS(DDS_SolOrtho_Unripe_DEGBySpecies, "DEGAnalysis/RNA-seq/DDS_SolOrtho_Unripe_DEGBySpecies.rds")
                   return(DDS_SolOrtho_Unripe_DEGBySpecies)})

# Test for DEGs with diff patterns by fruit type
DDS_SolOrtho_Unripe_DEGByFruit <- tryCatch(readRDS("DEGAnalysis/RNA-seq/DDS_SolOrtho_Unripe_DEGByFruit.rds"),
                                           error=function(e){
                                             DDS_SolOrtho_Unripe_DEGByFruit <- DESeqSpline(Expt_SolOrtho_Unripe,
                                                                                           CaseCtlVar = "Fruit",
                                                                                           timeVar = "Stage",
                                                                                           CollapseTechRep = TRUE)
                                             # columns have duplicate names, which would be changed and mess up metadata mapping
                                             colnames(DDS_SolOrtho_Unripe_DEGByFruit) <- NULL 
                                             saveRDS(DDS_SolOrtho_Unripe_DEGByFruit, "DEGAnalysis/RNA-seq/DDS_SolOrtho_Unripe_DEGByFruit.rds")
                                             return(DDS_SolOrtho_Unripe_DEGByFruit)})

# Find the genes with a common pattern across species
DDS_SolOrtho_Unripe_Noise <- tryCatch(readRDS("DEGAnalysis/RNA-seq/DDS_SolOrtho_Unripe_Noise.rds"),
                                      error=function(e){
                                        DDS_SolOrtho_Unripe_Noise <- DESeqSpline(Expt_SolOrtho_Unripe,
                                                                                 vsNoise = TRUE,
                                                                                 timeVar="Stage",
                                                                                 CollapseTechRep = TRUE)
                                        colnames(DDS_SolOrtho_Unripe_Noise) <- NULL 
                                        saveRDS(DDS_SolOrtho_Unripe_Noise, "DEGAnalysis/RNA-seq/DDS_SolOrtho_Unripe_Noise.rds")
                                        return(DDS_SolOrtho_Unripe_Noise)})



# Test for DEGs in Dry Fruited only
# DDS_DryOrtho <- tryCatch(readRDS("DEGAnalysis/RNA-seq/DDS_DryOrtho.rds"),
#                          error=function(e){
#                            DDS_DryOrtho <- DESeqSpline(Expt_Dry_Ortho,
#                                                        CaseCtlVar = "Species",
#                                                        timeVar = "Stage")
#                            # columns have duplicate names, which would be changed and mess up metadata mapping
#                            colnames(DDS_DryOrtho) <- NULL 
#                            saveRDS(DDS_DryOrtho, "DEGAnalysis/RNA-seq/DDS_DryOrtho.rds")
#                            return(DDS_DryOrtho)})

# Test for DE orthogenes at Ripening for Nobt only
DDS_NobtRipe_Ortho <- tryCatch(readRDS("DEGAnalysis/RNA-seq/DDS_NobtRipe_Ortho.rds"),
                               error=function(e){
                                 DDS_NobtRipe_Ortho <- DESeqDataSetFromMatrix(countData = assays(Expt_NobtRipe_Ortho)$counts,
                                                        colData = colData(Expt_NobtRipe_Ortho),
                                                        design = ~ Stage)
                                 # Collapse technical reps
                                 DDS_NobtRipe_Ortho$Accession <- sub("\\.\\d", "", DDS_NobtRipe_Ortho$Accession) #converted to char FYI
                                 DDS_NobtRipe_Ortho <- collapseReplicates(DDS_NobtRipe_Ortho,
                                                                DDS_NobtRipe_Ortho$Accession)
                                 # Carry On
                                 DDS_NobtRipe_Ortho <- estimateSizeFactors(DDS_NobtRipe_Ortho)
                                 DDS_NobtRipe_Ortho <- DESeq(DDS_NobtRipe_Ortho,
                                                             test="LRT",
                                                             reduced = ~1)
                                 saveRDS(DDS_NobtRipe_Ortho, "DEGAnalysis/RNA-seq/DDS_NobtRipe_Ortho.rds")
                                 return(DDS_NobtRipe_Ortho)})
# Test for DE orthogenes at Ripening for Slyc Only
DDS_Slyc1545_Ortho <- tryCatch(readRDS("DEGAnalysis/RNA-seq/DDS_Slyc1545_Ortho.rds"),
                               error=function(e){
                                 DDS_Slyc1545_Ortho <- DESeqDataSetFromMatrix(countData = assays(Expt_Slyc1545_Ortho)$counts,
                                                        colData = colData(Expt_Slyc1545_Ortho),
                                                        design = ~ Stage)
                                 DDS_Slyc1545_Ortho <- estimateSizeFactors(DDS_Slyc1545_Ortho)
                                 DDS_Slyc1545_Ortho <- DESeq(DDS_Slyc1545_Ortho,
                                                             test="LRT",
                                                             reduced = ~1)
                                 saveRDS(DDS_Slyc1545_Ortho, "DEGAnalysis/RNA-seq/DDS_Slyc1545_Ortho.rds")
                                 return(DDS_Slyc1545_Ortho)})
DDS_Slyc1535_Ortho <- tryCatch(readRDS("DEGAnalysis/RNA-seq/DDS_Slyc1535_Ortho.rds"),
                               error=function(e){
                                 DDS_Slyc1535_Ortho <- DESeqDataSetFromMatrix(countData = assays(Expt_Slyc1535_Ortho)$counts,
                                                                              colData = colData(Expt_Slyc1535_Ortho),
                                                                              design = ~ Stage)
                                 DDS_Slyc1535_Ortho <- estimateSizeFactors(DDS_Slyc1535_Ortho)
                                 DDS_Slyc1535_Ortho <- DESeq(DDS_Slyc1535_Ortho,
                                                             test="LRT",
                                                             reduced = ~1)
                                 saveRDS(DDS_Slyc1535_Ortho, "DEGAnalysis/RNA-seq/DDS_Slyc1535_Ortho.rds")
                                 return(DDS_Slyc1535_Ortho)})
DDS_Slyc3545_Ortho <- tryCatch(readRDS("DEGAnalysis/RNA-seq/DDS_Slyc3545_Ortho.rds"),
                               error=function(e){
                                 DDS_Slyc3545_Ortho <- DESeqDataSetFromMatrix(countData = assays(Expt_Slyc3545_Ortho)$counts,
                                                                              colData = colData(Expt_Slyc3545_Ortho),
                                                                              design = ~ Stage)
                                 DDS_Slyc3545_Ortho <- estimateSizeFactors(DDS_Slyc3545_Ortho)
                                 DDS_Slyc3545_Ortho <- DESeq(DDS_Slyc3545_Ortho,
                                                             test="LRT",
                                                             reduced = ~1)
                                 saveRDS(DDS_Slyc3545_Ortho, "DEGAnalysis/RNA-seq/DDS_Slyc3545_Ortho.rds")
                                 return(DDS_Slyc3545_Ortho)})
 
# Clustering --------------------------------------------------------------
# For all genes, this is best run noninteractively with the 4_NoninteractiveClustering.sh script
Cluster_Nobt <- tryCatch(readRDS("DEGAnalysis/RNA-seq/Cluster_Nobt.rds"),
                         error=function(e){
                           Cluster_Nobt <- DESeqCluster(DDS_NobtSE, numGenes = "all")
                           saveRDS(Cluster_Nobt, "DEGAnalysis/RNA-seq/Cluster_Nobt.rds")
                           return(Cluster_Nobt)})

Cluster_Melon <- tryCatch(readRDS("DEGAnalysis/RNA-seq/Cluster_Melon.rds"),
                         error=function(e){
                          Cluster_Melon <- DESeqCluster(DDS_Melon, numGenes = "all")
                          saveRDS(Cluster_Melon, "DEGAnalysis/RNA-seq/Cluster_Melon.rds")
                          return(Cluster_Melon)})


Cluster_Slyc <- tryCatch(readRDS("DEGAnalysis/RNA-seq/Cluster_Slyc.rds"),
                         error=function(e){
                           Cluster_Slyc <- DESeqCluster(DDS_Slyc, numGenes = "all")
                           saveRDS(Cluster_Slyc, "DEGAnalysis/RNA-seq/Cluster_Slyc.rds")
                           return(Cluster_Slyc)})

Cluster_SlycSRA <- tryCatch(readRDS("DEGAnalysis/RNA-seq/Cluster_SlycSRA.rds"),
                            error=function(e){
                              Cluster_SlycSRA <- DESeqCluster(DDS_SlycSRA, numGenes = "all")
                              saveRDS(Cluster_SlycSRA, "DEGAnalysis/RNA-seq/Cluster_SlycSRA.rds")
                              return(Cluster_SlycSRA)})

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

Cluster_TAIR <- tryCatch(readRDS("DEGAnalysis/RNA-seq/Cluster_TAIR.rds"),
                         error=function(e){
                           Cluster_TAIR <- DESeqCluster(DDS_TAIR, numGenes = "all")
                           saveRDS(Cluster_TAIR, "DEGAnalysis/RNA-seq/Cluster_TAIR.rds")
                           return(Cluster_TAIR)})


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
              Cluster_AllOrtho_Noise <- DESeqCluster(DDS_AllOrtho_Noise,
                                            numGenes = "all",
                                            timeVar = "Stage")
              saveRDS(Cluster_AllOrtho_Noise, "DEGAnalysis/RNA-seq/Cluster_AllOrtho_Noise.rds")
                                  return(Cluster_AllOrtho_Noise)})

Cluster_SolOrtho_DEGBySpecies <- tryCatch(readRDS("DEGAnalysis/RNA-seq/Cluster_SolOrtho_DEGBySpecies.rds"),
            error=function(e){
              Cluster_SolOrtho_DEGBySpecies <- DESeqCluster(DDS_SolOrtho_DEGBySpecies,
                                             numGenes = "all",
                                             CaseCtlVar = "Species",    
                                             timeVar = "Stage")
              saveRDS(Cluster_SolOrtho_DEGBySpecies, "DEGAnalysis/RNA-seq/Cluster_SolOrtho_DEGBySpecies.rds")
              return(Cluster_SolOrtho_DEGBySpecies)})

Cluster_SolOrtho_DEGByFruit <- tryCatch(readRDS("DEGAnalysis/RNA-seq/Cluster_SolOrtho_DEGByFruit.rds"),
           error=function(e){
             Cluster_SolOrtho_DEGByFruit <- DESeqCluster(DDS_SolOrtho_DEGByFruit,
                                             numGenes = "all",
                                             CaseCtlVar = "Fruit",
                                             timeVar = "Stage")
             saveRDS(Cluster_SolOrtho_DEGByFruit, "DEGAnalysis/RNA-seq/Cluster_SolOrtho_DEGByFruit.rds")
             return(Cluster_SolOrtho_DEGByFruit)})
Cluster_SolOrtho_Noise <- tryCatch(readRDS("DEGAnalysis/RNA-seq/Cluster_SolOrtho_Noise.rds"), 
           error=function(e){
             Cluster_SolOrtho_Noise <- DESeqCluster(DDS_SolOrtho_Noise,
                                                    numGenes = "all",
                                                    timeVar = "Stage")
             saveRDS(Cluster_SolOrtho_Noise, "DEGAnalysis/RNA-seq/Cluster_SolOrtho_Noise.rds")
             return(Cluster_SolOrtho_Noise)})

Cluster_AllOrtho_Unripe_DEGBySpecies <- tryCatch(readRDS("DEGAnalysis/RNA-seq/Cluster_AllOrtho_Unripe_DEGBySpecies.rds"),
                error=function(e){
                  Cluster_AllOrtho_Unripe_DEGBySpecies <- DESeqCluster(DDS_AllOrtho_Unripe_DEGBySpecies,
                        numGenes = "all",
                        CaseCtlVar = "Species",    
                        timeVar = "Stage")
                  saveRDS(Cluster_AllOrtho_Unripe_DEGBySpecies, "DEGAnalysis/RNA-seq/Cluster_AllOrtho_Unripe_DEGBySpecies.rds")
                  return(Cluster_AllOrtho_Unripe_DEGBySpecies)})

Cluster_AllOrtho_Unripe_DEGByFruit <- tryCatch(readRDS("DEGAnalysis/RNA-seq/Cluster_AllOrtho_Unripe_DEGByFruit.rds"),
                error=function(e){
                  Cluster_AllOrtho_Unripe_DEGByFruit <- DESeqCluster(DDS_AllOrtho_Unripe_DEGByFruit,
                    numGenes = "all",
                    CaseCtlVar = "Fruit",
                    timeVar = "Stage")
                  saveRDS(Cluster_AllOrtho_Unripe_DEGByFruit, "DEGAnalysis/RNA-seq/Cluster_AllOrtho_Unripe_DEGByFruit.rds")
                  return(Cluster_AllOrtho_Unripe_DEGByFruit)})

Cluster_AllOrtho_Unripe_Noise <- tryCatch(readRDS("DEGAnalysis/RNA-seq/Cluster_AllOrtho_Unripe_Noise.rds"), 
                                   error=function(e){
                                     Cluster_AllOrtho_Unripe_Noise <- DESeqCluster(DDS_AllOrtho_Unripe_Noise,
                                                                            numGenes = "all",
                                                                            timeVar = "Stage")
                                     saveRDS(Cluster_AllOrtho_Unripe_Noise, "DEGAnalysis/RNA-seq/Cluster_AllOrtho_Unripe_Noise.rds")
                                     return(Cluster_AllOrtho_Unripe_Noise)})

Cluster_SolOrtho_Unripe_DEGBySpecies <- tryCatch(readRDS("DEGAnalysis/RNA-seq/Cluster_SolOrtho_Unripe_DEGBySpecies.rds"),
                error=function(e){
                  Cluster_SolOrtho_Unripe_DEGBySpecies <- DESeqCluster(DDS_SolOrtho_Unripe_DEGBySpecies,
                                                                numGenes = "all",
                                                                CaseCtlVar = "Species",    
                                                                timeVar = "Stage")
                  saveRDS(Cluster_SolOrtho_Unripe_DEGBySpecies, "DEGAnalysis/RNA-seq/Cluster_SolOrtho_Unripe_DEGBySpecies.rds")
                  return(Cluster_SolOrtho_Unripe_DEGBySpecies)})

Cluster_SolOrtho_Unripe_DEGByFruit <- tryCatch(readRDS("DEGAnalysis/RNA-seq/Cluster_SolOrtho_Unripe_DEGByFruit.rds"),
                error=function(e){
                  Cluster_SolOrtho_Unripe_DEGByFruit <- DESeqCluster(DDS_SolOrtho_Unripe_DEGByFruit,
                                                              numGenes = "all",
                                                              CaseCtlVar = "Fruit",
                                                              timeVar = "Stage")
                  saveRDS(Cluster_SolOrtho_Unripe_DEGByFruit, "DEGAnalysis/RNA-seq/Cluster_SolOrtho_Unripe_DEGByFruit.rds")
                  return(Cluster_SolOrtho_Unripe_DEGByFruit)})
Cluster_SolOrtho_Unripe_Noise <- tryCatch(readRDS("DEGAnalysis/RNA-seq/Cluster_SolOrtho_Unripe_Noise.rds"), 
               error=function(e){
                 Cluster_SolOrtho_Unripe_Noise <- DESeqCluster(DDS_SolOrtho_Unripe_Noise,
                                                        numGenes = "all",
                                                        timeVar = "Stage")
                 saveRDS(Cluster_SolOrtho_Unripe_Noise, "DEGAnalysis/RNA-seq/Cluster_SolOrtho_Unripe_Noise.rds")
                 return(Cluster_SolOrtho_Unripe_Noise)})

Cluster_DryOrtho <- tryCatch(readRDS("DEGAnalysis/RNA-seq/Cluster_DryOrtho.rds"),
                             error=function(e){
                               Cluster_DryOrtho <- DESeqCluster(DDS_DryOrtho,
                                                                numGenes = "all",
                                                                CaseCtlVar = "Species",
                                                                timeVar = "Stage")
                               saveRDS(Cluster_DryOrtho, "DEGAnalysis/RNA-seq/Cluster_DryOrtho.rds")
                               return(Cluster_DryOrtho)})

# Play around with individual genes ---------------------------------------
Exampledds <- DDS_AllOrtho_DEGByFruit #assign one dds as the example to streamline code
ExampleRes <- results(Exampledds) #get results
ExampleResSig <- subset(ExampleRes, padj < 0.05) #subset by FDR
head(ExampleResSig[order(ExampleResSig$padj ), ]) #see best fitting genes for spline model
#Examine an individual Gene
topGene <- rownames(ExampleRes)[which.min(ExampleRes$padj)]
colData(Exampledds)$DAP <- as.factor(colData(Exampledds)$DAP)
colData(Exampledds)$Stage <- as.factor(colData(Exampledds)$Stage)
plotCounts(Exampledds, gene=topGene, intgroup=c("Species", "Stage"), normalized = T) #plot best fitting gene

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
FULgenes<-c(euFULI="OG0004494",
            euFULII="OG0007327")
# Plot FUL Genes
for (i in 1:length(FULgenes)) {
  pdf(file=paste0("DEGAnalysis/RNA-seq/Plots/Slyc_", names(FULgenes[i]), ".pdf"), # _SRA v _IH on Slyc
      width=6,
      height=4)
  plotCounts(Exampledds,
             gene=FULgenes[i],
             intgroup=c("Stage"), # Remove Genotype for nonSlycSRA
             main=paste0(names(FULgenes[i]), 
                         " p=",
                         formatC(ExampleRes$padj[rownames(ExampleRes)==FULgenes[i]], format = "e", digits = 1)),
             xlab="Stae",
             normalized=T,
             replace=T)
  dev.off()
}

# Plot Cluster Profiles ---------------------------------------------------
# Use the ggplot plotting scripts for figures. Use this for exploring
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

# Save Cluster genes to a file --------------------------------------------
X <- split(Cluster_AllOrtho_Noise$df, Cluster_AllOrtho_Noise$df$cluster)
for (i in 1:length(X)) {
  write.table(row.names(X[[i]]), 
              file=paste0("DEGAnalysis/RNA-seq/Lists/AllOrtho_Noise_Cluster_", max(X[[i]]$cluster), ".txt"),
              row.names = FALSE,
              quote = FALSE,
              col.names = FALSE)
}

