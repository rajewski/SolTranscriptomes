

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
NobtBamfiles <- BamFileList(metadata$Path[metadata$Species=="Tobacco"], yieldSize=2000000)
TAIR10BamFiles <- BamFileList(metadata$Path[metadata$Species=="Arabidopsis"], yieldSize=2000000)
SlycSRABamFiles <- BamFileList(metadata$Path[metadata$Species=="Tomato" & metadata$PE==0], yieldSize=2000000)
SlycIHBamFiles <- BamFileList(metadata$Path[metadata$Species=="Tomato" & metadata$PE==0], yieldSize=2000000)
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

library("GenomicAlignments")
#To my knowledge the SRA experiments were not strand-specific
TAIR10Expt <- summarizeOverlaps(features=TAIR10genes,
                                reads=TAIR10BamFiles,
                                mode="Union",
                                singleEnd=TRUE,
                                ignore.strand=TRUE)
SlycSRAExpt <- summarizeOverlaps(features=Slycgenes,
                                 reads=SlycSRABamFiles,
                                 mode="Union",
                                 singleEnd=TRUE,
                                 ignore.strand=TRUE)
SlycIHExpt <- summarizeOverlaps(feature=Slycgenes,
                                reads=SlycIHBamFiles,
                                mode="Union",
                                singleEnd=FALSE,
                                ignore.strand=FALSE,
                                fragments=TRUE)





