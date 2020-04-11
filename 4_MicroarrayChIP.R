library("Ringo")
# Import Annotation Information --------------------------------------------------
# Make a probe mapping environment from the exonerate data
tryCatch(SL4Remapping <- readRDS("ChIPAnalysis/ChIP-chip/SL4MRemapping.rds"), error=function(e){
  TomTargetMatches <- read.table("ExternalData/Microarray/GSE49125/GPL15968-24228.final.txt")
  colnames(TomTargetMatches) <- c("ProbeID", "Chromosome", "StartPos", "StopPos", "Strand")
  TomTargetMatches$Length <- abs(TomTargetMatches$StopPos - TomTargetMatches$StartPos)
  TomTargetMatches$Position <- pmin(TomTargetMatches$StartPos,TomTargetMatches$StopPos) #probes have no inherent direction
  TomTargetMatches$ProbeID <- as.character(TomTargetMatches$ProbeID)
  TomTargetMatches$Chromosome <- as.character(TomTargetMatches$Chromosome)
  SL4Remapping <- posToProbeAnno(pos = TomTargetMatches,
                                 chrNameColumn = "Chromosome",
                                 probeColumn = "ProbeID",
                                 lengthColumn = "Length",
                                 chrPositionColumn = "Position",
                                 genome="SL4.0",
                                 microarrayPlatform = "GPL15968",
                                 verbose=TRUE)
})

# Import the GFF
tryCatch(SL4gff.df <- readRDS("ChIPAnalysis/ChIP-chip/SL4gff.rds"), error=function(e){
  library(GenomicRanges)
  library(rtracklayer)
  SL4gff <- import("SlycDNA/ITAG4.0_gene_models.gff")
  SL4gff.df <- as.data.frame(SL4gff)
  SL4gff.df <- SL4gff.df[,c("seqnames", "start", "end", "strand", "type", "Name")]
  SL4gff.df <- SL4gff.df[SL4gff.df$type=="gene",]
  colnames(SL4gff.df) <- c("chr", "start", "end", "strand", "type", "name")
  saveRDS(SL4gff.df, file="ChIPAnalysis/ChIP-chip/SL4gff.rds")
})

# Read in Raw Data -----------------------------------------------
# All Data
# The 3rd rep for each is a dye swap, so I simply flipped the file names.
# I could have written in the actual swap on the cy3 and cy5 channels, but that might complicate things
tryCatch(TomRG <-readRDS("ChIPAnalysis/ChIP-chip/GSE49125_RG.rds"), error=function(e){
  TomRG <- readNimblegen("ChIPAnalysis/ChIP-chip/TomDesign.txt","/rhome/arajewski/R/x86_64-pc-linux-gnu-library/3.6/Ringo/exData/spottypes.txt", path=NULL)
  saveRDS(TomRG, file="ChIPAnalysis/ChIP-chip/GSE49125_RG.rds")
})

# Just FUL1
tryCatch(FUL1RG <- readRDS("ChIPAnalysis/ChIP-chip/GSE49125_FUL1_RG.rds"), error=function(e){
  FUL1RG <- readNimblegen("ChIPAnalysis/ChIP-chip/FUL1Design.txt", "/rhome/arajewski/R/x86_64-pc-linux-gnu-library/3.6/Ringo/exData/spottypes.txt", path=NULL)
  saveRDS(FUL1RG, file="ChIPAnalysis/ChIP-chip/GSE49125_FUL1_RG.rds")
})

# Assign one dataset to this variable
RG <- TomRG

#check for autocorrelation
Autocorr <- autocor(RG, probeAnno=SL4Remapping, chrom="SL4.0ch09", lag.max=1000)
plot(Autocorr)

#preprocess the intensities
ChipEset <- preprocess(RG, method="nimblegen")
sampleNames(ChipEset) <- with(RG$targets, paste(Cy5,"vs",Cy3, Rep,sep="_"))

#smooth the peaks in 800bp windows
ChipEsetSmooth <- computeRunningMedians(ChipEset, 
                                        probeAnno=SL4Remapping,
                                        modColumn = "Cy5",
                                        winHalfSize = 300,
                                        combineReplicates=TRUE)
sampleNames(ChipEsetSmooth) <- paste(sampleNames(ChipEset),"smoothed")


#try a plot for TAGL1
SL4gff.df[grep("Solyc07g055920", SL4gff.df$name),] #TAGL1 then add ~4kb upstream for promoter
SL4gff.df[grep("Solyc02g077920", SL4gff.df$name),] #CNR then add ~4kb upstream for promoter
SL4gff.df[grep("Solyc02g086930", SL4gff.df$name),] #HB-1 then add ~4kb upstream for promoter

chipAlongChrom1(ChipEset,
                SL4Remapping, 
                gff=SL4gff.df,
                ylim=c(0,2),
                #chrom="SL4.0ch07", #TAGL1
                #xlim=c(63756000, 63760000)) #TAGL1
                chrom="SL4.0ch02", #CNR or HB-1
                #xlim=c(40730780, 40733000)) #CNR
                xlim=c(47530500, 47533000))


#make a hisotgram of reporter intensities to get a threshold
(y0 <- apply(exprs(ChipEsetSmooth),2,upperBoundNull))

chersX <- findChersOnSmoothed(ChipEsetSmooth, probeAnno=SL4Remapping, threshold=y0)
chersX <- relateChers(chersX, gff=SL4gff.df, upstream=3000)
chersXD <- as.data.frame(chersX)
head(chersXD[order(chersXD$maxLevel, decreasing=TRUE),])
#save it as a bed object for IGV
chersToBED(chersX, file="ChIPAnalysis/ChIP-chip/chers.bed")

