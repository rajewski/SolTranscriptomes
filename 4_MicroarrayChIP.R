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
  #mcols(SL4gff) <- mcols(SL4gff)[,c(2,7)]
  #names(mcols(SL4gff)) <- c("type", "name")
  SL4gff.df <- as.data.frame(SL4gff)
  SL4gff.df <- SL4gff.df[,c("seqnames", "start", "end", "strand", "type", "Name")]
  SL4gff.df <- SL4gff.df[SL4gff.df$type=="mRNA",]
  colnames(SL4gff.df) <- c("chr", "start", "end", "strand", "type", "name")
  saveRDS(SL4gff.df, file="ChIPAnalysis/ChIP-chip/SL4gff.rds")
})

# Read in Raw Data -----------------------------------------------
# All Data
# The 3rd rep for each is a dye swap, so I simply flipped the file names.
# I could have written in the actual swap on the cy3 and cy5 channels, but that might complicate things
tryCatch(TomRG <-readRDS("ChIPAnalysis/ChIP-chip/GSE49125_Reduced.RG.rds"), error=function(e){
  # Ive decided to remove the second rep of FUL1 bc it was very odd and had low binding
  TomRG <- readNimblegen("ChIPAnalysis/ChIP-chip/TomDesignReduced.txt","/rhome/arajewski/R/x86_64-pc-linux-gnu-library/3.6/Ringo/exData/spottypes.txt", path=NULL)
  saveRDS(TomRG, file="ChIPAnalysis/ChIP-chip/GSE49125_Reduced.RG.rds")
})

#check for autocorrelation
Autocorr <- autocor(TomRG, probeAnno=SL4Remapping, chrom="SL4.0ch09", lag.max=1000)
plot(Autocorr)

#preprocess the intensities
TomRG <- backgroundCorrect(TomRG)
ChipEset <- preprocess(TomRG, method="vsn")
sampleNames(ChipEset) <- with(TomRG$targets, paste(Cy5,"vs",Cy3, Rep,sep="_"))

#smooth the peaks in sliding windows
ChipEsetSmooth <- computeRunningMedians(ChipEset, 
                                        probeAnno=SL4Remapping,
                                        modColumn = "Cy5",
                                        winHalfSize = 200,
                                        combineReplicates=TRUE,
                                        min.probes=3)

#try a plot for TAGL1
SL4gff.df[grep("Solyc07g055920", SL4gff.df$name),] #TAGL1 then add ~4kb upstream for promoter
SL4gff.df[grep("Solyc02g077920", SL4gff.df$name),] #CNR then add ~4kb upstream for promoter
SL4gff.df[grep("Solyc02g086930", SL4gff.df$name),] #HB-1 then add ~4kb upstream for promoter

chipAlongChrom1(ChipEsetSmooth,
                SL4Remapping, 
                gff=SL4gff.df,
                ylim=c(0,3),
                chrom="SL4.0ch07", #TAGL1
                xlim=c(63756000, 63760000)) #TAGL1
                #chrom="SL4.0ch02", #CNR or HB-1
                #xlim=c(40730780, 40733000)) #CNR
                #xlim=c(47530500, 47533000))


#make a hisotgram of reporter intensities to get a threshold
(y0 <- apply(exprs(ChipEsetSmooth),2,upperBoundNull))
hist(exprs(ChipEsetSmooth)[,2])
(y0G <- apply(exprs(ChipEsetSmooth),2,twoGaussiansNull))


#actually find chip enriched regions
chersX <- findChersOnSmoothed(ChipEsetSmooth,
                              probeAnno=SL4Remapping,
                              threshold=y0G,
                              distCutOff=200)
chersX <- relateChers(chersX, gff=SL4gff.df, upstream=3000)
chersXD <- as.data.frame(chersX)
head(chersXD[order(chersXD$maxLevel, decreasing=TRUE),])

#save it as a bed object for IGV
chersToBED(chersX, file="ChIPAnalysis/ChIP-chip/FULchers.bed")

#Find distance from cher to TSS
library(tidyr)
Dist2TSS <- data.frame(.id=NA, V1=NA, antibody=NA, maxlevel=NA, score=NA)
for (i in 1:length(chersX)){
  if (length(chersX[[i]]@extras$distMid2TSS)>0) {
    tmpcher <- cbind(plyr::ldply(chersX[[i]]@extras$distMid2TSS),
                     antibody=chersX[[i]]@antibody,
                     maxlevel=chersX[[i]]@maxLevel,
                     score=chersX[[i]]@score)
    Dist2TSS <- rbind(Dist2TSS, tmpcher)
  } 
}
names(Dist2TSS) <- c("Transcript", "Distance", "Antibody", "Max_Level", "Score")
Dist2TSS <- Dist2TSS[-1,]
Dist2TSS <- Dist2TSS[Dist2TSS$Distance<=2000,] #array was only designed w/in 2kb of TSS
Dist2TSS$Bins <- paste0(strsplit(Dist2TSS$Antibody, ".sm"),
                        "_",
                        cut(Dist2TSS$Distance, 
                            breaks=c(-2,50,100,200,500,1000,1500,2000),
                            labels=c("50","100","200","500","1000","1500","2000")))

write.csv(Dist2TSS,
          "ChIPAnalysis/ChIP-chip/TomatoFULBinding.tsv",
          row.names = F,
          quote = F)

pdf(file="ChIPAnalysis/ChIP-chip/Plot_FUL1_Dist2TSS.pdf")
  hist(Dist2TSS$Distance[grep("FUL1",Dist2TSS$Antibody)], 
       breaks=20,
       xlab="Distance from TSS (bp)",
       main=expression(paste(alpha,"-FUL1 Binding")))
dev.off()

pdf(file="ChIPAnalysis/ChIP-chip/Plot_FUL2_Dist2TSS.pdf")
  hist(Dist2TSS$Distance[grep("FUL2",Dist2TSS$Antibody)], 
       breaks=20,
       xlab="Distance from TSS (bp)",
       main=expression(paste(alpha,"-FUL2 Binding")))
dev.off()
