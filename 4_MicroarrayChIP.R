library("Ringo")

# Tomato ChIP-chip -----------------------------------------------
# Read in Raw Data
tryCatch(loadRDS("DEGAnalysis/Microarray/GSE49125_RG.rds"), error=function(e){
  TomRG <- readNimblegen("ChIPAnalysis/ChIP-chip/TomDesign.txt","/rhome/arajewski/R/x86_64-pc-linux-gnu-library/3.6/Ringo/exData/spottypes.txt", path=NULL)
  saveRDS(TomRG, file="DEGAnalysis/Microarray/GSE49125_RG.rds")
})

#Make a probe mapping environment from the exonerate data
TomTargetMatches <- read.table("ExternalData/Microarray/GSE49125/GPL15968-24228.final.txt")
colnames(TomTargetMatches) <- c("ProbeID", "Chromosome", "StartPos", "StopPos", "Strand")
TomTargetMatches$Length <- abs(TomTargetMatches$StopPos - TomTargetMatches$StartPos)
TomTargetMatches$Position <- pmin(TomTargetMatches$StartPos,TomTargetMatches$StopPos)
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

TomChipEset <- preprocess(TomRG, method = "nimblegen")

library(GenomicRanges)
library(rtracklayer) #to import GFF
SL4gff <- import("SlycDNA/ITAG4.0_gene_models.gff")
mcols(SL4gff)









system.file("exData",package="Ringo")
TomTargets$Genotype <- factor(TomTargets$Genotype, levels=c("WT", "KD"))
TomTargets$Line <- factor(TomTargets$Line, levels=c("10","30"))

# Read in the data files
TomEset <- read.maimages(TomTargets$FileName,
                         source="agilent",
                         green.only = TRUE,
                         names = paste(TomTargets$Genotype,TomTargets$Line,TomTargets$Rep, sep="_"),
                         other.columns = "gIsWellAboveBG")

# Background correct and normalize between arrays
TomNorm <- backgroundCorrect(TomEset, method="normexp")
TomNorm <- normalizeBetweenArrays(TomNorm, method="quantile")

#start removing bad probes
TomControl <- TomNorm$genes$ControlType==1L
TomIsExpr <- rowSums(TomNorm$other$gIsWellAboveBG > 0) >= 1
TomGeneNames <- read.table("ExternalData/Microarray/GSE41560/Agilent022270.final.txt", stringsAsFactors = F)
TomBadMap <- TomNorm$genes$ProbeName %notin% TomGeneNames$V1
TomNormFilt <- TomNorm[!TomControl & !TomBadMap & TomIsExpr, ]

#Add Gene Names from exonerate output and remove useless columns
TomNormFilt$genes$GeneName <- TomGeneNames$V2[match(TomNormFilt$genes$ProbeName, TomGeneNames$V1)]
TomNormFilt$genes<- TomNormFilt$genes[,c("ProbeName", "GeneName")]

# Fit Model and test DEs
TomDesign <- model.matrix(~TomTargets$Genotype+TomTargets$Line)
TomWeights <- arrayWeights(TomNormFilt, design=TomDesign)
TomFit <- lmFit(TomNormFilt, TomDesign, weights=TomWeights)
TomFit <- eBayes(TomFit,trend=TRUE,robust=TRUE)
summary(decideTests(TomFit[,-1], p.value = 0.05))
topTable(TomFit, coef = 2)
#Are FUL1 and FUL2 *actually* knocked down?
topTable(TomFit, coef=2, number=Inf)[topTable(TomFit, coef=2, number=Inf)$GeneName=="Solyc06g069430.3.1",]
topTable(TomFit, coef=2, number=Inf)[topTable(TomFit, coef=2, number=Inf)$GeneName=="Solyc03g114830.3.1",] #there are 4 probes that match this
#There is no probe specific to MBP10, but what about MBP20?
topTable(TomFit, coef=2, number=Inf)[topTable(TomFit, coef=2, number=Inf)$GeneName=="Solyc02g089210.4.1",]
#Output a list of DEGs
write.csv(topTable(TomFit, coef = 2, p.value = 0.05, number=Inf), file = "DEGAnalysis/Microarray/GSE41560_DEGs.csv")

