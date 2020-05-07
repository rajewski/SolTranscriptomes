library("limma")
#library("affy")
#library("makecdfenv")
`%notin%` <- Negate(`%in%`)

# Arabidopsis Affymetrix DEX Expt -----------------------------------------
# I have analyzed this data and find that FUL is not differentially expressed by the 
# DEX induction. I have tried both the CDF from NCBI and one from MBNI, but neither makes a difference.
# Incorporation of array weights does not help increase the signal to noise ratio. Switching from 
# processing the CEL files myself to using the non-log-2 normalized data from NCBI also does not change
# the conclusion. I suspect that this is why the data was never publish, but then why was it uplaoded?
#Make an environment with the probe mappings from NCBI
GPL14926Env <- make.cdf.env(filename = "GPL14926_agronomics1_TAIR9_gene.CDF",
                            cdf.path = "ExternalData/Microarray/")
#try another array CDF from http://mbni.org/customcdf/18.0.0/tairg.download/AGRONOMICS1_At_TAIRG_18.0.0.zip
#AGRONOMICS1 <- make.cdf.env(filename = "AGRONOMICS1_At_TAIRG.cdf", cdf.path = "ExternalData/Microarray/")

#Read in the experimental design
DexTargets <- readTargets("DEGAnalysis/Microarray/DexDesign.tsv") 
rownames(DexTargets) <- basename(DexTargets$FileName)
DexTargets$Treatment <- factor(DexTargets$Treatment, levels = c("Mock", "DEX"))
DexTargets$Genotype <- factor(DexTargets$Genotype, levels = c("ful-1", "ful-1_FUL-GR"))

#Read in and normalize the data
DexEset <- ReadAffy(filenames = DexTargets$FileName,
                    cdfname="GPL14926Env",
                    phenoData=DexTargets[,2:4])
DexRMA <- affy::rma(DexEset)

#ust some checks on data quality
head(exprs(DexRMA))
plotMDS(DexRMA) 

#design the model matrix https://support.bioconductor.org/p/108827/
DexDesign <- model.matrix(~DexRMA@phenoData@data$Genotype * DexRMA@phenoData@data$Treatment)
colnames(DexDesign) <- c("Intercept", "Genotype", "Dex", "Ixn")

#Incorporate array weights
DexArrayw <- arrayWeights(DexRMA)
barplot(DexArrayw, xlab="Array", ylab="Weight", col="white", las=2)
abline(h=1, lwd=1, lty=2)

#Fit the model
Dexfit <- lmFit(DexRMA, DexDesign, weights = DexArrayw)
Dexfit <- eBayes(Dexfit, trend=TRUE, robust=TRUE)

#Show the top genes by p-value when only looking at the interaction of Dex and Genotype (Coef=4)
topTable(Dexfit, coef=4)
# Look into results by a couple of positive control genes FUL and SHP1
# AT5G60910.X.Chr5.minus.24502482.24506143 #FUL
# AT3G58780.X.Chr3.plus.21738460.21741907 #SHP1
topTable(Dexfit, coef=4, number=Inf)["AT5G60910.X.Chr5.minus.24502482.24506143",]
#Set a significance threshold based on the p-value of non-DEG control genes prefixed with AFFX
Dexi <- grep("AFFX",rownames(DexRMA))
DexResults <- decideTests(Dexfit, method="nestedF", p.value = (min(Dexfit$p.value[Dexi, "Ixn"])))
summary(DexResults)
DexSigRes <- classifyTestsF(Dexfit, p.value = (min(Dexfit$p.value[Dexi, "Ixn"])))
summary(DexsigRes)
DexsigRes["AT5G60910.X.Chr5.minus.24502482.24506143",]
DexsigRes["AT3G58780.X.Chr3.plus.21738460.21741907",]

#Try with the normalized non-log2 data from NCBI
DPM <- read.csv("ExternalData/Microarray/GSE79553_Normalized_data_with_all_controls.txt",
                       header=TRUE,
                       sep="\t")
DPM <- DPM[!duplicated.array(DPM),]
rownames(DPM) <- DPM[,1]
DPM <- DPM[,-1]
DPMw <- arrayWeights(DPM)
DPMfit <- lmFit(DPM, DexDesign, weights=DPMw)
DPMfit <- eBayes(DPMfit)
topTable(DPMfit, coef=4)
DPMi <- grep("AFFX",rownames(DPM))
DPMSigRes <- classifyTestsF(DPMfit, p.value = (min(DPMfit$p.value[DPMi, "Ixn"])))
summary(DPMSigRes)
topTable(DPMfit, coef=4, number=Inf)["AT5G60910.X.Chr5.minus.24502482.24506143",]
DPMSigRes["AT5G60910.X.Chr5.minus.24502482.24506143",]
DPMSigRes["AT3G58780.X.Chr3.plus.21738460.21741907",]




# Tomato Agilent Timecourse -----------------------------------------------
# Import Expt design
TomTargets <- read.table("DEGAnalysis/Microarray/TomatoDesign.tsv", header=T)
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

