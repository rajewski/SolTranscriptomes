
# Arabidopsis Affymetrix DEX Expt -----------------------------------------
# I have analyzed this data and find that FUL is not differentially expressed by the 
# DEX induction. I have tried both the CDF from NCBI and one from MBNI, but neither makes a difference.
# Incorporation of array weights does not help increase the signal to noise ratio. Switching from 
# processing the CEL files myself to using the non-log-2 normalized data from NCBI also does not change
# the conclusion. I suspect that this is why the data was never publish, but then why was it uplaoded?
library("limma")
library("affy")
library("makecdfenv")
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
design <- model.matrix(~DexRMA@phenoData@data$Genotype * DexRMA@phenoData@data$Treatment)
colnames(design) <- c("Intercept", "Genotype", "Dex", "Ixn")

#Incorporate array weights
arrayw <- arrayWeights(DexRMA)
barplot(arrayw, xlab="Array", ylab="Weight", col="white", las=2)
abline(h=1, lwd=1, lty=2)

#Fit the model
fit <- lmFit(DexRMA, design, weights = arrayw)
fit <- eBayes(fit, trend=TRUE, robust=TRUE)

#Show the top genes by p-value when only looking at the interaction of Dex and Genotype (Coef=4)
topTable(fit, coef=4)
# Look into results by a couple of positive control genes FUL and SHP1
# AT5G60910.X.Chr5.minus.24502482.24506143 #FUL
# AT3G58780.X.Chr3.plus.21738460.21741907 #SHP1
topTable(fit, coef=4, number=Inf)["AT5G60910.X.Chr5.minus.24502482.24506143",]
#Set a significance threshold based on the p-value of non-DEG control genes prefixed with AFFX
i <- grep("AFFX",rownames(DexRMA))
results <- decideTests(fit, method="nestedF", p.value = (min(fit$p.value[i, "Ixn"])))
summary(results)
sigRes <- classifyTestsF(fit, p.value = (min(fit$p.value[i, "Ixn"])))
summary(sigRes)
sigRes["AT5G60910.X.Chr5.minus.24502482.24506143",]
sigRes["AT3G58780.X.Chr3.plus.21738460.21741907",]

#Try with the normalized non-log2 data from NCBI
DexPreNorm <- read.csv("ExternalData/Microarray/GSE79553_Normalized_data_with_all_controls.txt",
                       header=TRUE,
                       sep="\t")
DexPreNorm <- DexPreNorm[!duplicated.array(DexPreNorm),]
rownames(DexPreNorm) <- DexPreNorm[,1]
DexPreNorm <- DexPreNorm[,-1]
DexPreNormw <- arrayWeights(DexPreNorm)
fitDexPreNorm <- lmFit(DexPreNorm, design, weights=DexPreNormw)
fitDexPreNorm <- eBayes(fitDexPreNorm)
topTable(fitDexPreNorm, coef=4)
iDPM <- grep("AFFX",rownames(DexPreNorm))
sigResDPM <- classifyTestsF(fitDexPreNorm, p.value = (min(fitDexPreNorm$p.value[iDPM, "Ixn"])))
summary(sigResDPM)
topTable(fitDexPreNorm, coef=4, number=Inf)["AT5G60910.X.Chr5.minus.24502482.24506143",]
sigResDPM["AT5G60910.X.Chr5.minus.24502482.24506143",]
sigResDPM["AT3G58780.X.Chr3.plus.21738460.21741907",]

