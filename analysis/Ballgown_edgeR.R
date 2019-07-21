#This is an analysis script for the data that is output by Ballgown (mapped by stringtie).
#This is not to be confused with the analysis done by Kallisto, which is not my current method

#the methods for this are taken from https://bioconductor.org/packages/release/bioc/vignettes/ballgown/inst/doc/ballgown.html
#with some adjustments from https://rpubs.com/kapeelc12/Ballgown

#BiocManager::install("ballgown")
library(ballgown) #v 2.16.0

#Read in Ballgown experimental data  tables from StringTie
bg <- readRDS("bg.rds")
#bg <- ballgown(dataDir = "~/bigdata/Slycopersicum/Results_Stringtie/Ballgown/",samplePattern = "AC*", meas="all")
#pData(bg) = data.frame(id=sampleNames(bg), condition=c(rep("15DPA",3), rep("1DPA", 3)) ,rep=rep(c(1,2,3),2))
#saveRDS(bg, file="bg.rds")

#filter out lowly expressed genes by stdev
library(genefilter)
bg_filt = subset(bg,"rowVars(texpr(bg)) >1",genomesubset=TRUE)

#Confirm the libraries are normal
gene_expression <- bg_filt@expr$trans[,c(10,12,14,16,18,20,22)]
boxplot(log2(gene_expression[,2:7]+1), 
        names=sampleNames(bg_filt), las=2, ylab="log2(FPKM)", 
        main="Distribution of FPKMs for all 6 libraries")

#make an MDS to show clustering of libraries
gene_expression[,"sum"]=apply(gene_expression[,2:7], 1, sum)
r=cor(gene_expression[gene_expression$sum>5,2:7], use="pairwise.complete.obs", method="pearson")
mds=cmdscale(1-r, k=2, eig=TRUE)
plot(mds$points, type="n", xlab="", ylab="", main="MDS distance plot (all non-zero genes) for all libraries")
points(mds$points[,1], mds$points[,2], col="grey", cex=2, pch=16)
text(mds$points[,1], mds$points[,2], sampleNames(bg), col=c("tomato1", "tomato1", "tomato1", "green", "green", "green"))

#Calculate differential gene expression
stattest_results <- stattest(bg_filt, feature="gene", meas="FPKM", covariate="condition", getFC = T)


#do I even need the replist.tst file?