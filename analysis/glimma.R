library(limma)
library(edgeR)
so<-readRDS("Sleuth_so.rds")
library(tidyr)
data <- so$obs_raw[,c(1,2,6)] %>% spread(sample,est_counts)
rownames(data) <- data[,1] #give row names
data <- data[,-1] #drop names
data <- data[,c(4:9,1:3,10:15)] #reorder columns by time

#make list of genes
dge1 <- DGEList(counts=data, group=c(1,1,1,3,3,3,15,15,15,32,32,32,40,40,40))
dge1 <- dge1[filterByExpr(dge1),] #filter out low-expressed genes
dge1 <- calcNormFactors(dge1) 
countDF1 <- voom(dge1, normalize.method = "quantile") 

##PCA
PCA<-as.data.frame(t(countDF1$E))
SAMP=rep(c(1,2,3),5)
TREAT=c(1,1,1,3,3,3,15,15,15,32,32,32,40,40,40)
Genotype <- factor(SAMP)
Treatment <- factor(TREAT) 
library(devtools)
install_github("vqv/ggbiplot")
library(ggbiplot)
ir.Genotype <- Genotype
ir.Treatment <-Treatment
r<-PCA[,apply(countDF1, 2, var, na.rm=TRUE) != 0] 
ir.pca <- prcomp(r,center = T,scale. = F, tol=0) 
#print PCA
print(ggbiplot(ir.pca, obs.scale = 1, var.scale = 1, groups = ir.Genotype, labels=rownames(ir.pca$x), 
               labels.size=0, var.axes=0)
      +geom_point(aes(shape=ir.Treatment),size=3)
      +geom_point(aes(colour=ir.Genotype),size=2))

##Get DEGs and make clustered heatmap
#try a natural spline model
library(splines)
full_design <- model.matrix(formula(~ ns(c(1,1,1,3,3,3,15,15,15,32,32,32,40,40,40), df = 3)))
colnames(full_design) <- c("A", "B", "C", "D")
fit <- lmFit(countDF1$E, full_design) #this may actually need log counts
fit <- eBayes(fit)
results <- decideTests(fit, p.value=0.05,method='separate', adjust.method='BH')
#topTable(fit)
#summary(fit)
#summary(results)
write.fit(fit,'Sleuth_splines.txt', results=results, adjust='BH', digits=45)

#kludge code that works??
data <- read.delim("Sleuth_splines.txt",check.names=FALSE) 
####Select columns with gene names, fold change (LFC or Coef) and FDR scores (FDR or p.value.adj.)
data<-data[,c(2:5,19)]
#filter by F.p.value
data <- data["F.p.value">0.05,]
#calculate distance on LFC values and cluster
library(cluster)
no.clust <- 15 #arbitrary here
euc <- dist(data[,1:4], method = "euclidean")
pam <-pam(data, no.clust, diss = inherits(data, "dist"), metric = "euclidean", medoids = NULL, stand = FALSE, cluster.only = FALSE, do.swap = TRUE, pamonce = FALSE, trace.lev = 0)


rownames(DEG)



