library(limma)
library(edgeR)
library(tidyr)
so<-readRDS("Sleuth_so.rds")
data <- so$obs_norm[,c(1,3,2)] %>% spread(sample,est_counts)
rownames(data) <- data[,1] #give row names
data <- data[,-1] #drop names
data <- data[,c(4:9,1:3,10:15)] #reorder columns by time

#make list of genes
dge1 <- DGEList(counts=data, group=c(1,1,1,3,3,3,15,15,15,32,32,32,40,40,40))
dge1 <- dge1[filterByExpr(dge1),] #filter out low-expressed genes
dge1 <- calcNormFactors(dge1) 
countDF1 <- voom(dge1, normalize.method = "quantile") 

##PCA
#library(devtools)
#install_github("vqv/ggbiplot")
library(ggbiplot)
PCA<-as.data.frame(t(countDF1$E))
PCA<-PCA[,apply(countDF1, 2, var, na.rm=TRUE) != 0] 
ir.Genotype=factor(rep(c(1,2,3),5))
ir.Treatment=factor(c(1,1,1,3,3,3,15,15,15,32,32,32,40,40,40))
ir.pca <- prcomp(PCA,center = T,scale. = F, tol=0) 
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
fit <- lmFit(countDF1$E, full_design)
fit <- eBayes(fit)
results <- decideTests(fit, p.value=0.05,method='separate', adjust.method='BH')
#topTable(fit)
#summary(fit)
#summary(results)
write.fit(fit,'Sleuth_splines.txt', results=results, adjust='BH', digits=45)

##Clustering
library(cluster)
library(factoextra) #for fviz, remove if not necessary
library(e1071) # for cmeans
library(WGCNA)
#Read in the data
data <- read.delim("Sleuth_splines.txt",check.names=FALSE) 
data<-data[,c(2:5,19)] #Select gene names, LFC, F.p.value
data <- data[data$F.p.value<0.05,] #filter by F.p.value

#trying out several methods to find the optimum number of clusters (clara, cmeans, dbscan, and eventually WGCNA)
#find the optimum number of clusters
#so far each of these says 2-3 clusters, which is absurd. WGCNA says 20 but then finds 6, so im going with that for right now....
fviz_nbclust(data[,1:4], cmeans, method = "silhouette", k.max=20)+
  theme_classic()
gap_stat <- clusGap(data[,1:4], FUN=clara, K.max=10, B=10)
fviz_gap_stat(gap_stat)
no.clust <- 6 #arbitrary here, but see comment re WGCNA above

#Do clustering with clara, could sub in WGCNA here
clara <-clara(data, no.clust, samples=50, sampsize=200, metric = "euclidean", pamLike=TRUE)
saveRDS(clara, file="clara.RDS")
write.table(cbind(Cluster=clara$clustering,data[,-5]), paste0("DEGs_filtered_clara_k",no.clust,".txt"), quote=FALSE, sep="\t", col.names = NA)

#Summarize clusters
mat_data <- data.matrix(read.delim(paste0("DEGs_filtered_clara_k",no.clust,".txt"), row.names=1))
colnames(mat_data) <- c("Cluster", "FC.1.3", "FC.3.15", "FC.15.br", "FC.br.rr")
#find media pattern for each cluster
clust.med=cbind(1:no.clust, t(sapply(X=1:no.clust,FUN=function(i) {apply(mat_data[which(mat_data[,"Cluster"]==i),2:5], MARGIN=2,FUN=median)}))) #change apply index w/ diff timepoints?
clust.med=clust.med[order(clust.med[,2]),]  # sorts the clusters by the 1st timepoint
write.table(clust.med,paste0("DEGs_filtered_clara_k",no.clust,"_medians.txt",sep=""),sep="\t")
#resorts adds a column with the resorted cluster number
data.ord=matrix(ncol=ncol(mat_data)+1,nrow=0) # makes empty matrix
for(i in 1:no.clust) {
  cluster=cbind(subset(mat_data, mat_data[,1]==clust.med[i]), "clust.sort"=i)
  data.ord=rbind(data.ord,cluster)
} 
data.ord=data.ord[order(data.ord[,6],data.ord[,2]),] #sorts by new clust.sort column then by sort col.
write.table(data.ord[,-1],paste0("DEGs_filtered_clara_k",no.clust,"_reorder.txt"), sep="\t")
data <- read.delim(paste0("DEGs_filtered_clara_k",no.clust,"_reorder.txt"), header=TRUE, row.names=1)

#List SoLyc IDs for each cluster
for(i in 1:no.clust)
  write.table(rownames(data[which(data[,"clust.sort"]==i),]),paste0("DEGs_filtered_clara_k",no.clust,"_clust",as.character(i),"_GENE"),sep="\t",row.names=F,col.names=F)

##Make HeatMap, maybe offload this to a separate R script
library(gplots)
library(RColorBrewer)
pdf(paste0("DEGs_filtered_clara_k",no.clust,"heatmap.pdf"))
heatmap.2(as.matrix(data[,1:4]),
          main = "Tomato Pericarp RNA-Seq",
          notecol="black",
          col=colorRampPalette(c("blue","white", "red"))(n = 100),
          Colv=FALSE,
          Rowv=FALSE,
          labRow=FALSE,    
          cexRow=1,
          dendrogram="none",
          keysize=1,                                            #key fontsize
          cexCol=1,                                             #changes font size of the column labels
          trace="none",                                         # removes trace lines from heatmap
          density.info="none",                             # "histogram", "density" or "none" for the key
          margins=c(10,10),                                       # adjust space for 1-column labels, 2-row labels 
)
dev.off()