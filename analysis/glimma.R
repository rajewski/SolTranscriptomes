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
saveRDS(pam,file = "pam.RDS")
pam<-readRDS("pam.RDS")
save(pam,file=paste0("DEGs_filtered_pam_k",no.clust,"_euc"))
write.table(pam$clustering, file=paste0("DEGs_filtered_pam",no.clust,"_euc_clust.txt"), sep="\t")
pam1=read.delim("DEGs_filtered_pam15_euc_clust.txt")
colnames(pam1)="Cluster"
pam2<-cbind(pam1,data[,2:5])
write.table(pam2, "DEGs_filtered_pam15.txt", quote=FALSE, sep="\t", col.names = NA)

#Set outputs (wtf does this mean?)
FileTag=paste0("DEGs_filtered", no.clust)
HeatmapFileFormat=1  # 1=pdf or 2=tif 
HeatmapFile=paste(FileTag,"_heatmap",sep="") 
GO_IN=paste(FileTag,"_GO_IN.txt",sep="")

#Set parameters
sort.col="X"                    # column that will be used to order genes within clusters
col.clust.sort="FC.1.3"      # column on which to order the clusters based on their median value
col.use=""                          # DataFile columns to use (all="")

#prep data
mat_data <- data.matrix(read.delim("DEGs_filtered_pam15.txt", row.names=1))
colnames(mat_data) <- c("Cluster", "FC.1.3", "FC.3.15", "FC.15.br", "FC.br.rr")
#reorder the rows by cluster number. You'll need to specify the # of FC columns in the apply function
clust.med=cbind(1:no.clust, t(sapply(X=1:no.clust,FUN=function(i) {apply(mat_data[which(mat_data[,"Cluster"]==i),2:5], MARGIN=2,FUN=median)})))
clust.med=clust.med[order(clust.med[,col.clust.sort]),]   # sorts the median table and therefore clusters by the selected sample
clust.ord=rev(clust.med[,1])                              # makes vector with descending median of clusters..
clust.med[1:3,]                                           # checkpoint
clust.ord                                                 # checkpoint
write.table(clust.med,paste(FileTag,"_clust-med.txt",sep=""),sep="\t")
data.ord=matrix(ncol=ncol(mat_data)+1,nrow=0) # makes empty matrix
#resorts the data and adds a column with the resorted cluster order
for(i in 1:length(clust.ord)) {
  cluster=cbind(subset(mat_data, mat_data[,1]==clust.ord[i]), "clust.sort"=i)
  data.ord=rbind(data.ord,cluster)
  } 
data.ord=data.ord[order(data.ord[,5],data.ord[,"clust.sort"]),] #sorts by new clust.sort column then by sort col.
write.table(clust.ord,paste(FileTag,"_reorder.txt",sep=""), sep="\t")
data.ord[1:3,]                                            # checkpoint 
data=as.matrix(data.ord[,-c(1)]) #drop cluster numbers?
data[1:3,]
write.table(data, "DEGs_filtered_pam15reordered.txt", quote=FALSE, sep="\t", col.names = NA)
data<-read.delim("DEGs_filtered_pam15reordered.txt",header=TRUE, row.names=1) #read in cluster no. file

#For A GO Analysis, which I guess I'll try
#preps INPUT file for GO Analysis (3 columns: AGI, cluster number, number of members of the cluster) from cluster file
clustGO=read.delim("DEGs_filtered_pam15.txt")[1:2] #read in cluster no. file
clustGO<-clustGO[1:2]
colnames(clustGO)=c("GENE","Cluster") ### renaming columns
clustGO[1:3,] # checkpoint
clustGO_IN=clustGO
for(i in 1:length(clust.ord)) clustGO_IN[,"Cluster"]=replace(clustGO_IN[,"Cluster"],which(clustGO[,"Cluster"]==clust.ord[i]),i) #replaces cluster numbers with reordered cluster numbers
clustGO_IN[1:5,]
clustGO_IN=clustGO_IN[order(clustGO_IN[,"Cluster"]),] # sorts the file by cluster
clustGO_IN[1:5,]
clustGO_IN_sum=table(clustGO_IN[,2])                                   # makes vector with no. genes in each cluster
clustGO_IN[,3]=clustGO_IN_sum[clustGO_IN[,2]]                          # adds column to data with the no. genes in cluster as third column.
clustGO_IN[1:5,]
colnames(clustGO_IN)=c("GENE","Cluster","no.members")                   # renames the column names 
clustGO_IN[1:5,]                                                        # shows chunk of the table
write.table(clustGO_IN,file=GO_IN,row.name=F,sep="\t",quote=F)         # writes the cluster reorder sequence to a file

# makes vector with row numbers to specify cluster breaks position with horizontal white lines on the heatmap
clust.div=clustGO_IN_sum[1]
for(i in 2:length(clustGO_IN_sum)) clust.div=c(clust.div,clust.div[i-1]+clustGO_IN_sum[i])
names(clust.div)=paste("c",1:length(clustGO_IN_sum),sep="")
clust.div
#saveRDS(clust.div,file="clust.div.RDS") # save rds for making heatmap separately
no.members<-clustGO_IN$no.member

#makes a file with the SoLyc ID list for each cluster
for(i in 1:max(clustGO_IN[,"Cluster"]))
  write.table(clustGO_IN[which(clustGO_IN[,"Cluster"]==i),1],paste(FileTag,"_Clust",as.character(i),"GENE",sep=""),sep="\t",row.names=F,col.names=F)

##Make HeatMap
library(gplots)
library(RColorBrewer)
my_palette <- colorRampPalette(c("blue","white", "red"))(n = 299)

mats=data[order(data[,"clust.sort"]),] 
write.table(mats, "DEGsmatsordered.txt", quote=FALSE, sep="\t", col.names = NA)
data<-read.delim("DEGsmatsordered.txt",header=TRUE, row.names=1) #read in cluster no. file
exp<-data[,1:4] ## your fold change values are in these columns, change numbers if you have more or less columns!
cluster<-data[,5] ## your cluster numbers column, change column number if you have more or less columns!
data<-cbind(cluster,exp)
data=data[order(data[,"cluster"]),] ##
mats<-as.matrix(data[,2:5]) ###
saveRDS(mats,file="mats.RDS")

#Ive offloaded this to a separate script to run as a batch job
pdf("DEGs_filtered_pam15heatmap.pdf")
heatmap.2(mats,
          main = "Tomato Pericarp RNA-Seq",
          notecol="black",
          col=my_palette,            # this sets the color pallet made above
          Colv=FALSE,
          Rowv=FALSE,
          labRow=F,    
          cexRow=1,
          #labCol=fignames[col.use],                  # override column names for figure
          dendrogram="none",
          keysize=1,                                            #key fontsize
          cexCol=1,                                             #changes font size of the column labels
          trace="none",                                         # removes trace lines from heatmap
          density.info="none",                             # "histogram", "density" or "none" for the key
          #colsep=c(3,6,11,15,18,21,24,29),                     # column indices of vertical white lines 
          rowsep=clust.div,                                     # row indices of horizontal white lines
          sepcolor="black",
          sepwidth=c(0.001,0.001),                                  # sets thickness of white lines
          margins=c(10,10),                                       # adjust space for 1-column labels, 2-row labels 
          #lwid=c(0.1,1),                                       # 2x2 relat. widths - Left-key & Rdendo, Right-Cdendo & heatmap
          #lhei=c(0.1,1),                                       # 2x2 relat. heights - Top-key & Cdendo, Bott.-Rdendo & heatmap
          breaks=seq(-3,3,by=6/length(my_palette))             # heat range and bin size - must match no. of colors
)
dev.off()
