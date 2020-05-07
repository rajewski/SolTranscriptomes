#source("https://bioconductor.org/biocLite.R")
#biocLite("GENIE3")
setwd("~/bigdata/Slycopersicum/analysis/")
library(GENIE3)
library(doRNG)
library(tidyr)
library(igraph)

#--------Make external data the first time
# # # #Read in data from sleuth or a csv
# # # so <- readRDS("Sleuth_so.rds") #read in sleuth object with all sample expression
# # # data <- so$obs_raw[,c(1,2,6)] %>% spread(sample,est_counts) #reformat for analysis
# # # rownames(data) <- data[,1] #give row names
# # # data <- data[,-1] #drop names
# # # data <- data[,c(4:9,1:3,10:15)] #reorder columns by time
# # # write.csv(data, file="GenieInput.csv")
# # data <- read.csv("GenieInput.csv", row.names = 1, check.names = FALSE)
# #Make external list of links above a very low threshold
# for (i in 1:15) {
#   clust <- read.csv(paste0("DEGs_filtered15_Clust",i,"GENE"), header=F)
#   clustdata <- data[row.names(data) %in% clust$V1,] 
#   network <- GENIE3(as.matrix(clustdata), nCores=24, verbose = TRUE)
#   linklist <- getLinkList(network)
#   linklist <- linklist[linklist$weight>=0.001,]
#   write.csv(linklist, file=paste0("Clust",i,"_LinkList.csv"))
# }

##Cluster 1 that contains FUL1
linklist1 <- read.csv("Clust1_LinkList.csv")[,-c(1)]
linklist1 <- linklist1[linklist1$weight>=linklist1[1,3]*.2,]
net1 <- graph_from_data_frame(d=linklist1, directed=T) 
net1 <- simplify(net1, remove.multiple = T, remove.loops = T)
#Plot FUL2 network
net1_FUL2 <- induced.subgraph(graph=net1,vids=unlist(neighborhood(graph=net1,order=1,nodes="Solyc03g114830.3.1")))
V(net1_FUL2)["Solyc03g114830.3.1"]$color<-"red"
pdf(file="Clust1_FUL2_Network.pdf", width=8, height=8)
  plot(net1_FUL2, edge.arrow.size=0.4, vertex.label=NA, layout=layout_with_fr(net1_FUL2, niter=10000, grid = "nogrid"))
dev.off()
#plot top network
net1_top <- induced.subgraph(graph=net1,vids=unlist(neighborhood(graph=net1,order=1,nodes=as.character(linklist1[1,1]))))
V(net1_top)[as.character(linklist1[1,1])]$color<-"red"
pdf(file="Clust1_Top_Network.pdf", width=8, height=8)
  plot(net1_top, edge.arrow.size=0.4, vertex.label=NA, layout=layout_with_fr(net1_top, niter=10000, grid = "nogrid"))
dev.off()

##Cluster 2 that contains MC and MBP20
linklist2 <- read.csv("Clust2_LinkList.csv")[,-c(1)]
linklist2 <- linklist2[linklist2$weight>=linklist2[1,3]*.2,]
net2 <- graph_from_data_frame(d=linklist2, directed=T) 
net2 <- simplify(net2, remove.multiple = T, remove.loops = T)
#Plot MC network
net2_MC <- induced.subgraph(graph=net2,vids=unlist(neighborhood(graph=net2,order=1,nodes="Solyc05g056620.2.1")))
V(net2_MC)["Solyc05g056620.2.1"]$color<-"red"
pdf(file="Clust2_MC_Network.pdf", width=8, height=8)
  plot(net2_MC, edge.arrow.size=0.4, vertex.label=NA, layout=layout_with_fr(net2_MC, niter=10000))
dev.off()
#Plot MBP20 network
net2_MBP20 <- induced.subgraph(graph=net2,vids=unlist(neighborhood(graph=net2,order=1,nodes="Solyc02g089210.3.1")))
V(net2_MBP20)["Solyc02g089210.3.1"]$color<-"red"
pdf(file="Clust2_MBP20_Network.pdf", width=8, height=8)
  plot(net2_MBP20, edge.arrow.size=0.4, vertex.label=NA, layout=layout_with_fr(net2_MBP20, niter=10000))
dev.off()
#plot top network
net2_top <- induced.subgraph(graph=net2,vids=unlist(neighborhood(graph=net2,order=1,nodes=as.character(linklist2[1,1]))))
V(net2_top)[as.character(linklist2[1,1])]$color<-"red"
pdf(file="Clust2_Top_Network.pdf", width=8, height=8)
  plot(net2_top, edge.arrow.size=0.4, vertex.label=NA, layout=layout_with_fr(net2_top, niter=10000))
dev.off()

##Cluster 3 that contains FUL1
linklist3 <- read.csv("Clust3_LinkList.csv")[,-c(1)]
linklist3 <- linklist3[linklist3$weight>=linklist3[1,3]*.2,]
net3 <- graph_from_data_frame(d=linklist3, directed=T) 
net3 <- simplify(net3, remove.multiple = T, remove.loops = T)
#Plot FUL1 network
net3_FUL1 <- induced.subgraph(graph=net3,vids=unlist(neighborhood(graph=net3,order=1,nodes="Solyc06g069430.3.1")))
V(net3_FUL1)["Solyc06g069430.3.1"]$color<-"red"
pdf(file="Clust3_FUL1_Network.pdf", width=8, height=8)
  plot(net3_FUL1, edge.arrow.size=0.4, vertex.label=NA, layout=layout_with_fr(net3_FUL1, niter=10000))
dev.off()
#plot top network
net3_top <- induced.subgraph(graph=net3,vids=unlist(neighborhood(graph=net3,order=1,nodes=as.character(linklist3[1,1]))))
V(net3_top)[as.character(linklist3[1,1])]$color<-"red"
pdf(file="Clust3_Top_Network.pdf", width=8, height=8)
  plot(net3_top, edge.arrow.size=0.4, vertex.label=NA, layout=layout_with_fr(net3_top, niter=10000))
dev.off()

###
#Skipping clusters 4-7 because they have no interesting genes
###

##Cluster 8 that contains JOINTLESS
linklist8 <- read.csv("Clust8_LinkList.csv")[,-c(1)]
linklist8 <- linklist8[linklist8$weight>=linklist8[1,3]*.2,]
net8 <- graph_from_data_frame(d=linklist8, directed=T) 
net8 <- simplify(net8, remove.multiple = T, remove.loops = T)
#Plot JOINTLESS network
net8_J <- induced.subgraph(graph=net8,vids=unlist(neighborhood(graph=net8,order=1,nodes="Solyc11g010570.2.1")))
V(net8_J)["Solyc11g010570.2.1"]$color<-"red"
pdf(file="Clust8_J_Network.pdf", width=8, height=8)
  plot(net8_J, edge.arrow.size=0.4, vertex.label=NA, layout=layout_with_fr(net8_J, niter=10000))
dev.off()
#plot top network
net8_top <- induced.subgraph(graph=net8,vids=unlist(neighborhood(graph=net8,order=1,nodes=as.character(linklist8[1,1]))))
V(net8_top)[as.character(linklist8[1,1])]$color<-"red"
pdf(file="Clust8_Top_Network.pdf", width=8, height=8)
  plot(net8_top, edge.arrow.size=0.4, vertex.label=NA, layout=layout_with_fr(net8_top, niter=10000))
dev.off()

##Repeat with Cluster 14 that contains CNR
linklist14 <- read.csv("Clust14_LinkList.csv")[,-c(1)]
linklist14 <- linklist14[linklist14$weight>=linklist14[1,3]*.2,]
net14 <- graph_from_data_frame(d=linklist14, directed=T) 
net14 <- simplify(net14, remove.multiple = T, remove.loops = T)
#Make CNR subgraph
net14_CNR <- induced.subgraph(graph=net14,vids=unlist(neighborhood(graph=net14,order=1,nodes="Solyc02g077920.3.1")))
V(net14_CNR)["Solyc02g077920.3.1"]$color<-"red"
pdf(file="Clust14_CNR_Network.pdf", width=8, height=8)
  plot(net14_CNR, vertex.label=NA, edge.arrow.size=.4, layout=layout_with_fr(net14_CNR, niter=10000))
dev.off()
#make top vertex graph
net14_top <- induced.subgraph(graph=net14,vids=unlist(neighborhood(graph=net14,order=1,nodes=as.character(linklist14[1,1]))))
V(net14_top)[as.character(linklist14[1,1])]$color<-"red"
pdf(file="Clust14_Top_Network.pdf", width=8, height=8)
  plot(net14_top, vertex.label=NA, edge.arrow.size=.4, layout=layout_with_fr(net14_top, niter=10000))
dev.off()

#make a network of just cluster 15 from the glimma script
linklist15 <- read.csv("Clust15_LinkList.csv")[,-c(1)]
linklist15 <- linklist15[linklist15$weight>=linklist15[1,3]*.2,]
net15 <- graph_from_data_frame(d=linklist15, directed=T) 
net15 <- simplify(net15, remove.multiple = T, remove.loops = T)
#Plot NOR gene network
net15_NOR <- induced.subgraph(graph=net15,vids=unlist(neighborhood(graph=net15,order=1,nodes="Solyc10g006880.3.1")))
V(net15_NOR)["Solyc10g006880.3.1"]$color<-"red"
pdf(file="Clust15_NOR_Network.pdf", width=8, height=8)
  plot(net15_NOR, vertex.label=NA, edge.arrow.size=.4, layout=layout_with_fr(net15_NOR, niter=10000, grid = "nogrid"))
dev.off()
#Plot RIN gene network
net15_RIN <- induced.subgraph(graph=net15,vids=unlist(neighborhood(graph=net15,order=1,nodes="Solyc05g012020.3.1")))
V(net15_RIN)["Solyc05g012020.3.1"]$color<-"red"
pdf(file="Clust15_RIN_Network.pdf", width=8, height=8)
  plot(net15_RIN, vertex.label=NA, edge.arrow.size=.4, layout=layout_with_fr(net15_RIN, niter=10000))
dev.off()
#Plot top gene network
net15_top <- induced.subgraph(graph=net15,vids=unlist(neighborhood(graph=net15,order=1,nodes=c(as.character(linklist15[1,1])))))
V(net15_top)[as.character(linklist15[1,1])]$color<-"red"
pdf(file="Clust15_Top_Network.pdf", width=8, height=8)
  plot(net15_top, vertex.label=NA, edge.arrow.size=.4, layout=layout_with_fr(net15_top, niter=10000))
dev.off()

#try with all genes but using FUL genes as regulators
fulgenes <- c("Solyc06g069430.3.1", "Solyc03g114830.3.1", "Solyc02g065730.2.1", "Solyc02g089210.3.1")
#networkFUL <- GENIE3(as.matrix(data), regulators=fulgenes, nCores = 24, verbose = TRUE)
#linkListFUL <- getLinkList(networkFUL)
#write.csv(linkListFUL, file="FRUITFULL_targets.csv")
linklistFUL <- read.csv("FRUITFULL_targets.csv")[,-c(1)]
linklistFUL <- linklistFUL[linklistFUL$weight>=linklistFUL[1,3]*.9,]
netFUL <- graph_from_data_frame(d=linklistFUL, directed=T) 
netFUL <- simplify(netFUL, remove.multiple = T, remove.loops = T)
V(netFUL)[fulgenes]$color<-"red"
netFUL_all <- induced.subgraph(graph=netFUL,vids=unlist(neighborhood(graph=netFUL,order=1,nodes=as.character(fulgenes))))
l <- layout_with_dh(netFUL_all)
pdf(file="All_FUL_Network.pdf", width=8, height=8)
  plot(netFUL_all, edge.arrow.size=0.04, vertex.label=NA, layout=l)
dev.off()
