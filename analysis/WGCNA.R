#snippet of code for hte WGCNA clustering, this should be added to the glimma script
#I should also remande the glimma script and combine it with the sleuth script to make it more consistent

#this is adapted from http://pklab.med.harvard.edu/scw2014/WGCNA.html
library(WGCNA)
allowWGCNAThreads()
data <- read.delim("Sleuth_splines.txt",check.names=FALSE) 
data <-data[,c(2:5,19)]
data <- data["F.p.value">0.05,] #filter by F.p.value
data <- data[1:500,]
#WGCNA wants genes as columns
data <- as.data.frame(t(data))

##Pick a power threshold
powers = c(c(1:10), seq(from = 12, to=40, by=2));
sft=pickSoftThreshold(data,dataIsExpr = TRUE, powerVector = powers, corFnc = cor, corOptions = list(use = 'p'), networkType = "unsigned")
# Plot the results
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n", main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");
# Red line corresponds to using an R^2 cut-off
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
#The best power is 20 because the fit curve levels off there.
#@set the power and make the adjacency matri
softPower=20
adjacency = adjacency(data, power = softPower)

# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
##maybe merge similar clusters
# Calculate eigengenes
MEList = moduleEigengenes(data, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

##Based on this I want to merge things with a 0.8 similarity, so I cut at 0.2
MEDissThres = 0.2
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(data, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

##plot the dendro again after the merge
sizeGrWindow(12, 9)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

####I think from here I can just cbind the gene names, cluster/color assignment, and FC values together to get something for heatmapping
####I can probably also use the MEs object to get the "average" for each cluster for plotting
####I'm a little concerned that the FC values are negative sometimes, I may need ot add something to them to make the positive


#Blockwise 
bwnet = blockwiseModules(data, maxBlockSize = 500, power = 20, TOMType = "unsigned", minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.15, numericLabels = TRUE, saveTOMs = TRUE, saveTOMFileBase = "SlycSubsetTOM-blockwise",verbose = 3)
##Plot the results
# Relabel blockwise modules
bwLabels = matchLabels(bwnet$colors, moduleLabels);
# Convert labels to colors for plotting
bwModuleColors = labels2colors(bwLabels)
plotDendroAndColors(bwnet$dendrograms[[1]],bwnet$blockGenes[[1]])




