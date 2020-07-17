# Install libraries
library(topGO)
library(dplyr)
library(tidyr)
library(Rgraphviz)
source("X_Functions.R")

for (i in as.numeric(gsub("\\D",
                          "",
                          list.files(path="DEGAnalysis/RNA-seq/Lists/",
                                     pattern="^DryOrtho_Cluster_*")))) {
Table <- GOEnrich(gene2go="DEGAnalysis/Pfam/Slyc.gene2go.tsv",
                      GOIs=paste0("DEGAnalysis/RNA-seq/Lists/Solanum_3DF_Noise_Cluster_",i,".txt"))
Plot <- GOPlot(Table, Title=paste0("Solanum Conserved Cluster ", i))
if(is.null(Plot)) {
  next
}
ggsave(paste0("DEGAnalysis/Pfam/Plots/Solanum_3DF_Noise_GO_AllGenes.pdf"), height=8, width=14)
}

Table <- GOEnrich(gene2go = "DEGAnalysis/Pfam/Slyc.gene2go.tsv",
                  GOIs="DEGAnalysis/RNA-seq/Lists/Solanum_3DF_Noise_AllGenes.txt")
Plot <- GOPlot(Table, Title="Conserved Solanum Genes")




# Summarize the test with a hierarchical graph
dev.off()
showSigOfNodes(SlycGOData, score(SlycGOResults), firstSigNodes = 10, useInfo = "all")

# Make a file of orthogroups to GO terms, via tomato
Orthogroups <- read.table("Orthofinder/OrthoFinder/Results_Apr23/Orthogroups/Orthogroups.tsv",
                          sep="\t",
                          stringsAsFactors = F,
                          header=T)
Orthogroups <- Orthogroups %>% filter_all(all_vars(!grepl(',',.))) #Remove multiples
Orthogroups <- Orthogroups[,c(1,4)]
Orthogroups <- Orthogroups %>% filter_all(all_vars(!grepl("^$",.))) # Remove empties
GO <- read.table("DEGAnalysis/Pfam/Slyc.gene2go.tsv", stringsAsFactors = F)
GO <- separate_rows(as.data.frame(GO[,c(1,2)]), 2, sep="\\|")
GO <- GO %>% 
  distinct() %>% 
  group_by(V1) %>% 
  mutate(V2 = paste0(V2, collapse = "|")) %>%
  distinct()
Ortho2Go <- merge(GO, Orthogroups, by.x="V1", by.y="Solanum")[,c(3,2)]
write.table(Ortho2Go, "DEGAnalysis/Pfam/Ortho.gene2go.tsv", sep="\t", quote=F, col.names = F, row.names = F)
