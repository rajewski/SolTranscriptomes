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

# Make a file of orthogroups to GO terms, via arabidopsis
Orthogroups <- read.table("Orthofinder/OrthoFinder/Results_Aug31/Orthogroups/Orthogroups.tsv",
                          sep="\t",
                          stringsAsFactors = F,
                          header=T)
Orthogroups <- Orthogroups %>% filter_all(all_vars(!grepl(',',.))) #Remove multiples
Orthogroups <- Orthogroups %>% filter_all(all_vars(!grepl("^$",.))) # Remove empties
Orthogroups$Cucumis <- strtrim(Orthogroups$Cucumis,12)
GOSlyc <- read.table("DEGAnalysis/Pfam/Slyc.gene2go.tsv", stringsAsFactors = F)
GOSlyc <- merge(GOSlyc, Orthogroups, by.x=1, by.y="Solanum")[,c("Orthogroup", "V2")]
GOTAIR <- read.table("DEGAnalysis/Pfam/TAIR10.gene2go.tsv", stringsAsFactors = F)
GOTAIR <- merge(GOTAIR, Orthogroups, by.x=1, by.y="Arabidopsis")[,c("Orthogroup", "V2")]
GONobt <- read.table("DEGAnalysis/Pfam/Nobt.gene2go.tsv", stringsAsFactors = F)
GONobt <- merge(GONobt, Orthogroups, by.x=1, by.y="Nicotiana")[,c("Orthogroup", "V2")]
GOMelon <- read.table("DEGAnalysis/Pfam/Melon.gene2go.tsv", stringsAsFactors = F, sep="\t")[,1:2]
GOMelon <- merge(GOMelon, Orthogroups, by.x=1, by.y="Cucumis")[,c("Orthogroup", "V2")]
GO <- rbind(GOSlyc, GOTAIR, GONobt, GOMelon) #make superlist of of all GO terms across species
GO <- separate_rows(as.data.frame(GO[,c(1,2)]), 2, sep="\\|")
GO <- GO %>% 
  distinct() %>% 
  group_by(Orthogroup) %>% 
  mutate(V2 = paste0(V2, collapse = "|")) %>%
  distinct()
write.table(GO, "DEGAnalysis/Pfam/Ortho.831All.gene2go.tsv", sep="\t", quote=F, col.names = F, row.names = F)
