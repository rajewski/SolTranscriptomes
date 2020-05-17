# Install libraries
library(topGO)
library(dplyr)
library(tidyr)
library(Rgraphviz)

# Function to do the enrichment
GOEnrich <- function(gene2go="",
                     GOIs="",
                     GOCategory="BP",
                     NumCategories=20) {
  # Clean the lists of GO terms from the bash script and convert to a named list object
  GO <- read.table(gene2go, stringsAsFactors = F)
  GO <- separate_rows(as.data.frame(GO[,c(1,2)]), 2, sep="\\|")
  GO <- GO %>% 
    distinct() %>% 
    group_by(V1) %>% 
    mutate(V2 = paste0(V2, collapse = "|")) %>%
    distinct()
  GO <- setNames(strsplit(x = unlist(GO$V2), split = "\\|"), GO$V1)
  # Read in the genes of interest to be tested for GO enrichment
  GOI <- scan(file=GOIs,
              what = list(""))[[1]]
  GOI <- factor(as.integer(names(GO) %in% GOI))
  names(GOI) <- names(GO)
  # Create a TopGO object with the necessary data, also get a summary with
  GOData <- new("topGOdata",
                description = "Slyc Cluster 1",
                ontology = GOCategory,
                allGenes = GOI,
                nodeSize = 5,
                annot = annFUN.gene2GO,
                gene2GO=GO)
  # Do the enrichment test
  GOResults <- runTest(GOData,
                       algorithm="weight01",
                       statistic = "fisher")
  # Summarize the test with a table
  GOTable <- GenTable(GOData,
                      Fisher = GOResults,
                      orderBy = "Fisher",
                      ranksOf = "Fisher",
                      topNodes = NumCategories)
  return(GOTable)
}

# Function to make a cute plot
GOPlot <- function(GenTable=X,
                   Title="") 
  {
  # Summarize with a prettier graph from https://www.biostars.org/p/350710/
  GenTable$Fisher <- gsub("^< ", "", GenTable$Fisher) #remove too low pvals
  GoGraph <- GenTable[as.numeric(GenTable$Fisher)<0.05, c("GO.ID", "Term", "Fisher", "Significant")]
  if(dim(GoGraph)[1]==0){
    return(NULL) #stop empty graphs
  }
  GoGraph$Term <- gsub(" [a-z]*\\.\\.\\.$", "", GoGraph$Term) #clean elipses
  GoGraph$Term <- gsub("\\.\\.\\.$", "", GoGraph$Term)
  GoGraph$Term <- paste(GoGraph$GO.ID, GoGraph$Term, sep=", ")
  GoGraph$Term <- factor(GoGraph$Term, levels=rev(GoGraph$Term))
  GoGraph$Fisher <- as.numeric(GoGraph$Fisher)
  require(ggplot2)
  pal <- c("#74A9CF", "#3690C0", "#0570B0", "#045A8D", "#023858")
  GoPlot <- ggplot(GoGraph, aes(x=Term, y=-log10(Fisher), fill=Significant)) +
    stat_summary(geom = "bar", fun = mean, position = "dodge") +
    xlab("Biological Process") +
    ylab("Log Fold Enrichment") +
    scale_fill_gradientn(colours = pal) +
    ggtitle(Title) +
    scale_y_continuous(breaks = round(seq(0, max(-log10(GoGraph$Fisher),2), by = 2), 1)) +
    theme_bw(base_size=24) +
    theme(
      panel.grid = element_blank(),
      legend.position=c(0.9,0.2),
      legend.background=element_blank(),
      legend.key=element_blank(),     #removes the border
      legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
      legend.text=element_text(size=18),  #Text size
      legend.title=element_blank(),
      plot.title=element_text(angle=0, size=24, face="bold", vjust=1),
      axis.text.x=element_text(angle=0, size=18, face="bold", hjust=1.10),
      axis.text.y=element_text(angle=0, size=18, face="bold", vjust=0.5),
      axis.title=element_text(size=24, face="bold"),
      title=element_text(size=18)) +
    guides(fill=guide_colorbar(ticks=FALSE, label.position = 'left')) +
    coord_flip()
  print(GoPlot)
}

for (i in as.numeric(gsub("\\D",
                          "",
                          list.files(path="DEGAnalysis/RNA-seq/Lists/",
                                     pattern="^Solanum_Cluster_*")))) {
Table <- GOEnrich(gene2go="DEGAnalysis/Pfam/Slyc.gene2go.tsv",
                      GOIs=paste0("DEGAnalysis/RNA-seq/Lists/Solanum_Cluster_",i,".txt"))
Plot <- GOPlot(Table, Title=paste0("Solanum spp.  Cluster ", i))
if(is.null(Plot)) {
  next
}
ggsave(paste0("DEGAnalysis/Pfam/Plots/Solanum_GO_Cluster_", i, ".pdf"), height=8, width=14)
}

Table <- GOEnrich(gene2go = "DEGAnalysis/Pfam/Ortho.gene2go.tsv",
                  GOIs="DEGAnalysis/RNA-seq/Lists/AllOrtho_DEGBySpecies_Cluster_1.txt")
Plot <- GOPlot(OrthoTable, Title="DEGs by Fruit Type Cluster 1")




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
