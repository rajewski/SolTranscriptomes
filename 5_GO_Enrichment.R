#Install libraries
#source("https://bioconductor.org/biocLite.R")
#biocLite("Rgraphviz")
library(topGO)
library(dplyr)
library("tidyr")
library(Rgraphviz)

pal <- c("#74A9CF", "#3690C0", "#0570B0", "#045A8D", "#023858")

# Clean the lists of GO terms from the bash script and convert to a named list object
GOSlyc <- read.table("DEGAnalysis/Pfam/Slyc.gene2go.tsv",
                     stringsAsFactors = F)
GOSlyc <- separate_rows(as.data.frame(GOSlyc[,c(1,2)]), 2, sep="\\|")
GOSlyc <- GOSlyc %>% 
  distinct() %>% 
  group_by(V1) %>% 
  mutate(V2 = paste0(V2, collapse = "|")) %>%
  distinct()
GOSlyc <- setNames(strsplit(x = unlist(GOSlyc$V2), split = "\\|"), GOSlyc$V1)

# Read in the genes of interest o be tested for GO enrichment
SlycGOI <- scan(file="DEGAnalysis/RNA-seq/SlycIH_3Stage_Cluster_1.txt",
                  what = list(""))[[1]]
SlycGOI <- factor(as.integer(names(GOSlyc) %in% names(SlycGOI)))
names(SlycGOI) <- names(GOSlyc)

# Create a TopGO object with the necessary data, also get a summary with ()
(SlycGOData <- new("topGOdata", 
                  description = "Slyc Cluster 1",
                  ontology = "BP",
                  allGenes = SlycGOI,
                  nodeSize = 5,
                  annot = annFUN.gene2GO,
                  gene2GO=GOSlyc))
# Do the enrichment test
(SlycGOResults <- runTest(SlycGOData, algorithm="parentchild", statistic = "fisher"))
# Summarize the test with a table
(SylcGoTable <- GenTable(SlycGOData,
                          Fisher = SlycGOResults,
                          orderBy = "Fisher",
                          ranksOf = "Fisher",
                          topNodes = 30))
# Summarize the test with a hierarchical graph
dev.off()
showSigOfNodes(SlycGOData, score(SlycGOResults), firstSigNodes = 10, useInfo = "all")

# Summarize with a prettier graph from https://www.biostars.org/p/350710/
SylcGoGraph <- SylcGoTable[as.numeric(SylcGoTable$Fisher)<0.05,
                           c("GO.ID", "Term", "Fisher", "Significant")]
SylcGoGraph$Term <- gsub(" [a-z]*\\.\\.\\.$", "", SylcGoGraph$Term) #clean elipses
SylcGoGraph$Term <- gsub("\\.\\.\\.$", "", SylcGoGraph$Term)
SylcGoGraph$Term <- paste(SylcGoGraph$GO.ID, SylcGoGraph$Term, sep=", ")
SylcGoGraph$Term <- factor(SylcGoGraph$Term, levels=rev(SylcGoGraph$Term))
SylcGoGraph$Fisher <- as.numeric(SylcGoGraph$Fisher)
require(ggplot2)
SlycGoPlot <- ggplot(SylcGoGraph, aes(x=Term, y=-log10(Fisher), fill=Significant)) +
  stat_summary(geom = "bar", fun.y = mean, position = "dodge") +
  xlab("Biological Process") +
  ylab("Log Fold Enrichment") +
  scale_fill_gradientn(colours = pal) +
  ggtitle("Tomato Cluster 1") +
  scale_y_continuous(breaks = round(seq(0, max(-log10(SylcGoGraph$Fisher)), by = 2), 1)) +
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
SlycGoPlot

