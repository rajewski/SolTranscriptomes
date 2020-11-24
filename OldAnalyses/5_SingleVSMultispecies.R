library("DESeq2")
library("tibble")
library("dplyr")
library("UpSetR")
source("X_Functions.R")

#Find DEGs that are lost by doing a multispecies analysis


# Read in multispecies datasets
DDS_Nobt <- readRDS("DEGAnalysis/RNA-seq/DDS_NobtSE.rds")
DDS_TAIR <- readRDS("DEGAnalysis/RNA-seq/DDS_TAIR.rds")
DDS_Slyc <- readRDS("DEGAnalysis/RNA-seq/DDS_Slyc.rds")
DDS_Spimp <- readRDS("DEGAnalysis/RNA-seq/DDS_Spimp.rds")
DDS_MS_Model1 <- readRDS("DEGAnalysis/RNA-seq/DDS_AllOrtho_Noise.rds")
DDS_MS_Model2 <- readRDS("DEGAnalysis/RNA-seq/DDS_AllOrtho_DEGByFruit.rds")
DDS_MS_Model3 <- readRDS("DEGAnalysis/RNA-seq/DDS_AllOrtho_DEGBySpecies.rds")

# Read in orthgroup mapping
Orthologs <- read.table("Orthofinder/OrthoFinder/Results_May17/Orthogroups/Orthogroups.tsv",
                        header = T, sep="\t")

# Filter by FDR and add in gene names of orthogroups
Results_Model1 <- results(DDS_MS_Model1) %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < 0.01)
Results_Model1 <- merge(Results_Model1, Orthologs, by.x="gene", by.y="Orthogroup")

Results_Model2 <- results(DDS_MS_Model2) %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < 0.01)
Results_Model2 <- merge(Results_Model2, Orthologs, by.x="gene", by.y="Orthogroup")

Results_Model3 <- results(DDS_MS_Model3) %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < 0.01)
Results_Model3 <- merge(Results_Model3, Orthologs, by.x="gene", by.y="Orthogroup")

# Filter the individual species data sets
Results_Nobt <- results(DDS_Nobt) %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < 0.01)

Results_TAIR <- results(DDS_TAIR) %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < 0.01)

Results_Slyc <- results(DDS_Slyc) %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < 0.01)

Results_Spimp <- results(DDS_Spimp) %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < 0.01)

DEGOverlap <- function(ResultsObject,
                       longName=""){
  List_Single <- list(as.data.frame(unique(ResultsObject$gene)))
  List_Model1 <- list(as.data.frame(unique(Results_Model1[longName])))
  List_Model2 <- list(as.data.frame(unique(Results_Model2[longName])))
  List_Model3 <- list(as.data.frame(unique(Results_Model3[longName])))
  Overlap <- Map(list,List_Single, List_Model1, List_Model2, List_Model3)
  GeneList <- as.data.frame(unique(unlist(Overlap)))
  for (GeneTable in c(List_Single, List_Model1, List_Model2, List_Model3)) {
    Common<-as.data.frame(GeneList[,1] %in% as.data.frame(GeneTable)[,1])
    Common[,1][Common[,1]==TRUE]<- 1
    GeneList <-cbind(GeneList, Common)
  }
  colnames(GeneList) <- c("ID","Single Species","Model 1", "Model 2", "Model 3")
  upsetplot <- upset(GeneList[c(5:2)], 
        sets = colnames(GeneList)[c(5:2)],
        mainbar.y.label = "Number of DEGs", 
        sets.x.label = NULL, 
        keep.order = T,
        text.scale = 1.5,
        point.size = 3, 
        line.size = 1,
        nintersects = 100)
  return(upsetplot)
}

Nobtupset <- DEGOverlap(Results_Nobt, longName = "Nicotiana")
TAIRupset <- DEGOverlap(Results_TAIR, longName = "Arabidopsis")
Slycupset <- DEGOverlap(Results_Slyc, longName = "Solanum")
Spimpupset <- DEGOverlap(Results_Spimp, longName = "Solanum")


