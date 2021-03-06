# You need to run the 5_PFAM_Enrichment.sh first because it generates input files for this script
source("X_Functions.R")




# IPR on RNA seq clusters -------------------------------------------------

# Work through the Slyc IH data
for ( i in as.numeric(gsub("\\D",
                           "",
                           list.files(path="DEGAnalysis/RNA-seq/Lists/",
                                      pattern="^Slyc_Cluster_*")))) {
  PfamEnrichment(AllGenesFile = "DEGAnalysis/Pfam/Slyc.protein.names.txt",
                 AllIPRFile = "DEGAnalysis/Pfam/Slyc.gene2ipr.tsv",
                 IPRDescFile = "DEGAnalysis/Pfam/Slyc.ipr2desc.tsv",
                 ExcludedGenesFile = "DEGAnalysis/Pfam/Slyc.nopfam.tsv",
                 TopGenesFile = paste0("DEGAnalysis/RNA-seq/Lists/Slyc_Cluster_", i, ".txt"),
                 OutputFile = paste0("DEGAnalysis/Pfam/Lists/Slyc_IPR_Cluster_", i, ".txt"))
}

# Work through the Spimp IH data
for ( i in as.numeric(gsub("\\D",
                           "",
                           list.files(path="DEGAnalysis/RNA-seq/Lists/",
                                      pattern="^Spimp_Cluster_*")))) {
  PfamEnrichment(AllGenesFile = "DEGAnalysis/Pfam/Slyc.protein.names.txt",
                 AllIPRFile = "DEGAnalysis/Pfam/Slyc.gene2ipr.tsv",
                 IPRDescFile = "DEGAnalysis/Pfam/Slyc.ipr2desc.tsv",
                 ExcludedGenesFile = "DEGAnalysis/Pfam/Slyc.nopfam.tsv",
                 TopGenesFile = paste0("DEGAnalysis/RNA-seq/Lists/Spimp_Cluster_", i, ".txt"),
                 OutputFile = paste0("DEGAnalysis/Pfam/Lists/Spimp_IPR_Cluster_", i, ".txt"))
}

# Work through the Spimp vs Slyc interaction data
for ( i in as.numeric(gsub("\\D",
                           "",
                           list.files(path="DEGAnalysis/RNA-seq/Lists/",
                                      pattern="^Solanum_Cluster_*")))) {
for (i in c(1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18,19,20,21)){
    PfamEnrichment(AllGenesFile = "DEGAnalysis/Pfam/Slyc.protein.names.txt",
                 AllIPRFile = "DEGAnalysis/Pfam/Slyc.gene2ipr.tsv",
                 IPRDescFile = "DEGAnalysis/Pfam/Slyc.ipr2desc.tsv",
                 ExcludedGenesFile = "DEGAnalysis/Pfam/Slyc.nopfam.tsv",
                 TopGenesFile = paste0("DEGAnalysis/RNA-seq/Lists/Solanum_3DF_Noise_Cluster_", i, ".txt"),
                 OutputFile = paste0("DEGAnalysis/Pfam/Lists/Solanum_3DF_Noise_IPR_Cluster_", i, ".txt"))
}


# Work through the TAIR data
for ( i in as.numeric(gsub("\\D",
                           "",
                           list.files(path="DEGAnalysis/RNA-seq/Lists/",
                                      pattern="^TAIR_Cluster_*")))) {
  PfamEnrichment(AllGenesFile = "DEGAnalysis/Pfam/TAIR10.protein.names.txt",
                 AllIPRFile = "DEGAnalysis/Pfam/TAIR10.gene2ipr.tsv",
                 IPRDescFile = "DEGAnalysis/Pfam/TAIR10.ipr2desc.tsv",
                 ExcludedGenesFile = "DEGAnalysis/Pfam/TAIR10.nopfam.tsv",
                 TopGenesFile = paste0("DEGAnalysis/RNA-seq/Lists/TAIR_Cluster_", i, ".txt"),
                 OutputFile = paste0("DEGAnalysis/Pfam/Lists/TAIR_IPR_Cluster_", i, ".txt"))
}

# And for Nobt
for ( i in as.numeric(gsub("\\D",
                      "",
                      list.files(path="DEGAnalysis/RNA-seq/Lists/",
                                 pattern="^Nobt_Cluster_*")))) {
  PfamEnrichment(AllGenesFile = "DEGAnalysis/Pfam/Nobt.protein.names.txt",
                 AllIPRFile = "DEGAnalysis/Pfam/Nobt.gene2ipr.tsv",
                 IPRDescFile = "DEGAnalysis/Pfam/Nobt.ipr2desc.tsv",
                 ExcludedGenesFile = "DEGAnalysis/Pfam/Nobt.nopfam.tsv",
                 TopGenesFile = paste0("DEGAnalysis/RNA-seq/Lists/Nobt_Cluster_", i, ".txt"),
                 OutputFile = paste0("DEGAnalysis/Pfam/Lists/Nobt_IPR_Cluster_", i, ".txt"))
}

# How to deal with the orthogroup data?
library("tidyr")
library("tibble")
library("dplyr")
Orthogroups <- read.table("Orthofinder/OrthoFinder/Results_May17/Orthogroups/Orthogroups.tsv",
                          sep="\t",
                          stringsAsFactors = F,
                          header=T)
Orthogroups <- Orthogroups %>% filter_all(all_vars(!grepl(',',.))) #Remove multiples
Orthogroups <- Orthogroups[,c(1,4)] 
Orthogroups <- Orthogroups %>% filter_all(all_vars(!grepl("^$",.))) # Remove empties
write.table(Orthogroups[,1], "DEGAnalysis/Pfam/Ortho.protein.names.txt", sep="\t", quote=F, row.names = F, col.names = F) #List the gene universe in orthogene format

# Just map the orthogroups to Slyc IPR hits because it is the most complete
SlycIPR <- read.table(file="DEGAnalysis/Pfam/Slyc.gene2ipr.tsv")
OrthoIPR <- merge(SlycIPR, Orthogroups, by.x="V1", by.y="Solanum")[,c(3,2)]
write.table(OrthoIPR, "DEGAnalysis/Pfam/Ortho.gene2ipr.tsv", sep="\t", quote=F, row.names = F, col.names = F)
SlycNoIPR <- read.table(file="DEGAnalysis/Pfam/Slyc.nopfam.tsv")
OrthoNoIPR <- merge(SlycNoIPR, Orthogroups, by.x="V1", by.y="Solanum")
write.table(OrthoNoIPR$Orthogroup, "DEGAnalysis/Pfam/Ortho.nopfam.tsv", quote=F, row.names = F, col.names = F)

# Map tdry fruited orthos to TAIR
Orthogroups <- read.table("Orthofinder/OrthoFinder/Results_May17/Orthogroups/Orthogroups.tsv",
                          sep="\t",
                          stringsAsFactors = F,
                          header=T)
Orthogroups <- Orthogroups[,1:3] 
Orthogroups <- Orthogroups %>% filter_all(all_vars(!grepl(',',.))) #Remove multiples
Orthogroups <- Orthogroups %>% filter_all(all_vars(!grepl("^$",.))) # Remove empties
write.table(Orthogroups[,1], "DEGAnalysis/Pfam/DryOrtho.protein.names.txt", sep="\t", quote=F,row.names = F, col.names = F) 

TAIRIPR <- read.table(file="DEGAnalysis/Pfam/TAIR10.gene2ipr.tsv")
OrthoIPR <- merge(TAIRIPR, Orthogroups, by.x="V1", by.y="Arabidopsis")[,c(3,2)]
write.table(OrthoIPR, "DEGAnalysis/Pfam/DryOrtho.gene2ipr.tsv", sep="\t", quote=F, row.names = F, col.names = F)
TAIRNoIPR <- read.table(file="DEGAnalysis/Pfam/TAIR10.nopfam.tsv")
OrthoNoIPR <- merge(TAIRNoIPR, Orthogroups, by.x="V1", by.y="Arabidopsis")
write.table(OrthoNoIPR$Orthogroup, "DEGAnalysis/Pfam/DryOrtho.nopfam.tsv", quote=F, row.names = F, col.names = F)

# And for Orthos
for ( i in as.numeric(gsub("\\D",
                           "",
                           list.files(path="DEGAnalysis/RNA-seq/Lists/",
                                      pattern="^AllOrtho_DEGBySpecies_Cluster_*")))) {
  PfamEnrichment(AllGenesFile = "DEGAnalysis/Pfam/Ortho.protein.names.txt",
                 AllIPRFile = "DEGAnalysis/Pfam/Ortho.gene2ipr.tsv",
                 IPRDescFile = "DEGAnalysis/Pfam/Slyc.ipr2desc.tsv",
                 ExcludedGenesFile = "DEGAnalysis/Pfam/Ortho.nopfam.tsv",
                 TopGenesFile = paste0("DEGAnalysis/RNA-seq/Lists/AllOrtho_DEGBySpecies_Cluster_", i, ".txt"),
                 OutputFile = paste0("DEGAnalysis/Pfam/Lists/AllOrtho_DEGBySpecies_IPR_Cluster_", i, ".txt"))
}

#For Dry orthos
for ( i in as.numeric(gsub("\\D",
                           "",
                           list.files(path="DEGAnalysis/RNA-seq/Lists/",
                                      pattern="^DryOrtho_Cluster_*")))) {
  PfamEnrichment(AllGenesFile = "DEGAnalysis/Pfam/DryOrtho.protein.names.txt",
                 AllIPRFile = "DEGAnalysis/Pfam/DryOrtho.gene2ipr.tsv",
                 IPRDescFile = "DEGAnalysis/Pfam/TAIR10.ipr2desc.tsv",
                 ExcludedGenesFile = "DEGAnalysis/Pfam/DryOrtho.nopfam.tsv",
                 TopGenesFile = paste0("DEGAnalysis/RNA-seq/Lists/DryOrtho_Cluster_", i, ".txt"),
                 OutputFile = paste0("DEGAnalysis/Pfam/Lists/DryOrtho_IPR_Cluster_", i, ".txt"))
}




# ChIP hits on RNA seq clusters ---------------------------------------------------
# this is to test if a given cluster is enriched for genes with promoters bound by FUL1/2
FULBinding <- read.csv("ChIPAnalysis/ChIP-chip/TomatoFULBinding.tsv",
                       stringsAsFactors = F)
# Table of transcripts and antibody hits
write.table(FULBinding[,c("Transcript", "Antibody")],
          file="ChIPAnalysis/ChIP-chip/TomatoFULTargets.csv",
          row.names = F,
          col.names= F,
          quote = F,
          sep="\t")

# Descriptions of the possible hits
write.table(rbind(c("FUL1.sm", "FUL1_bound"),
                  c("FUL2.sm", "FUL2_bound"),
),
            file="ChIPAnalysis/ChIP-chip/TomatoFULdescription.csv",
            row.names = F,
            col.names= F,
            quote = F,
            sep="\t")
# Table of the zero genes to exclude
write.table(NULL, file="ChIPAnalysis/ChIP-chip/DummyChip.txt")


# Slyc IH data on general antibody binding
for (i in as.numeric(gsub("\\D",
                          "",
                          list.files(path="DEGAnalysis/RNA-seq/",
                                     pattern="^SlycIH_Cluster_*")))){
  PfamEnrichment(AllGenesFile = "DEGAnalysis/Pfam/Slyc.protein.names.txt",
                 AllIPRFile = "ChIPAnalysis/ChIP-chip/TomatoFULTargets.csv",
                 IPRDescFile = "ChIPAnalysis/ChIP-chip/TomatoFULdescription.csv",
                 ExcludedGenesFile = "ChIPAnalysis/ChIP-chip/DummyChip.txt",
                 TopGenesFile = paste0("DEGAnalysis/RNA-seq/SlycIH_Cluster_", i, ".txt"),
                 OutputFile = paste0("ChIPAnalysis/ChIP-chip/SlycIH_Cluster", i, ".txt"))
}

# Slyc SRA data on general antibody binding
for (i in as.numeric(gsub("\\D",
                          "",
                          list.files(path="DEGAnalysis/RNA-seq/",
                                     pattern="^SlycSRA_Cluster_*")))) {
  PfamEnrichment(AllGenesFile = "DEGAnalysis/Pfam/Slyc.protein.names.txt",
                 AllIPRFile = "ChIPAnalysis/ChIP-chip/TomatoFULTargets.csv",
                 IPRDescFile = "ChIPAnalysis/ChIP-chip/TomatoFULdescription.csv",
                 ExcludedGenesFile = "ChIPAnalysis/ChIP-chip/DummyChip.txt",
                 TopGenesFile = paste0("DEGAnalysis/RNA-seq/SlycSRA_Cluster_", i, ".txt"),
                 OutputFile = paste0("ChIPAnalysis/ChIP-chip/SlycSRA_Cluster", i, ".txt"))
}


# Plot distance histograms if the cluster is enriched
PlotChIPHist <- function(DistFile="ChIPAnalysis/ChIP-chip/TomatoFULBinding.tsv",
                         ClusterFile="") {
  Dist2TSS <- read.csv(DistFile)
  Cluster <- read.table(ClusterFile,
                          sep="\t",
                          header=T,
                          stringsAsFactors = F)
  if (length(grep("FUL",Cluster$Description))>1){
    par(mfrow=c(1,2))
  }
  if (length(grep("FUL1",Cluster$Description))>0){
    D2T_1 <- Dist2TSS[(Dist2TSS$Transcript %in% strsplit(Cluster[grep("FUL1",Cluster$Description),"Gene_IDs"], split=",")[[1]] & Dist2TSS$Antibody=="FUL1.sm"),]
    hist(D2T_1$Distance, 
         xlab="Distance from TSS (bp)",
         main=expression(paste(alpha,"-FUL1 Binding")))
  }
  if (length(grep("FUL2",Cluster$Description))>0){
    D2T_2 <- Dist2TSS[(Dist2TSS$Transcript %in% strsplit(Cluster[grep("FUL2",Cluster$Description),"Gene_IDs"], split=",")[[1]] & Dist2TSS$Antibody=="FUL2.sm"),]
    hist(D2T_2$Distance, 
         xlab="Distance from TSS (bp)",
         main=expression(paste(alpha,"-FUL2 Binding")))
  }
}

#Make Slyc IH Histos
for (i in as.numeric(gsub("\\D",
                          "",
                          list.files(path="ChIPAnalysis/ChIP-chip/",
                                     pattern="^SlycIH_Cluster*")))) {
  if(file.size(paste0("ChIPAnalysis/ChIP-chip/SlycIH_Cluster",i, ".txt"))<10){
    next
  }
  pdf(file=paste0("ChIPAnalysis/ChIP-chip/Plot_SlycIH_ChIPHisto_Cluster", i,".pdf"),
      height=6,
      width=11)
  PlotChIPHist(ClusterFile = paste0("ChIPAnalysis/ChIP-chip/SlycIH_Cluster",i, ".txt"))
  dev.off()
}

#Make Slyc SRA Histos
for (i in as.numeric(gsub("\\D",
                          "",
                          list.files(path="ChIPAnalysis/ChIP-chip/",
                                     pattern="^SlycSRA_Cluster*")))) {
  if(file.size(paste0("ChIPAnalysis/ChIP-chip/SlycSRA_Cluster",i, ".txt"))<10){
    next
  }
  pdf(file=paste0("ChIPAnalysis/ChIP-chip/Plot_SlycSRA_ChIPHisto_Cluster", i,".pdf"),
      height=6,
      width=11)
  PlotChIPHist(ClusterFile = paste0("ChIPAnalysis/ChIP-chip/SlycSRA_Cluster",i, ".txt"))
  dev.off()
}





# IPR on ChIP hits --------------------------------------------------------

# Try to find the enrichment of promoters bound by FUL

chiphits <- read.csv("ChIPAnalysis/ChIP-chip/TomatoFULBinding.tsv")
chiphits <- chiphits[chiphits$Distance<100,]
chiphits$Transcript <- as.character(chiphits$Transcript)
chiphits <- split(chiphits$Transcript, chiphits$Antibody)
write.table(chiphits[[1]], 
            file="ChIPAnalysis/ChIP-chip/FUL1hits.tsv",
            col.names=F,
            row.names=F,
            quote=F)
write.table(chiphits[[2]], 
            file="ChIPAnalysis/ChIP-chip/FUL2hits.tsv",
            col.names=F,
            row.names=F,
            quote=F)

PfamEnrichment(AllGenesFile = "DEGAnalysis/Pfam/Slyc.protein.names.txt",
                 AllIPRFile = "DEGAnalysis/Pfam/Slyc.gene2ipr.tsv",
                 IPRDescFile = "DEGAnalysis/Pfam/Slyc.ipr2desc.tsv",
                 ExcludedGenesFile = "DEGAnalysis/Pfam/Slyc.nopfam.tsv",
                 TopGenesFile = "ChIPAnalysis/ChIP-chip/FUL1hits.tsv",
                 OutputFile = paste0("ChIPAnalysis/ChIP-chip/FUL1-bound_IPREnrichment.txt"))

# Plot Enrichments --------------------------------------------------------
# I want to use the typical schema for plotting GO enrichments as for plotting these various enrichments. I'll start with IPR domains and then generalize to ChIP and GO
PlotEnrichment <- function(ClusterTable="",
                           Title="",
                           TopGenes=15){
  library(ggplot2)
  library(stringr)
  pal <- c("#74A9CF", "#3690C0", "#0570B0", "#045A8D", "#023858")
  # Read in cluster data file
  ClusterData <- read.table(ClusterTable,
                            sep="\t",
                            header=T,
                            stringsAsFactors = F)
  # Sort and filter gene number
    if(dim(ClusterData)[1]<TopGenes){
      warning("Fewer categories than selected.")
    }
  ClusterData <- na.omit(ClusterData[order(ClusterData$Enrichment_pval)[1:TopGenes],])
  # Make a shortened description
  ClusterData$ShortDescr <- paste0(str_trim(strtrim(gsub("\\(IPR[0-9]*\\)","", ClusterData$Description, perl=T), 30)),"...",substr(ClusterData$Description, nchar(ClusterData$Description)-9, nchar(ClusterData$Description)-1))
  # Make the Plot
  Plot <- ggplot(ClusterData, 
                 aes(x=reorder(ShortDescr, -Enrichment_pval), 
                     y=-log10(Enrichment_pval), 
                     fill=Hits)) +
    stat_summary(geom = "bar", fun = mean, position = "dodge") +
    xlab("Interpro Domain") +
    ylab("Log Fold Enrichment") +
    scale_fill_gradientn(colours = pal) +
    ggtitle(Title) +
    scale_y_continuous(breaks = round(seq(0, max(-log10(ClusterData$Enrichment_pval)), by = 2), 1)) +
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
  return(Plot)
}

# Nobt IPR
for ( i in as.numeric(gsub("\\D",
                           "",
                           list.files(path="DEGAnalysis/Pfam/Lists",
                                      pattern="^Nobt_IPR_Cluster_*")))) {
  PlotEnrichment(ClusterTable = paste0("DEGAnalysis/Pfam/Lists/Nobt_IPR_Cluster_",i,".txt"),
                 Title=paste0("Tobacco Cluster ",i))
  ggsave(filename=paste0("DEGAnalysis/Pfam/Plots/Nobt_IPR_Cluster_",i,".pdf"),
         width=12,
         height=8)
}

# TAIR IPR
for ( i in as.numeric(gsub("\\D",
                           "",
                           list.files(path="DEGAnalysis/Pfam/",
                                      pattern="^TAIR_Cluster_*")))) {
  if(file.size(paste0("DEGAnalysis/Pfam/TAIR_Cluster",i,".txt"))<10) {
    next
  }
  PlotEnrichment(ClusterTable = paste0("DEGAnalysis/Pfam/TAIR_Cluster",i,".txt"),
                 Title=paste0("Arabidopsis Cluster ",i))
  ggsave(filename=paste0("DEGAnalysis/Pfam/Plot_TAIR_Cluster",i,"_IPR_Enrichment.pdf"),
         width=12,
         height=8)
}

# Slyc IH IPR
for ( i in as.numeric(gsub("\\D",
                           "",
                           list.files(path="DEGAnalysis/Pfam/Lists",
                                      pattern="^Slyc_IPR_Cluster_*")))) {
  if(file.size(paste0("DEGAnalysis/Pfam/Lists/Slyc_IPR_Cluster_",i,".txt"))<10){
    next
  }
  PlotEnrichment(ClusterTable = paste0("DEGAnalysis/Pfam/Lists/Slyc_IPR_Cluster_",i,".txt"),
                 Title=paste0("Tomato Cluster ",i))
  ggsave(filename=paste0("DEGAnalysis/Pfam/Plots/Slyc_IPR_Cluster_",i,".pdf"),
         width=12,
         height=8)
}

# Slyc SRA IPR
for ( i in as.numeric(gsub("\\D",
                           "",
                           list.files(path="DEGAnalysis/Pfam/",
                                      pattern="^SlycSRA_Cluster_*")))) {
  if(file.size(paste0("DEGAnalysis/Pfam/SlycSRA_Cluster",i,".txt"))<10){
    next
  }
  PlotEnrichment(ClusterTable = paste0("DEGAnalysis/Pfam/SlycSRA_Cluster",i,".txt"),
                 Title=paste0("Tomato Cluster ",i))
  ggsave(filename=paste0("DEGAnalysis/Pfam/Plot_SlycSRA_Cluster",i,"_IPR_Enrichment.pdf"),
         width=12,
         height=8)
}

# Solanum IPR
for (i in c(1,5,7,8,10,12)) {
  PlotEnrichment(ClusterTable = paste0("DEGAnalysis/Pfam/Lists/Solanum_3DF_Noise_IPR_Cluster_",i,".txt"),
                 Title=paste0("Solanum Cluster ",i))
ggsave(filename=paste0("DEGAnalysis/Pfam/Plots/Solanum_3DF_Noise_Cluster_",i,"_IPR_Enrichment.pdf"),
         width=12,
         height=8)
}


