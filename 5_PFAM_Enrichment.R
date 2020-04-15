# You need to run the 5_PFAM_Enrichment.sh first because it generates input files for this script
# Again this script borrows heavily from the supplment of https://doi.org/10.1104/pp.108.132985


PfamEnrichment <- function(AllGenesFile = "",
                           AllIPRFile = "",
                           IPRDescFile = "",
                           ExcludedGenesFile = "",
                           TopGenesFile = "",
                           OutputFile = "",
                           p_cut_off = 0.1,
                           WriteToFile=TRUE) {
  # Define a useful function
  `%notin%` <- Negate(`%in%`)
  # Read in files to lists
  AllIPR <- scan(file=AllIPRFile,
                  what = list("",""),
                  sep = "\t")
  IPRDesc <- scan(file=IPRDescFile,
                   what = list("",""),
                   sep = "\t",
                   quote = "")
  ExcludedGenes <- scan(file=ExcludedGenesFile,
                        what = list(""))[[1]]
  AllGenes <- scan(file=AllGenesFile,
                   what=list(""))[[1]]
  top_genes <- scan(file=TopGenesFile,
                    what = list(""))[[1]] 
  output_file = OutputFile
  # Clean lists
  AllGenes <- AllGenes[AllGenes %notin% ExcludedGenes]
  top_genes <- top_genes[top_genes %notin% ExcludedGenes]
  top_genes <- top_genes[top_genes %in% AllGenes]
  # Split list of domains
  domain_grouping <- split(AllIPR[[1]], AllIPR[[2]])
  # Set up hypergeometric test
  ## q: Genes with given domain that are in the top set are the white balls drawn from the urn.
  ## m: Genes with given domain are white balls.
  ## n: Genes without given domain are black balls.
  ## k: The balls drawn are the genes in the top set.
  p_value_vector = c() 
  enrichment_vector = c() 
  q_vector = c()
  m_vector = c()
  # Do hypergeometric test
  for (index in 1:length(domain_grouping)){
    domain = names(domain_grouping)[index] 
    genes_with_domain = domain_grouping[[index]]
    q = sum(genes_with_domain %in% top_genes)
    m = length(genes_with_domain)
    n = length(AllGenes) - m
    k = length(top_genes)
    p_value = phyper(q - 1, m, n, k, lower.tail = FALSE) 
    p_value_vector[index] = p_value
    enrichment = q/(m*k/(n+m)) 
    enrichment_vector[index] = enrichment 
    q_vector[index] = q
    m_vector[index] = m
  }
  # Adjust p-values
  library(multtest)
  adjusted_p <- mt.rawp2adjp(p_value_vector, proc = c("BH"))
  # Write result
  significant_indices = adjusted_p$index[adjusted_p$adjp[,2] < p_cut_off]
  if (length(significant_indices) == 0 ) {
    warning("No significantly enriched domains in this group. Skipping")
    if (WriteToFile){
        write.table(NULL, file = output_file, sep = "\t", quote=FALSE, row.names = FALSE)
    }
    return(NULL)
  }
  significant_enrichments = enrichment_vector[significant_indices] 
  indices_by_enrichment = significant_indices[sort(significant_enrichments, decreasing = TRUE, index.return = TRUE)[[2]]]
  sorted_descriptions = IPRDesc[[2]][match(names(domain_grouping)[indices_by_enrichment], IPRDesc[[1]])]
  sorted_IPR_ID = names(domain_grouping)[indices_by_enrichment] 
  sorted_full_descriptions = paste0(sorted_descriptions, " (", sorted_IPR_ID, ")") 
  sorted_enrichment = enrichment_vector[indices_by_enrichment]
  sorted_p_values = adjusted_p$adjp[match(indices_by_enrichment, adjusted_p$index),2] 
  sorted_q = q_vector[indices_by_enrichment]
  sorted_m = m_vector[indices_by_enrichment]
  sorted_gene_groups = c()
  for (index in 1:length(indices_by_enrichment)){
    temp_genes = top_genes[as.vector(na.omit(match(domain_grouping[[indices_by_enrichment[index]]], top_genes)))]
    sorted_gene_groups[index] = paste(temp_genes, sep = ",", collapse = ",") 
  }
  output = cbind(Description=sorted_full_descriptions, 
                 Enrichment=round(sorted_enrichment, digits = 1), 
                 Enrichment_pval=as.character(format(signif(sorted_p_values, digits = 1), scientific = TRUE)), 
                 Hits=sorted_q, 
                 Total=sorted_m, 
                 Gene_IDs=sorted_gene_groups)
  if(WriteToFile) {
  write.table(output, file = output_file, sep = "\t", quote=FALSE, row.names = FALSE)
  }
 return(output)
}


# IPR on RNA seq clusters -------------------------------------------------

# Work through the Slyc IH data
for ( i in 1:length(list.files(path="DEGAnalysis/RNA-seq/", pattern="SlycIH_All_Cluster*"))) {
  PfamEnrichment(AllGenesFile = "DEGAnalysis/Pfam/Slyc.protein.names.txt",
                 AllIPRFile = "DEGAnalysis/Pfam/Slyc.gene2ipr.tsv",
                 IPRDescFile = "DEGAnalysis/Pfam/Slyc.ipr2desc.tsv",
                 ExcludedGenesFile = "DEGAnalysis/Pfam/Slyc.nopfam.tsv",
                 TopGenesFile = paste0("DEGAnalysis/RNA-seq/SlycIH_All_Cluster_", i, ".txt"),
                 OutputFile = paste0("DEGAnalysis/Pfam/SlycIH_All_Cluster", i, ".txt"))
}

# Work through the Slyc SRA data
for ( i in 1:length(list.files(path="DEGAnalysis/RNA-seq/", pattern="SlycSRA_All_Cluster*"))) {
  PfamEnrichment(AllGenesFile = "DEGAnalysis/Pfam/Slyc.protein.names.txt",
                 AllIPRFile = "DEGAnalysis/Pfam/Slyc.gene2ipr.tsv",
                 IPRDescFile = "DEGAnalysis/Pfam/Slyc.ipr2desc.tsv",
                 ExcludedGenesFile = "DEGAnalysis/Pfam/Slyc.nopfam.tsv",
                 TopGenesFile = paste0("DEGAnalysis/RNA-seq/SlycSRA_All_Cluster_", i, ".txt"),
                 OutputFile = paste0("DEGAnalysis/Pfam/SlycSRA_All_Cluster", i, ".txt"))
}

# Work through the TAIR data
for ( i in 1:length(list.files(path="DEGAnalysis/RNA-seq/", pattern="TAIR10_Cluster*"))) {
  PfamEnrichment(AllGenesFile = "DEGAnalysis/Pfam/TAIR10.protein.names.txt",
                 AllIPRFile = "DEGAnalysis/Pfam/TAIR10.gene2ipr.tsv",
                 IPRDescFile = "DEGAnalysis/Pfam/TAIR10.ipr2desc.tsv",
                 ExcludedGenesFile = "DEGAnalysis/Pfam/TAIR10.nopfam.tsv",
                 TopGenesFile = paste0("DEGAnalysis/RNA-seq/TAIR10_Cluster_", i, ".txt"),
                 OutputFile = paste0("DEGAnalysis/Pfam/TAIR10_All_Cluster", i, ".txt"))
}

# And for Nobt
for ( i in 1:length(list.files(path="DEGAnalysis/RNA-seq/", pattern="Nobt_All_Cluster*"))) {
  PfamEnrichment(AllGenesFile = "DEGAnalysis/Pfam/Nobt.protein.names.txt",
                 AllIPRFile = "DEGAnalysis/Pfam/Nobt.gene2ipr.tsv",
                 IPRDescFile = "DEGAnalysis/Pfam/Nobt.ipr2desc.tsv",
                 ExcludedGenesFile = "DEGAnalysis/Pfam/Nobt.nopfam.tsv",
                 TopGenesFile = paste0("DEGAnalysis/RNA-seq/Nobt_All_Cluster_", i, ".txt"),
                 OutputFile = paste0("DEGAnalysis/Pfam/Nobt_All_Cluster", i, ".txt"))
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
# Table of transcripts and binned antibody hits
write.table(FULBinding[,c("Transcript", "Bins")],
            file="ChIPAnalysis/ChIP-chip/TomatoFULBinnedTargets.csv",
            row.names = F,
            col.names= F,
            quote = F,
            sep="\t")
# Descriptions of the possible hits
write.table(rbind(c("FUL1.sm", "FUL1_bound"),
                  c("FUL2.sm", "FUL2_bound"),
                  c("FUL1_50", "FUL1 binding between 0 and 50bp of TSS"),
                  c("FUL1_100", "FUL1 binding between 510 and 100bp of TSS"),
                  c("FUL1_200", "FUL1 binding between 101 and 200bp of TSS"),
                  c("FUL1_500", "FUL1 binding between 201 and 500bp of TSS"),
                  c("FUL1_1000", "FUL1 binding between 501 and 1000bp of TSS"),
                  c("FUL1_1500", "FUL1 binding between 1001 and 1500bp of TSS"),
                  c("FUL1_2000", "FUL1 binding between 1501 and 2000bp of TSS"),
                  c("FUL2_50", "FUL2 binding between 0 and 50bp of TSS"),
                  c("FUL2_100", "FUL2 binding between 510 and 100bp of TSS"),
                  c("FUL2_200", "FUL2 binding between 101 and 200bp of TSS"),
                  c("FUL2_500", "FUL2 binding between 201 and 500bp of TSS"),
                  c("FUL2_1000", "FUL2 binding between 501 and 1000bp of TSS"),
                  c("FUL2_1500", "FUL2 binding between 1001 and 1500bp of TSS"),
                  c("FUL2_2000", "FUL2 binding between 1501 and 2000bp of TSS")),
            file="ChIPAnalysis/ChIP-chip/TomatoFULdescription.csv",
            row.names = F,
            col.names= F,
            quote = F,
            sep="\t")
# Table of the zero genes to exclude
write.table(NULL, file="ChIPAnalysis/ChIP-chip/DummyChip.txt")


# Slyc IH data on general antibody binding
for ( i in 1:length(list.files(path="DEGAnalysis/RNA-seq/", pattern="SlycIH_All_Cluster*"))) {
  PfamEnrichment(AllGenesFile = "DEGAnalysis/Pfam/Slyc.protein.names.txt",
                 AllIPRFile = "ChIPAnalysis/ChIP-chip/TomatoFULTargets.csv",
                 IPRDescFile = "ChIPAnalysis/ChIP-chip/TomatoFULdescription.csv",
                 ExcludedGenesFile = "ChIPAnalysis/ChIP-chip/DummyChip.txt",
                 TopGenesFile = paste0("DEGAnalysis/RNA-seq/SlycIH_All_Cluster_", i, ".txt"),
                 OutputFile = paste0("ChIPAnalysis/ChIP-chip/SlycIH_All_Cluster", i, ".txt"))
}

# Slyc SRA data on general antibody binding
for ( i in 1:length(list.files(path="DEGAnalysis/RNA-seq/", pattern="SlycSRA_All_Cluster*"))) {
  PfamEnrichment(AllGenesFile = "DEGAnalysis/Pfam/Slyc.protein.names.txt",
                 AllIPRFile = "ChIPAnalysis/ChIP-chip/TomatoFULTargets.csv",
                 IPRDescFile = "ChIPAnalysis/ChIP-chip/TomatoFULdescription.csv",
                 ExcludedGenesFile = "ChIPAnalysis/ChIP-chip/DummyChip.txt",
                 TopGenesFile = paste0("DEGAnalysis/RNA-seq/SlycSRA_All_Cluster_", i, ".txt"),
                 OutputFile = paste0("ChIPAnalysis/ChIP-chip/SlycSRA_All_Cluster", i, ".txt"))
}

# Slyc IH data on binned antibody binding
for ( i in 1:length(list.files(path="DEGAnalysis/RNA-seq/", pattern="SlycIH_All_Cluster*"))) {
  PfamEnrichment(AllGenesFile = "DEGAnalysis/Pfam/Slyc.protein.names.txt",
                 AllIPRFile = "ChIPAnalysis/ChIP-chip/TomatoFULBinnedTargets.csv",
                 IPRDescFile = "ChIPAnalysis/ChIP-chip/TomatoFULdescription.csv",
                 ExcludedGenesFile = "ChIPAnalysis/ChIP-chip/DummyChip.txt",
                 TopGenesFile = paste0("DEGAnalysis/RNA-seq/SlycIH_All_Cluster_", i, ".txt"),
                 OutputFile = paste0("ChIPAnalysis/ChIP-chip/SlycIH_All_Binned_Cluster", i, ".txt"))
}

# Slyc SRA data on binned antibody binding
for ( i in 1:length(list.files(path="DEGAnalysis/RNA-seq/", pattern="SlycSRA_All_Cluster*"))) {
  PfamEnrichment(AllGenesFile = "DEGAnalysis/Pfam/Slyc.protein.names.txt",
                 AllIPRFile = "ChIPAnalysis/ChIP-chip/TomatoFULBinnedTargets.csv",
                 IPRDescFile = "ChIPAnalysis/ChIP-chip/TomatoFULdescription.csv",
                 ExcludedGenesFile = "ChIPAnalysis/ChIP-chip/DummyChip.txt",
                 TopGenesFile = paste0("DEGAnalysis/RNA-seq/SlycSRA_All_Cluster_", i, ".txt"),
                 OutputFile = paste0("ChIPAnalysis/ChIP-chip/SlycSRA_All_Binned_Cluster", i, ".txt"))
}


enriched <- read.table("ChIPAnalysis/ChIP-chip/SlycIH_All_Binned_Cluster3.txt", sep="\t")

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

for (i in 1:length(list.files(path="ChIPAnalysis/ChIP-chip/", pattern="SlycIH_All_Cluster*"))) {
  if(file.size(paste0("ChIPAnalysis/ChIP-chip/SlycIH_All_Cluster",i, ".txt"))<10){
    next
  }
  pdf(file=paste0("ChIPAnalysis/ChIP-chip/Plot_SlycIH_ChIPHisto_Cluster", i,".pdf"),
      height=6,
      width=11)
  PlotChIPHist(ClusterFile = paste0("ChIPAnalysis/ChIP-chip/SlycIH_All_Cluster",i, ".txt"))
  dev.off()
}








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
for ( i in 1:length(list.files(path="DEGAnalysis/Pfam/", pattern="Nobt_All_Cluster*"))) {
  PlotEnrichment(ClusterTable = paste0("DEGAnalysis/Pfam/Nobt_All_Cluster",i,".txt"),
                 Title=paste0("Tobacco Cluster ",i))
  ggsave(filename=paste0("DEGAnalysis/Pfam/Plot_Nobt_Cluster",i,"_IPR_Enrichment.pdf"),
         width=12,
         height=8)
}

# TAIR IPR
for ( i in 1:length(list.files(path="DEGAnalysis/Pfam/", pattern="TAIR10_All_Cluster*"))) {
  if(file.size(paste0("DEGAnalysis/Pfam/TAIR10_All_Cluster",i,".txt"))<10){
    next
  }
  PlotEnrichment(ClusterTable = paste0("DEGAnalysis/Pfam/TAIR10_All_Cluster",i,".txt"),
                 Title=paste0("Arabidopsis Cluster ",i))
  ggsave(filename=paste0("DEGAnalysis/Pfam/Plot_TAIR_Cluster",i,"_IPR_Enrichment.pdf"),
         width=12,
         height=8)
}

# Slyc IH IPR
for ( i in 1:length(list.files(path="DEGAnalysis/Pfam/", pattern="SlycIH_All_Cluster*"))) {
  if(file.size(paste0("DEGAnalysis/Pfam/SlycIH_All_Cluster",i,".txt"))<10){
    next
  }
  PlotEnrichment(ClusterTable = paste0("DEGAnalysis/Pfam/SlycIH_All_Cluster",i,".txt"),
                 Title=paste0("Tomato Cluster ",i))
  ggsave(filename=paste0("DEGAnalysis/Pfam/Plot_SlycIH_Cluster",i,"_IPR_Enrichment.pdf"),
         width=12,
         height=8)
}

# Slyc SRA IPR
for ( i in 1:length(list.files(path="DEGAnalysis/Pfam/", pattern="SlycSRA_All_Cluster*"))) {
  if(file.size(paste0("DEGAnalysis/Pfam/SlycSRA_All_Cluster",i,".txt"))<10){
    next
  }
  PlotEnrichment(ClusterTable = paste0("DEGAnalysis/Pfam/SlycSRA_All_Cluster",i,".txt"),
                 Title=paste0("Tomato Cluster ",i))
  ggsave(filename=paste0("DEGAnalysis/Pfam/Plot_SlycSRA_Cluster",i,"_IPR_Enrichment.pdf"),
         width=12,
         height=8)
}

# Try the function with ChIP data of Slyc IH Binned
#it works but messes up the labels a bit and might not be the best way to plot
for ( i in 1:length(list.files(path="ChIPAnalysis/ChIP-chip/", pattern="SlycIH_All_Binned_Cluster*"))) {
  if(file.size(paste0("ChIPAnalysis/ChIP-chip/SlycIH_All_Binned_Cluster",i,".txt"))<10){
    next
  }
  PlotEnrichment(ClusterTable = paste0("ChIPAnalysis/ChIP-chip/SlycIH_All_Binned_Cluster",i,".txt"),
                 Title=paste0("Tomato Cluster ",i))
  ggsave(filename=paste0("ChIPAnalysis/ChIP-chip/Plot_SlycIH_Cluster",i,"_ChIP_Enrichment.pdf"),
         width=12,
         height=8)
}


