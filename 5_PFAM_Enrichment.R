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
    sorted_gene_groups[index] = paste(temp_genes, sep = ", ", collapse = ", ") 
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







