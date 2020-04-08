# You need to run the 5_PFAM_Enrichment.sh first because it generates input files for this script
# Again this script borrows heavily from the supplment of https://doi.org/10.1104/pp.108.132985


# Do a hypergeometric test for enrichment of each domain among top_genes 
## q: Genes with given domain that are in the top set are the white balls drawn from the urn.
## m: Genes with given domain are white balls.
## n: Genes without given domain are black balls.
## k: The balls drawn are the genes in the top set.
## loop over each domain 
PfamEnrichment <- function(AllGenesFile = "",
                           AllPfamFile = "",
                           PfamDescFile = "",
                           ExcludedGenesFile = "",
                           TopGenesFile = "",
                           OutputFile = "") {
  # Read in files to lists
  AllPfam <- scan(file=AllPfamFile,
                  what = list("",""),
                  sep = "\t")
  PfamDesc <- scan(file=PfamDescFile,
                   what = list("",""),
                   sep = "\t",
                   quote = "")
  ExcludedGenes <- scan(file=ExcludedGenesFile,
                        what = list(""))[[1]]
  AllGenes <- scan(file=AllGenesFile,
                   what=list(""))[[1]]
  top_genes <- scan(file=TopGenesFile,
                    what = list(""),
                    skip = 1)[[1]] 
  output_file = OutputFile
  # Clean lists
  AllGenes <- AllGenes[-c(AllGenes %in% ExcludedGenes)]
  top_genes <- top_genes[-c(top_genes %in% ExcludedGenes)]
  top_genes <- top_genes[top_genes %in% AllGenes]
  # Split list of domains
  domain_grouping <- split(AllPfam[[1]], AllPfam[[2]])
  # Set up hypergeometric test
  p_value_vector = c() 
  enrichment_vector = c() 
  q_vector = c()
  m_vector = c()
  # Do hypergeometric test
  for (index in 1:length(domain_grouping)){
    domain = names(domain_grouping)[index] 
    genes_with_domain = domain_grouping[[index]]
    q = length(as.vector(na.omit(match(genes_with_domain, top_genes)))) 
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
}



## adjust p-values
library(multtest)
adjusted_p <- mt.rawp2adjp(p_value_vector, proc = c("BH"))

## write result
p_cut_off = 0.1
significant_indices = as.vector(na.omit(ifelse(adjusted_p$adjp[,2] < p_cut_off, adjusted_p$index, NA)))
significant_enrichments = enrichment_vector[significant_indices] 
indices_by_enrichment = significant_indices[sort(significant_enrichments, decreasing = TRUE, index.return = TRUE)[[2]]]
sorted_descriptions = PfamDesc[[2]][match(names(domain_grouping)[indices_by_enrichment], PfamDesc[[1]])]
sorted_Pfam_ID = names(domain_grouping)[indices_by_enrichment] 
sorted_full_descriptions = paste(sorted_descriptions, " (", sorted_Pfam_ID, ")", sep = "") 
sorted_enrichment = enrichment_vector[indices_by_enrichment]
sorted_p_values = adjusted_p$adjp[match(indices_by_enrichment, adjusted_p$index),2] 
sorted_q = q_vector[indices_by_enrichment]
sorted_m = m_vector[indices_by_enrichment]

sorted_gene_groups = c()
for (index in 1:length(indices_by_enrichment)){
  temp_genes = top_genes[as.vector(na.omit(match(domain_grouping[[indices_by_enrichment[index]]], top_genes)))]
  sorted_gene_groups[index] = paste(temp_genes, sep = ", ", collapse = ", ") }

output = cbind(sorted_full_descriptions, round(sorted_enrichment, digits = 1), as.character(format(signif(sorted_p_values, digits = 1), scientific = TRUE)), sorted_q, sorted_m, sorted_gene_groups)
write(t(output), file = output_file, sep = "\t", ncolumns = 6)
