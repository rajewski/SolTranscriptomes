
#Read in the Sample list
metadata <- read.table("DEGAnalysis/SampleList.txt", header=T, sep="\t")
# Add Path to BAM file
metadata$Path <- NA
metadata$Path[metadata$Species=="Arabidopsis"] <- paste0("DEGAnalysis/STAR/TAIR10/",metadata$Accession[metadata$Species=="Arabidopsis"],".Aligned.sortedByCoord.out.bam")
metadata$Path[metadata$Species=="Tomato"] <- paste0("DEGAnalysis/STAR/Slyc/",metadata$Accession[metadata$Species=="Tomato"],".Aligned.sortedByCoord.out.bam")
metadata$Path[metadata$Species=="Tobacco"] <- paste0("DEGAnalysis/STAR/Nobt/",metadata$Accession[metadata$Species=="Tobacco"],".Aligned.sortedByCoord.out.bam")

#Read in BAM files
#BiocManager::install("Rsamtools")
library("Rsamtools")
bamfiles <- BamFileList(metadata$Path, yieldSize=2000000)
seqinfo(bamfiles[13]) #check that it worked

#BiocManager::install("GenomicFeatures")
library("GenomicFeatures")
Slyctxdb <- makeTxDbFromGFF("SlycDNA/ITAG4.0_gene_models.gff", organism="Solanum lycopersicum")

