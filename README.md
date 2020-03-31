# Summary

In this repository I am conducting a DEG analysis of *Solanum lycopersicum* cv. Ailsa Craig, *Nicotiana obtusifolia*, and *Arabidopsis thaliana*. I will supplement this DEG analysis with publicly available ChIP-chip and ChIP-seq data from *S. lycopersicum* and *A. thaliana* as well as publicly available RNA-seq data from *FRUITFULL* mutants of both speices.

# Methods

## DEG Analysis

For RNA-seq, this will consist of STAR as the read aligner and DEseq2 as the differential expression testing software. For the microarray analyses, as a preprocessing step, the probe sequences were mapped back to transcriptomic targets, and nonmapping or multimapping probes were excluded from the analyes. Differential expression testing will be done with limma/voom. In all caases, the data are fit to a spline regression model with one fewer degree of freedom than there are timepoints in the specific study. (See supplement of [Sander et al, 2017](https://www.ncbi.nlm.nih.gov/pubmed/27797772) for an implementation of this with DESeq2.)

### In-House Data
 - I am working on uploading these to NCBI so that the whole thing will be with publicly available data

Our lab has generated several transcriptome libraries that will be used for this analysis. They are symlinked to in the `SlycRNA` and `NobtRNA` for tomato and tobacco, respecitively

### Public Data

#### Gene Expression Data

The RNA seq files are downloaded from SRA with the `1_GetExternalSRAData.sh` script and placed in the `ExternalData/RNAseq` directory.

| Species | NCBI BioProject ID | Description |
| ------- | ------------------ | ----------- |
| Tomato | [PRJNA213528](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA213528)  | RNA-seq from a developmental series of fruits using wild-type and *ful1*/*ful2* knockdown lines. From Fujusawa et al, 2014 |
| Tomato | [PRJNA177429](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA177429) | Microarray from a developmental series of fruits using wild-type and *ful1*/*ful2* knockdown lines. From Bemer et al, 2010 |
| Arabidopsis | [PRJEB25745](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB25745) | RNA-seq of wild-type fruit valve tissue. From Mizzotti et al, 2018 | 
| Arabidopsis | [PRJNA316153](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA316153) | Microarray of *ful-1* FUL-GR siliques under both mock and DEX treatment. From Bemer et al, 2017 |

#### Genomes

I am currently using the [SL4.0 genome and the ITAG4.0 annotation for tomato](https://solgenomics.net/organism/Solanum_lycopersicum/genome). For tobacco (*N. obtusifolia*) I am using the publicly available [genome and annotation](http://nadh.ice.mpg.de/NaDH/download/overview) but with a slight modification to add in the the MBP20 gene, which was missed in the original annotation. For *Arabidopsis* I'll be using the TAIR10 genome and annotation. In the case of tomato and tobacco these files are symlinked to copies we already have using the `SlyDNA` and `NobtDNA` directories, respecitively. 

For Arabidopsis, the data is downloaded with the `1_GetExternalGenome.sh` to `ExternalData/TAIR10`, but I have excluded those files from the git to save space.

# ChIP Analysis

As a preprocessing step for the microarray study, array probe sequences were mapped to the SL4.0 genome, and probe locations were updated. Nonmapping or multimapping probes were excluded from the analysis.


- How tf am I going to do a chip analysis from raw data?
- How to integrate chipseq with chip-chip? a simple list of enriched regions/promoters

### Public Data

| Species | NCBI BioProject ID | Description |
| ------- | ------------------ | ----------- |
| Tomato | [PRJNA213106](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA213106) | ChIP-chip of wild-type fruits using anti-FUL1 or anti-FUL2 antibodies. From Fujisawa et al, 2014 |
| Arabidopsis | [PRJNA316152](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA316152) | ChIP-seq of pistil/siliques of *ful-1 pFUL*:*FUL-GFP* using an anti-GFP antibody. From Bemer et al, 2017 |
| Arabidopsis | [PRJNA427320](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA427320) | ChIP-seq of pistil/siliques of *ful-1 pFUL*:*FUL-GFP* using an anti-GFP antibody. From Balanza et al 2018 |



I am editing a new branch of this repo as I go and will add my new methods as nedded

module load stringtie/1.3.3b
prepDE.py -i Ballgown
```

In the case of tobacco, the transcript names caused an error and a modified version of the python script (included in this repo) was used along with a list of samples (also included for demonstration):

```bash
python prepDE.py  -i samplelist.txt
```

## Differential Expression Testing

The main workhorse script in this git is `DEseq2.R`. It takes as input the count matrix generated by the python script mentioned earlier. The script also creates a design matrix that has the conditions for each sample. The statistical analysis is straight forward and has relatively few options. Currently the cutoff for DE is a log fold change threshold of 1 and an alpha of 0.05. This script then outputs two files for each species being tested. One for genes that are upregulated and another for genes that are downregulated. The up or down determination is made based on setting one of the conditions as the arbitrary baseline. Unfortunately it seems to alphabetize them, which is not always the best solution. The direction of the comparison can be found using the `summary(res)` command.

The output .csv files from this script can then be fed into the [GO Enrichment git](https://github.com/rajewski/GoAnalysis) I have.

## Orthology

Because my project is comparing two difeferent species that do not share a common transcriptome for RNA-seq read mapping, I have to find some way to compare gene expression across species. There are a couple ways I've though of to do this, but currently I am restricting myself to looking at 1:1 orthologous genes between the two species. This certianly leaves out many many genes that are important, so I am trying to polish this analysis to get more power.

This section of the repo is meant to be run at the same time as the [GO Enrichment pipeline](https://github.com/rajewski/GoAnalysis), and thus the `ReciprocalHits.sh` script takes as input some files generated in that pipeline. Note that any list of protein sequences could be used as an input though.

As is, the `ReciprocalHits.sh` script is what I have been using to find 1:1 orthologous genes. The search is based on the best reciprocal hits method. To determine orthology with this method, the protein sequence of a given gene from species A is searched against the proteome of species B. The top hit from this search is then used as the query for a second search against the proteome of species A. If the gene from the original query matches the top hit of the second search, then the genes are classified as 1:1 orthologous. 

The `ReciprocalHits.sh` script creates a custom DIAMOND database (instead of BLAST) for each species, and then conducts the search for each up- or down-regulated gene used as an input. The input protein sequence file (e.g. NobtUp.fasta) is used as the query to create an output file of best hits (e.g. NobtUptoSlyc.dmndo) against the other species' proteome.. The input protein sequence file can be generated during the GO Enrichment pipeline mentioned earlier. Using a couple of bash and awk one liners, the .dmndo is parsed into a fasta file of protein sequences for the top hits (e.g. NobtUptoSlycHits.fasta) and a list of query-hit pairs (e.g. NobtUptoSlycHits.tsv).

The fasta file of protein sequences for the top hits is then used as a query list for a search back against original species' proteome. This generates a DIAMOND output file (e.g. NobtUptoNobt.dmndo), which is then parsed with the same bash one-liner to generate a list of query-hit pairs (e.g. NobtUptoNobtHits.tsv).

The next step is done in R with `Orthology.R`. This script takes the list of query-hit pairs from both searches (e.g. NobtUptoSlycHits.tsv and NobtUptoNobtHits.tsv), cleans the names, merges the two lists, and then asks if the original query gene is the same as the final hit gene. This merged list is then aggregated by the gene from species B that acts as a "bridge" between the original query and the final hit. If any gene from species B every acts as a bridge between mismatched query/hit genes, it is thrown out. This step eliminates genes with questionable orthology. Later any gene from species B that acts as a bridge between multiple query/hit genes is also removed. This step eliminates genes that are not 1:1 orthologs.

This orthology search is repeated for up- and down-regulated genes in each species.The script then searches through the lists of DEGs to find genes in both sets that are present in 1:1 orthology and writes them to a .csv file (e.g. ConservedUp.csv). This csv file can then be used for a separate GO enrichment analysis.

## Other Analysis Files

I do have plans to integrate the analysis scripts from STAR and StringTie, so stay tuned as I make this a more integrated pipeline. Additionally there are a handful of other scripts in this git that are legacy from lder versions of the analysis that I had been doing. These include:
### Other DE gene testing scripts:

 * `Ballgown.R`

   In place of DEseq2 I tried using Ballgown, which is part of the same package as StringTie. I found the analysis with Ballgown was not as thorough or reliable as the one with DEseq2, but I am keeping this Ballgown script around in case I decide to use it later.
   
 * `sleuth.sh`

   Formerly this pipeline was run with the output of Kallisto, but I have since switched to STAR and StringTie. Ultimately this decision was arbitrary and I don't think that one is superior to the other, but I also wanted to aling my analysis pipeline with that of my labmates for easier troubleshooting. In the Kallisto-based pipeline, sleuth was used in place of DEseq2 for differential expression testing. DEseq2 also allows for more flexible analysis of the genes. That being said, the sleuth scripts are useful sources of ideas for analyses.

### Clustering Scripts
   
 * `glimma.R`

   This script is rather poorly design (part of the reason I'm not using it anymore) but works in concert with `sleuth.R` to cluster the profiles of gene expression. This script aslo makes a heatmap of gene expression profiles, and this takes a long time to actually plot. Because this plotting is so intense, I have created two other scripts to do it as a separate job on a cluster computer.

 * `heatmap.R`
 
 This script is simply the subsection of code from glimma.R that makes a heatmap from a saved R data object. It is meant to be run as a submitted job on a cluster.
 
 * `runr.sh`
 
 This script is simply a submission script for heatmap.R.
 
 * `WGCNA.R`
 
 In place of the glimma.R script, I also tried using WGCNA for clustering. I have some philosophical problems with WGCNA, but those are mostly moot.
 
 * `genie.R`
 
 This is technically more of a network creation script, but it operates off of clustering output files. It creates directed networks of genes in a cluster and outputs pdf graphs of the clusters.


