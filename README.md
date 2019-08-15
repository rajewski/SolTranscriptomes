# slyc-WT

## Summary

This repo is for the analysis of the transcriptome data from tomato and tobacco using a pipeline consisting of STAR, StringTie, and DEseq2. Currently this repo is starting with the output of Stringtie that has been condensed using the prepDE.py script from StringTie:

```bash
module load stringtie/1.3.3b
prepDE.py -i Ballgown
```

In the case of tobacco, the transcript names caused an error and a modified version of the python script (included in this repo) was used along with a list of samples (also included for demonstration):

```bash
python prepDE.py  -i samplelist.txt
```

## Other Analysis Files

I do have plans to integrate the analysis scripts from STAR and StringTie, so stay tuned as I make this a more integrated pipeline. Additionally there are a handful of other scripts in this git that are legacy from lder versions of the analysis that I had been doing. These include:
* Other DE gene testing scripts:

 * Ballgown.R

   In place of DEseq2 I tried using Ballgown, which is part of the same package as StringTie. I found the analysis with Ballgown was not as thorough or reliable as the one with DEseq2, but I am keeping this Ballgown script around in case I decide to use it later.
   
 * sleuth.sh

   Formerly this pipeline was run with the output of Kallisto, but I have since switched to STAR and StringTie. Ultimately this decision was arbitrary and I don't think that one is superior to the other, but I also wanted to aling my analysis pipeline with that of my labmates for easier troubleshooting. In the Kallisto-based pipeline, sleuth was used in place of DEseq2 for differential expression testing. DEseq2 also allows for more flexible analysis of the genes. That being said, the sleuth scripts are useful sources of ideas for analyses.

* Clustering Scripts
   
 * glimma.R

   This script is rather poorly design (part of the reason I'm not using it anymore) but works in concert with sleuth.R to cluster the profiles of gene expression. This script aslo makes a heatmap of gene expression profiles, and this takes a long time to actually plot. Because this plotting is so intense, I have created two other scripts to do it as a separate job on a cluster computer.

 * heatmap.R
 
 This script is simply the subsection of code from glimma.R that makes a heatmap from a saved R data object. It is meant to be run as a submitted job on a cluster.
 
 * runr.sh
 
 This script is simply a submission script for heatmap.R.



index.sh is designed to create an transcriptome index from the ITAG3.0_cDNA.fasta file and put it in the data folder. I chose the cDNA instead of the CDS because it is possible that some of my reads could hit against a UTR and that sequence wouldn't be present in the CDS file. Importanly, I cleaned the transcript headers of this file to just contain the gene ID and not the putative function.

quant.sh is designed to work as an arrary job, taking the names of stages and replicates from the replist.txt file. If executed properly, this should spawn 15 subjobs that analyze each of the 5 stages in 3 reps and output them to a separte results folder. I have run this again with kallisto using 100 bootstraps and attempting to generate a genomebam file using the ITAG3.0_gene_models.gff file and the another TSV file of the tomato chromosomes along with their lengths (generated from a one-liner I found athttp://www.danielecook.com/generate-fasta-sequence-lengths/)