# slyc-WT

This repo is for the analyzis of the wildtype transcriptome data from tomato using kallisto.

index.sh is designed to create an transcriptome index from the ITAG3.0_cDNA.fasta file and put it in the data folder. I chose the cDNA instead of the CDS because it is possible that some of my reads could hit against a UTR and that sequence wouldn't be present in the CDS file. Importanly, I cleaned the transcript headers of this file to just contain the gene ID and not the putative function.

quant.sh is designed to work as an arrary job, taking the names of stages and replicates from the replist.txt file. If executed properly, this should spawn 15 subjobs that analyze each of the 5 stages in 3 reps and output them to a separte results folder. I have run this again with kallisto using 100 bootstraps and attempting to generate a genomebam file using the ITAG3.0_gene_models.gff file and the another TSV file of the tomato chromosomes along with their lengths (generated from a one-liner I found athttp://www.danielecook.com/generate-fasta-sequence-lengths/)