#!/bin/bash -l
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --nodes=1
#SBATCH --time=1:00:00
#SBATCH --mail-user=araje002@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH -o /bigdata/littlab/arajewski/FULTranscriptomes/logs/GetExternalData-%A.out
set -e

cd  /rhome/arajewski/bigdata/FULTranscriptomes/ExternalData/TAIR10

#Fetch Data
if [ ! -e TAIR10.fa ]; then
    echo Downloading Arabidopsis genome from TAIR...
    curl https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas > TAIR10.fa
    echo Done.
else
    echo Arabidopsis genome already present.
fi

if [ ! -e TAIR10.gff3 ]; then
    echo Downloading Arabidopsis genome annotation from TAIR...
    curl https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff > TAIR10.gff3
    echo Done.
else
    echo Arabidopsis genome annotation already present.
fi

if [ ! -e TAIR10.proteins.fa ]; then
    echo Downloading Arabiodpsis proteins from TAIR...
    curl https://www.arabidopsis.org/download_files/Proteins/TAIR10_protein_lists/TAIR10_pep_20101214 > TAIR10.proteins.fa
    echo Done.
else
    echo Arabidopsis proteins already present.
fi
