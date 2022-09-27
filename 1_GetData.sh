#!/bin/bash -l
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --nodes=1
#SBATCH --time=1:00:00
#SBATCH --mail-user=araje002@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH -o ./logs/GetExternalData-%A_%a.out
set -e

# Define a list of necessary external data sources
# This  should be run as an array job with indices specified below
# Index 0 for Arabidopsis data
# Index 1-3 for microarray data # Removed code, but accessions remain in list
# Index 4-51 are RNAseq
# Index  52-58 are for ChIP data # Removed code, but accessions remain in list
# Index 59 is for the melon data
AllExternal=(TAIR10 GSE41560 GSE49125 GSE79553 SRR943813 SRR943814 SRR943815 SRR943816 SRR943817 SRR943818 SRR943825 SRR943826 SRR943827 SRR943828 SRR943829 SRR943830 ERR2809794 ERR2809795 ERR2809796 ERR2809797 ERR2809798 ERR2809799 ERR2809800 ERR2809801 ERR2809802 ERR2809803 ERR2809804 ERR2809805 ERR2809806 ERR2809807 ERR2809808 ERR2809809 ERR2809810 ERR2809811 ERR2809812 ERR2809813 ERR2809814 ERR2809815 ERR2809816 ERR2809817 SRR3199616 SRR3199617 SRR3199618 SRR3199656 SRR3199657 SRR3199622 SRR3199623 SRR3199624 SRR3199628 SRR3199629 SRR3199630 SRR3199635 SRR6412402 SRR6412403 SRR6412404 SRR3288009 SRR3288010 SRR3288011 SRR3288012 Melon )

####### Download Arabidopsis Data
if [ "$SLURM_ARRAY_TASK_ID" == 0 ]; then
  cd  ExternalData/${AllExternal[$SLURM_ARRAY_TASK_ID]}
  # Get genome fasta
  if [ ! -e ${AllExternal[$SLURM_ARRAY_TASK_ID]}.fa ]; then
      echo Downloading ${AllExternal[$SLURM_ARRAY_TASK_ID]} genome...
      curl https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas > ${AllExternal[$SLURM_ARRAY_TASK_ID]}.fa
      echo Done.
  else
      echo ${AllExternal[$SLURM_ARRAY_TASK_ID]} genome already present.
  fi
  # Get annotation gff
  if [ ! -e ${AllExternal[$SLURM_ARRAY_TASK_ID]}.gff3 ]; then
      echo Downloading ${AllExternal[$SLURM_ARRAY_TASK_ID]} genome annotation...
      curl https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff > ${AllExternal[$SLURM_ARRAY_TASK_ID]}.gff3
      echo Done.
  else
      echo ${AllExternal[$SLURM_ARRAY_TASK_ID]} genome annotation already present.
  fi
  # Get protein fasta
  if [ ! -e ${AllExternal[$SLURM_ARRAY_TASK_ID]}.proteins.fa ]; then
      echo Downloading ${AllExternal[$SLURM_ARRAY_TASK_ID]} proteins...
      curl https://www.arabidopsis.org/download_files/Proteins/TAIR10_protein_lists/TAIR10_pep_20110103_representative_gene_model > ${AllExternal[$SLURM_ARRAY_TASK_ID]}.proteins.fa
      echo Done.
  else
      echo ${AllExternal[$SLURM_ARRAY_TASK_ID]} proteins already present.
  fi
  cd ../../
fi

####### Array IDs 1-3 are removed as obsolete expts


####### Download SRA RNAseq Data
if [ "$SLURM_ARRAY_TASK_ID" -ge 4 ] && [ "$SLURM_ARRAY_TASK_ID" -le 51 ]; then
  mkdir -p ExternalData/RNAseq
  cd ExternalData/RNAseq
  # Prefetch Data
  if [ ! -d ${AllExternal[$SLURM_ARRAY_TASK_ID]} ]; then
      echo Downloading ${AllExternal[$SLURM_ARRAY_TASK_ID]} data from SRA...
      module load sratoolkit/2.10.0
      prefetch ${AllExternal[$SLURM_ARRAY_TASK_ID]}
      echo Done.
  else
      echo ${AllExternal[$SLURM_ARRAY_TASK_ID]} data already present.
  fi
  # Dump Fastq files for data
  if [ ! -e ${AllExternal[$SLURM_ARRAY_TASK_ID]}_1.fastq.gz ]; then
      echo Dumping fastq data for ${AllExternal[$SLURM_ARRAY_TASK_ID]}...
      module load sratoolkit/2.10.0
      fastq-dump --defline-seq '@$sn[_$rn]/$ri' --defline-qual '+$sn[_$rn]/$ri' --split-files --gzip -B ${AllExternal[$SLURM_ARRAY_TASK_ID]}
      echo Done.
  else
      echo Fastq data for ${AllExternal[$SLURM_ARRAY_TASK_ID]} already present.
  fi
  #Trim Reads
  if [ ! -e ${AllExternal[$SLURM_ARRAY_TASK_ID]}_1_trimmed.fq.gz ]; then
      echo Running Trim Galore on ${AllExternal[$SLURM_ARRAY_TASK_ID]}...
      module unload perl
      module swap miniconda2 miniconda3
      module swap python/2.7.5 python/3.6.0
      module load fastqc
      conda activate cutadaptenv
      export PATH=/bigdata/littlab/arajewski/Datura/software/TrimGalore-0.6.5:$PATH
      trim_galore \
	  --no_report_file \
	  -j $SLURM_CPUS_PER_TASK \
          ${AllExternal[$SLURM_ARRAY_TASK_ID]}_1.fastq.gz
      echo Done.
  else
      echo ${AllExternal[$SLURM_ARRAY_TASK_ID]} RNA seq already trimmed.
  fi
  cd ../../
fi

####### Array IDs 52-58 removed as obsolete

# Melon Genome and RNA seq data
if [ "$SLURM_ARRAY_TASK_ID" == 59 ]; then
    mkdir -p ExternalData/C_melo
    cd  ExternalData/C_melo
  # Get genome fasta
  if [ ! -e CM3.5.1_genome.fa ]; then
      echo Downloading C. melo genome fasta genome...
      curl ftp://cucurbitgenomics.org/pub/cucurbit/genome/melon/v3.5.1/CM3.5.1_genome.fa.gz > CM3.5.1_genome.fa.gz
      echo Done downloading, now unzipping...
      gunzip CM3.5.1_genome.fa.gz
      echo Done.
  else
      echo C. melo genome already present.
  fi
  # Get annotation gff                                                                                                                                                           
  if [ ! -e CM3.5.1_gene.gff ]; then
      echo Downloading C. melo genome annotation...
      curl ftp://cucurbitgenomics.org/pub/cucurbit/genome/melon/v3.5.1/CM3.5.1_gene.gff.gz > CM3.5.1_gene.gff.gz
      echo Done downloading, now unzipping...
      gunzip CM3.5.1_gene.gff.gz
      echo Done.
  else
      echo C. melo genome annotation already present.
  fi
  # Get protein fasta
  if [ ! -e CM3.5.1_protein.fa ]; then
      echo Downloading C. melo proteins...
      curl ftp://cucurbitgenomics.org/pub/cucurbit/genome/melon/v3.5.1/CM3.5.1_protein.fa.gz > CM3.5.1_protein.fa.gz
      echo Done downlaoding, now unzipping...
      gunzip CM3.5.1_protein.fa.gz
      echo Done.
  else
      echo C. melo proteins already present.
  fi
  # Get GO Terms
  if [ ! -e CM3.5.1_GO_anno.txt ]; then
      echo Downloading C. melo GO Terms...
      curl ftp://cucurbitgenomics.org/pub/cucurbit/genome/melon/v3.5.1/CM3.5.1_GO_anno.txt.gz > CM3.5.1_GO_anno.txt.gz
      echo Done downlaoding, now unzipping...
      gunzip CM3.5.1_GO_anno.txt
      echo Done.
  else
      echo C. melo proteins already present.
  fi
  cd ../../
fi
