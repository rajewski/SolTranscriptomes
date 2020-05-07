#!/bin/bash -l
#SBATCH --ntasks=16
#SBATCH --nodes=1
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --mail-user=araje002@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --out ../logs/ReciprocalHits-%A.out
set -eu

#The goal of this script is to take the DEG list from Nobt and Slyc produced by DESeq2 and 
#look for orthologous genes between them.

#Make Diamond search databased for Nobt and Slyc
module load diamond/0.9.24
if [ ! -e Nobt.dmnd ]; then
  echo $(date): Making Nicotiana obtusifolia Diamond database.
  diamond makedb \
    --in ~/bigdata/Nobtusifolia/Genome_Files/NIOBT_r1.0.proteins.fa \
    --db Nobt \
    -p $SLURM_NTASKS
  echo $(date): Done.
else
  echo $(date): N. obtusifolia Diamond database already present.
fi

if [ ! -e Slyc.dmnd ]; then
  echo $(date): Making Solanum lycopersicum Diamond database.
  diamond makedb \
    --in ~/bigdata/Slycopersicum/data/ITAG3.2_proteins.fasta \
    --db Slyc \
    -p $SLURM_NTASKS
  echo $(date): Done.
else
  echo $(date): S. lycopersicum Diamond database already present
fi

#Find the best hit interspecific for each of the top 100 genes in the Slyc/Nobt Up/Down datasets
#then search that hit against the original species and see if the input matches the output.
#if yes, then call it orthologous.

#########Nobt Up
#search Nobt DEGs against Slyc DB
if [ ! -e NobtUptoSlyc.dmndo ]; then
  echo $(date): Searching Up-regulated Nobt genes against Slyc DB with DIAMOND
  diamond blastp \
    --db Slyc.dmnd \
    -p $SLURM_NTASKS \
    --sensitive \
    --outfmt 6 qseqid sseqid evalue full_sseq \
    --max-target-seqs 1 \
    --query ~/bigdata/Nobtusifolia/RNA-seq/Results_Ballgown/GoAnalysis/NobtUp.fasta \
    -o NobtUptoSlyc.dmndo
  echo $(date): Done.
  cut -f1,2 NobtUptoSlyc.dmndo > NobtUptoSlycHits.tsv
  awk '{print ">"$2"\n"$4}' NobtUptoSlyc.dmndo > NobtUptoSlycHits.fasta
else
  echo $(date): Nobt to Slyc search already complete.
fi

#search Slyc hits back against Nobt DB
if [ ! -e NobtUptoNobt.dmndo ]; then
  echo $(date): Searching Up-regulated Nobt gene hits back against Nobt DB with DIAMOND
  diamond blastp \
    --db Nobt.dmnd \
    -p $SLURM_NTASKS \
    --sensitive \
    --outfmt 6 qseqid sseqid evalue \
    --max-target-seqs 1 \
    --query NobtUptoSlycHits.fasta \
    -o NobtUptoNobt.dmndo
  echo $(date): Done.
  cut -f1,2 NobtUptoNobt.dmndo > NobtUptoNobtHits.tsv
else
  echo $(date): Nobt to Nobt search already complete.
fi

#########Nobt Down
#search Nobt DEGs against Slyc DB
if [ ! -e NobtDowntoSlyc.dmndo ]; then
  echo $(date): Searching Down-regulated Nobt genes against Slyc DB with DIAMOND
  diamond blastp \
    --db Slyc.dmnd \
    -p $SLURM_NTASKS \
    --sensitive \
    --outfmt 6 qseqid sseqid evalue full_sseq \
    --max-target-seqs 1 \
    --query ~/bigdata/Nobtusifolia/RNA-seq/Results_Ballgown/GoAnalysis/NobtDown.fasta \
    -o NobtDowntoSlyc.dmndo
  echo $(date): Done.
  cut -f1,2 NobtDowntoSlyc.dmndo > NobtDowntoSlycHits.tsv
  awk '{print ">"$2"\n"$4}' NobtDowntoSlyc.dmndo > NobtDowntoSlycHits.fasta
else
  echo $(date): Nobt to Slyc search already complete.
fi

#search Slyc hits back against Nobt DB
if [ ! -e NobtDowntoNobt.dmndo ]; then
  echo $(date): Searching Down-regulated Nobt gene hits back against Nobt DB with DIAMOND
  diamond blastp \
    --db Nobt.dmnd \
    -p $SLURM_NTASKS \
    --sensitive \
    --outfmt 6 qseqid sseqid evalue \
    --max-target-seqs 1 \
    --query NobtDowntoSlycHits.fasta \
    -o NobtDowntoNobt.dmndo
  echo $(date): Done.
  cut -f1,2 NobtDowntoNobt.dmndo > NobtDowntoNobtHits.tsv
else
  echo $(date): Nobt to Nobt search already complete.
fi


#####Tomato#####

#########Nobt Up
#search Slyc DEGs against Nobt DB
if [ ! -e SlycUptoNobt.dmndo ]; then
  echo $(date): Searching Up-regulated Slyc genes against Nobt DB with DIAMOND
  diamond blastp \
    --db Nobt.dmnd \
    -p $SLURM_NTASKS \
    --sensitive \
    --outfmt 6 qseqid sseqid evalue full_sseq \
    --max-target-seqs 1 \
    --query ~/bigdata/Nobtusifolia/RNA-seq/Results_Ballgown/GoAnalysis/SlycUp.fasta \
    -o SlycUptoNobt.dmndo
  echo $(date): Done.
  cut -f1,2 SlycUptoNobt.dmndo > SlycUptoNobtHits.tsv
  awk '{print ">"$2"\n"$4}' SlycUptoNobt.dmndo > SlycUptoNobtHits.fasta
else
  echo $(date): Slyc to Nobt search already complete.
fi

#search Slyc hits back against Nobt DB
if [ ! -e SlycUptoSlyc.dmndo ]; then
  echo $(date): Searching Up-regulated Slyc gene hits back against Slyc DB with DIAMOND
  diamond blastp \
    --db Slyc.dmnd \
    -p $SLURM_NTASKS \
    --sensitive \
    --outfmt 6 qseqid sseqid evalue \
    --max-target-seqs 1 \
    --query SlycUptoNobtHits.fasta \
    -o SlycUptoSlyc.dmndo
  echo $(date): Done.
  cut -f1,2 SlycUptoSlyc.dmndo > SlycUptoSlycHits.tsv
else
  echo $(date): Nobt to Nobt search already complete.
fi

#########Slyc Down
#search Slyc DEGs against Nobt DB
if [ ! -e SlycDowntoNobt.dmndo ]; then
  echo $(date): Searching Down-regulated Slyc genes against Nobt DB with DIAMOND
  diamond blastp \
    --db Nobt.dmnd \
    -p $SLURM_NTASKS \
    --sensitive \
    --outfmt 6 qseqid sseqid evalue full_sseq \
    --max-target-seqs 1 \
    --query ~/bigdata/Nobtusifolia/RNA-seq/Results_Ballgown/GoAnalysis/SlycDown.fasta \
    -o SlycDowntoNobt.dmndo
  echo $(date): Done.
  cut -f1,2 SlycDowntoNobt.dmndo > SlycDowntoNobtHits.tsv
  awk '{print ">"$2"\n"$4}' SlycDowntoNobt.dmndo > SlycDowntoNobtHits.fasta
else
  echo $(date): Slyc to Nobt search already complete.
fi

#search Nobt hits back against Slyc DB
if [ ! -e SlycDowntoSlyc.dmndo ]; then
  echo $(date): Searching Down-regulated Slyc gene hits back against Slyc DB with DIAMOND
  diamond blastp \
    --db Slyc.dmnd \
    -p $SLURM_NTASKS \
    --sensitive \
    --outfmt 6 qseqid sseqid evalue \
    --max-target-seqs 1 \
    --query SlycDowntoNobtHits.fasta \
    -o SlycDowntoSlyc.dmndo
  echo $(date): Done.
  cut -f1,2 SlycDowntoSlyc.dmndo > SlycDowntoSlycHits.tsv
else
  echo $(date): Slyc to Slyc search already complete.
fi


