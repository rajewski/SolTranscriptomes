#!/bin/bash -l
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem-per-cpu=4G
#SBATCH --nodes=1
#SBATCH --time=1-00:00:00
#SBATCH --mail-user=araje002@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH -o /bigdata/littlab/arajewski/FULTranscriptomes/logs/Orthofinder-%A.out
set -e

# Make symlinks for all the protein fasta files
ln -s ../ExternalData/TAIR10/TAIR10.proteins.fa Arabidopsis.fa
ln -s ../SlycDNA/Slyc.proteins.fa Solanum.fa
ln -s ../NobtDNA/NIOBT_r1.0.proteins.fa Nicotiana.fa

module load orthofinder/2.3.7
module load mafft/7.427
module load IQ-TREE/1.6.12

orthofinder \
    -t $SLURM_CPUS_PER_TASK \
    -M msa \
    -T iqtree \
    -f ./
