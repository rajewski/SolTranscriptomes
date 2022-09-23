#!/bin/bash -l
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem-per-cpu=7G
#SBATCH --nodes=1
#SBATCH --time=1-00:00:00
#SBATCH --mail-user=araje002@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH -o /bigdata/littlab/arajewski/FULTranscriptomes/logs/Orthofinder-%A.out
set -e

# Make symlinks for all the protein fasta files
#ln -sf ../ExternalData/TAIR10/TAIR10.proteins.fa Arabidopsis.fa
ln -sf ../SlycDNA/Slyc.proteins.fa Solanum.fa
ln -sf ../NobtDNA/NIOBT_r1.0.proteins.fa Nicotiana.fa
#ln -sf ../ExternalData/C_melo/CM3.5.1_protein_longestIsoforms.fa Cucumis.fa

module load orthofinder/2.4.0
module load mafft/7.427
module load IQ-TREE/1.6.12

orthofinder \
    -t $SLURM_CPUS_PER_TASK \
    -f ./JustSol
