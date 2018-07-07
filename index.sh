#!/bin/bash -l
#SBATCH --ntasks=5
#SBATCH --nodes=1
#SBATCH --mem=40G
#SBATCH --time=2:00:00
#SBATCH --mail-user=araje002@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH -p batch

# $((SLURM_MEM_PER_NODE/1000))'G'
# $SLURM_NTASKS

module load kallisto/0.44.0
echo `date` "Loading Kallisto v 0.440"

echo `date` "Building index of transcripts"
kallisto index -i data/transcripts.idx data/ITAG3.2_cDNA.fasta
echo `date` "Done."
