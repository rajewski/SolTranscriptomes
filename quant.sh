#!/bin/bash -l
#SBATCH --ntasks=10
#SBATCH --nodes=1
#SBATCH --mem=40G
#SBATCH --time=3:00:00
#SBATCH --mail-user=araje002@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH -p batch

# $((SLURM_MEM_PER_NODE/1000))'G'
# $SLURM_NTASKS

echo `date` "loading Kallisto 0.44.0"
module load kallisto/0.44.0

#determine the rep to work on for this job (must be an array job)
RepList=~/bigdata/Slycopersicum/replist.txt
Sample=$(awk "NR==$SLURM_ARRAY_TASK_ID" $RepList)

echo `date` "Running Kallisto on $Sample"

kallisto quant -i data/transcripts.idx -b 100 -o "Results_"$Sample --genomebam --gtf data/ITAG3.2_gene_models.gff --chromosomes data/chromosomes.txt data/"AC"$Sample"LEFT" data/"AC"$Sample"RIGHT"
