#!/bin/bash -l
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=7G
#SBATCH --nodes=1
#SBATCH --time=2:00:00
#SBATCH --mail-user=araje002@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH -o /bigdata/littlab/arajewski/FULTranscriptomes/logs/STARIndex-%A.out
set -e

TAIRDIR=/rhome/arajewski/bigdata/FULTranscriptomes/ExternalData/TAIR10

#Change chromosome names for TAIR.fa to match TAIR.gff3
sed 's/>\([1-5]\).*/>Chr\1/' $TAIRDIR/TAIR10.fa | sed 's/>mito.*/>ChrM/' | sed 's/>chloro.*/>ChrC/' > $TAIRDIR/TAIR10.fa2
mv $TAIRDIR/TAIR10.fa2 $TAIRDIR/TAIR10.fa

module load STAR/2.5.3a
#Make index Files for STAR Alignments
STAR \
    --runThreadN $SLURM_CPUS_PER_TASK \
    --runMode genomeGenerate \
    --genomeDir $TAIRDIR/ \
    --genomeFastaFiles $TAIRDIR/TAIR10.fa \
    --sjdbGTFfile $TAIRDIR/TAIR10.gff3 \
    --sjdbOverhang 100 \
    --sjdbGTFtagExonParentTranscript Parent 

#Add if statement to prevent remaking. 

#Add Nobt and Slyc indcies
