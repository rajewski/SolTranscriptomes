#!/bin/bash -l
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --nodes=1
#SBATCH --time=01-00:00:00
#SBATCH --mail-user=araje002@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH -o /bigdata/littlab/arajewski/FULTranscriptomes/logs/MACS-%A.out
set -e

module load MACS2/2.2.6
INDIR=ChIPAnalysis/ChIPseq/STAR
OUTDIR=ChIPAnalysis/ChIPseq/MACS
mkdir -p $OUTDIR
# My plan is to call peaks for each sample separately and then combine them all later
# based on https://hbctraining.github.io/Intro-to-ChIPseq/lessons/07_handling-replicates-idr.html

# Analyze PRJNA316152 Expt
macs2 callpeak \
    -t $INDIR/SRR3288009.Aligned.sortedByCoord.out.bam  \
    -c $INDIR/SRR3288011.Aligned.sortedByCoord.out.bam  \
    -n PRJNA316152_rep1 \
    --outdir $OUTDIR \
    -g 1.2e8 \
    -B \
    --call-summits \
    -p 1e-3

macs2 callpeak \
    -t $INDIR/SRR3288010.Aligned.sortedByCoord.out.bam \
    -c $INDIR/SRR3288012.Aligned.sortedByCoord.out.bam \
    -n PRJNA316152_rep2 \
    --outdir $OUTDIR \
    -g 1.2e8 \
    -B \
    --call-summits \
    -p 1e-3

# Analyze PRJNA427320 Expt
macs2 callpeak \
    -t $INDIR/SRR6412402.Aligned.sortedByCoord.out.bam  \
    -c $INDIR/SRR6412404.Aligned.sortedByCoord.out.bam  \
    -n PRJNA427320_rep1 \
    --outdir $OUTDIR \
    -g 1.2e8 \
    -B \
    --call-summits \
    -p 1e-3

macs2 callpeak \
    -t $INDIR/SRR6412403.Aligned.sortedByCoord.out.bam \
    -c $INDIR/SRR6412404.Aligned.sortedByCoord.out.bam  \
    -n PRJNA427320_rep2 \
    --outdir $OUTDIR \
    -g 1.2e8 \
    -B \
    --call-summits \
    -p 1e-3


# Get the intersection of all peaks (this is very stringent, but I'll take it
module load bedtools
bedtools intersect \
    -a $OUTDIR/PRJNA316152_rep1_peaks.narrowPeak \
    -b $OUTDIR/PRJNA316152_rep2_peaks.narrowPeak $OUTDIR/PRJNA427320_rep1_peaks.narrowPeak $OUTDIR/PRJNA427320_rep2_peaks.narrowPeak  > $OUTDIR/TAIR10_FUL_overlaps.narrowPeak
