#!/bin/bash -l
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=8G
#SBATCH --nodes=1
#SBATCH --time=4:00:00
#SBATCH --mail-user=araje002@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH -o /bigdata/littlab/arajewski/FULTranscriptomes/logs/STARAlignment-%A_%a.out
set -e

# This is meant to be run as an array job with IDs between 0 and 7 to specify the experiment to map
Expt=(SlycIH SlycSRA Spimp Nobt TAIRRNA TAIRChIP Nobt Nobt_All_SE)
OUTDIR=STAR/${Expt[$SLURM_ARRAY_TASK_ID]}
mkdir -p $OUTDIR

# Set STAR variables based on array ID
case "$SLURM_ARRAY_TASK_ID" in
  "0")
    SampleList=( AC1DPA1 AC1DPA2 AC1DPA3 AC3DPA1 AC3DPA2 AC3DPA3 AC15DPA1 AC15DPA2 AC15DPA3 ACbreaker1 ACbreaker2 ACbreaker3 ACrr1 ACrr2 ACrr3 )
    INDEXDIR=SlycDNA
    INDIR=SlycRNA
    ISPE=1
    ;;
  "1")
    SampleList=( SRR943813 SRR943814 SRR943815 SRR943816 SRR943817 SRR943818 SRR943825 SRR943826 SRR943827 SRR943828 SRR943829 SRR943830 )
    INDEXDIR=SlycDNA
    INDIR=ExternalData/RNAseq
    ISPE=0
    ;;
  "2")
    SampleList=( PIMP1DPA1 PIMP1DPA2 PIMP1DPA3 PIMP3DPA1 PIMP3DPA2 PIMP3DPA3 PIMP15DPA1 PIMP15DPA2 PIMP15DPA3 PIMPbreaker1 PIMPbreaker2 PIMPbreaker3 PIMPrr1 PIMPrr2 PIMPrr3 )
    INDEXDIR=SlycDNA
    INDIR=SpimpRNA
    ISPE=1
    ;;
  "3")
    SampleList=( NobtPre1 NobtPre2 NobtPre3 Nobt3DPA1 Nobt3DPA2 Nobt3DPA3 Nobt6DPA1 Nobt6DPA2 Nobt6DPA3 )
    INDEXDIR=NobtDNA
    INDIR=NobtRNA
    ISPE=1
    ;;
  "4")
    SampleList=( ERR2809794 ERR2809795 ERR2809796 ERR2809797 ERR2809798 ERR2809799 ERR2809800 ERR2809801 ERR2809802 ERR2809803 ERR2809804 ERR2809805 ERR2809806 ERR2809807 ERR2809808 ERR2809809 ERR2809810 ERR2809811 ERR2809812 ERR2809813 ERR2809814 ERR2809815 ERR2809816 ERR2809817 )
    INDEXDIR=ExternalData/TAIR10
    INDIR=ExternalData/RNAseq
    ISPE=0
    ;;
  "5")
    SampleList=( SRR6412402 SRR6412403 SRR6412404 SRR3288009 SRR3288010 SRR3288011 SRR3288012 )
    INDEXDIR=ExternalData/TAIR10
    INDIR=ExternalData/ChIPseq
    ISPE=0
    ;;
  "6")
    SampleList=( Nobtbr1 Nobtbr2 Nobtbr3 Nobtbr1.1 Nobtbr2.1 Nobtbr3.1 )
    INDEXDIR=NobtDNA
    INDIR=NobtRNA
    ISPE=0
    ;;
  "7")
    SampleList=( NobtPre1 NobtPre2 NobtPre3 Nobt3DPA1 Nobt3DPA2 Nobt3DPA3 Nobt6DPA1 Nobt6DPA2 Nobt6DPA3 )
    INDEXDIR=NobtDNA
    INDIR=NobtRNA
    ISPE=0
    ;;
esac

# General STAR commands
module load STAR/2.5.3a
for i in ${SampleList[@]}; do
  if [ ! -e $OUTDIR/${i}.Aligned.sortedByCoord.out.bam ]; then
    echo Mapping $i
    #check if the mapping should be for paired end data
    if [ "$ISPE" -ge "1" ]; then
      STAR \
        --runThreadN $SLURM_CPUS_PER_TASK \
  	    --genomeDir $INDEXDIR/ \
  	    --outFileNamePrefix $OUTDIR/$i. \
  	    --outSAMtype BAM SortedByCoordinate \
  	    --readFilesIn $INDIR/${i}_val_*.fq.gz \
  	    --readFilesCommand zcat
    else
	    STAR \
        --runThreadN $SLURM_CPUS_PER_TASK \
  	    --genomeDir $INDEXDIR/ \
  	    --outFileNamePrefix $OUTDIR/$i. \
  	    --outSAMtype BAM SortedByCoordinate \
  	    --readFilesIn $INDIR/${i}_1_trimmed.fq.gz \
  	    --readFilesCommand zcat
	  fi
	  # remove useless logfiles and splice junction information
	  rm $OUTDIR/*.out $OUTDIR/*.tab
	  echo Done.
  else
    echo $i already mapped.
  fi
done
