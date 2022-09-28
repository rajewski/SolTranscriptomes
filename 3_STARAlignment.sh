#!/bin/bash -l
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=8G
#SBATCH --nodes=1
##SBATCH --time=4:00:00
#SBATCH --mail-user=rajewski23@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -o ./logs/STARAlignment-%A_%a.out
set -e

# This is meant to be run as an array job with IDs between 0 and 8 to specify the experiment to map
Expt=( SlycIH Spimp Nobt TAIRRNA Nobt Nobt_All_SE Slyc_SE Spimp_SE Melon QC )
OUTDIR=STAR/${Expt[$SLURM_ARRAY_TASK_ID]}
mkdir -p $OUTDIR

# Set STAR variables based on array ID
# I'm adding a value for ISPE. 0 means SE, 1 means true PE, and 2 means forced to be SE
case "$SLURM_ARRAY_TASK_ID" in
  "0")
    SampleList=( AC1DPA1 AC1DPA2 AC1DPA3 AC3DPA1 AC3DPA2 AC3DPA3 AC15DPA1 AC15DPA2 AC15DPA3 ACbreaker1 ACbreaker2 ACbreaker3 ACrr1 ACrr2 ACrr3 )
    INDEXDIR=ExternalData/Slyc
    INDIR=SRA
    ISPE=1
    REFBED=/bigdata/littlab/arajewski/FULTranscriptomes/ExternalData/Slyc/ITAG4.0_gene_models.bed
    REFGFF=/bigdata/littlab/arajewski/FULTranscriptomes/ExternalData/Slyc/ITAG4.0_gene_models.gff
    ;;
  "1")
    SampleList=( SP1DPA1 SP1DPA2 SP1DPA3 SP3DPA1 SP3DPA2 SP3DPA3 SP15DPA1 SP15DPA2 SP15DPA3 SPbreaker1 SPbreaker2 SPbreaker3 SPrr1 SPrr2 SPrr3 )
    INDEXDIR=ExternalData/Slyc
    INDIR=SRA
    ISPE=1
    REFBED=/bigdata/littlab/arajewski/FULTranscriptomes/ExternalData/Slyc/ITAG4.0_gene_models.bed
    REFGFF=/bigdata/littlab/arajewski/FULTranscriptomes/ExternalData/Slyc/ITAG4.0_gene_models.gff
    ;;
  "2")
    SampleList=( NobtPre1 NobtPre2 NobtPre3 Nobt3DPA1 Nobt3DPA2 Nobt3DPA3 Nobt6DPA1 Nobt6DPA2 Nobt6DPA3 )
    INDEXDIR=NobtDNA
    INDIR=SRA
    ISPE=1
    REFBED=/bigdata/littlab/shared/Nobtusifolia/Genome_Files/NIOBT_r1.0.update.bed
    REFGFF=/bigdata/littlab/shared/Nobtusifolia/Genome_Files/NIOBT_r1.0.update.gff
    ;;
  "3")
    SampleList=( ERR2809794 ERR2809795 ERR2809796 ERR2809797 ERR2809798 ERR2809799 ERR2809800 ERR2809801 ERR2809802 ERR2809803 ERR2809804 ERR2809805 ERR2809806 ERR2809807 ERR2809808 ERR2809809 ERR2809810 ERR2809811 ERR2809812 ERR2809813 ERR2809814 ERR2809815 ERR2809816 ERR2809817 )
    INDEXDIR=ExternalData/TAIR10
    INDIR=ExternalData/RNAseq
    ISPE=0
    REFBED=/bigdata/littlab/arajewski/FULTranscriptomes/ExternalData/TAIR10/TAIR10.genes.bed
    REFGFF=/bigdata/littlab/arajewski/FULTranscriptomes/ExternalData/TAIR10/TAIR10.gff
    ;;
  "4")
    SampleList=( Nobtbr1 Nobtbr2 Nobtbr3 Nobtbr1.1 Nobtbr2.1 Nobtbr3.1 )
    INDEXDIR=NobtDNA
    INDIR=SRA
    ISPE=0
    REFBED=/bigdata/littlab/shared/Nobtusifolia/Genome_Files/NIOBT_r1.0.update.bed
    REFGFF=/bigdata/littlab/shared/Nobtusifolia/Genome_Files/NIOBT_r1.0.update.gff
    ;;
  "5")
    SampleList=( NobtPre1 NobtPre2 NobtPre3 Nobt3DPA1 Nobt3DPA2 Nobt3DPA3 Nobt6DPA1 Nobt6DPA2 Nobt6DPA3 )
    INDEXDIR=NobtDNA
    INDIR=SRA
    ISPE=2
    REFBED=/bigdata/littlab/shared/Nobtusifolia/Genome_Files/NIOBT_r1.0.update.bed
    REFBED=/bigdata/littlab/shared/Nobtusifolia/Genome_Files/NIOBT_r1.0.update.gff
    ;;
  "6")
    SampleList=( AC1DPA1 AC1DPA2 AC1DPA3 AC3DPA1 AC3DPA2 AC3DPA3 AC15DPA1 AC15DPA2 AC15DPA3 ACbreaker1 ACbreaker2 ACbreaker3 ACrr1 ACrr2 ACrr3 )
    INDEXDIR=ExternalData/Slyc
    INDIR=SRA
    ISPE=2
    REFBED=/bigdata/littlab/arajewski/FULTranscriptomes/ExternalData/Slyc/ITAG4.0_gene_models.bed
    REFGFF=/bigdata/littlab/arajewski/FULTranscriptomes/ExternalData/Slyc/ITAG4.0_gene_models.gff
    ;;
  "7")
    SampleList=( SP1DPA1 SP1DPA2 SP1DPA3 SP3DPA1 SP3DPA2 SP3DPA3 SP15DPA1 SP15DPA2 SP15DPA3 SPbreaker1 SPbreaker2 SPbreaker3 SPrr1 SPrr2 SPrr3 )
    INDEXDIR=ExternalData/Slyc
    INDIR=SRA
    ISPE=2
    REFBED=/bigdata/littlab/arajewski/FULTranscriptomes/ExternalData/Slyc/ITAG4.0_gene_models.bed
    REFGFF=/bigdata/littlab/arajewski/FULTranscriptomes/ExternalData/Slyc/ITAG4.0_gene_models.gff
    ;;
  "8")
    SampleList=( SRR3199616 SRR3199617 SRR3199618 SRR3199656 SRR3199657 SRR3199622 SRR3199623 SRR3199624 SRR3199628 SRR3199629 SRR3199630 SRR3199635 )
    INDEXDIR=ExternalData/C_melo
    INDIR=ExternalData/RNAseq
    ISPE=0
    REFBED=/bigdata/littlab/arajewski/FULTranscriptomes/ExternalData/C_melo/CM3.5.1_gene.bed
    REFGFF=/bigdata/littlab/arajewski/FULTranscriptomes/ExternalData/C_melo/CM3.5.1_gene.gff
    ;;
  "9")
    module load singularity/3.9.3
    echo "Performing Final QC on Mapping"
    singularity exec SIFs/multiQC_1.13.sif multiqc \
	--force \
	--outdir ./ \
	--ignore "*_SE/*" \
	--ignore "*.SJ.out.tab" \
	STAR
    echo "Done"
    exit 0
esac

# General STAR commands
module load star/2.7.10a
for i in ${SampleList[@]}; do
  if [ ! -e $OUTDIR/${i}.Aligned.sortedByCoord.out.bam ]; then
    echo Mapping $i
    #check if the mapping should be for paired end data
    if [ "$ISPE" -eq "1" ]; then
	STAR \
	    --runThreadN $SLURM_CPUS_PER_TASK \
  	    --genomeDir $INDEXDIR/ \
  	    --outFileNamePrefix $OUTDIR/$i. \
  	    --outSAMtype BAM SortedByCoordinate \
  	    --readFilesIn $INDIR/${i}_*_val_*.fq.gz \
  	    --readFilesCommand zcat
    elif [ "$ISPE" -eq "0" ]; then
	STAR \
	    --runThreadN $SLURM_CPUS_PER_TASK \
  	    --genomeDir $INDEXDIR/ \
  	    --outFileNamePrefix $OUTDIR/$i. \
  	    --outSAMtype BAM SortedByCoordinate \
  	    --readFilesIn $INDIR/${i}_1_trimmed.fq.gz \
  	    --readFilesCommand zcat
    elif [ "$ISPE" -eq "2" ]; then
	# Just count the first read
	STAR \
	   --runThreadN $SLURM_CPUS_PER_TASK \
            --genomeDir $INDEXDIR/ \
            --outFileNamePrefix $OUTDIR/$i. \
            --outSAMtype BAM SortedByCoordinate \
            --readFilesIn $INDIR/${i}_1_val_1.fq.gz \
            --readFilesCommand zcat
    fi
    # remove useless logfiles and splice junction information
    # rm $OUTDIR/*.out $OUTDIR/*.tab # I need these logs for reviewers
    echo Done.
  else
      echo "$i already mapped."
      # Don't QC the stuff that was forced to SE
      if [ "$ISPE" -eq 2 ]; then
	echo "No QC for forced single-end experiments"
	exit 0
      fi
      module load singularity
      # Qualimap
      if [ ! -e "$OUTDIR/QC/${i}/BamQC_${i}.pdf" ]; then
      echo "Performing QC of BAM file for $i"
      singularity exec \
	-B /bigdata/littlab/shared/Nobtusifolia/Genome_Files/:/bigdata/littlab/shared/Nobtusifolia/Genome_Files/ \
	SIFs/qualimap_2.2.1.sif qualimap bamqc \
	-bam $OUTDIR/${i}.Aligned.sortedByCoord.out.bam \
	-gff $REFGFF \
	-outdir $OUTDIR/QC/${i}/ \
	-outfile BamQC_${i}.pdf \
	--java-mem-size=8G
      fi
      if [ ! -e "$OUTDIR/QC/${i}/RSeQC_${i}.out" ]; then
      echo "Indexing BAM file"
      module load samtools
      samtools index $OUTDIR/${i}.Aligned.sortedByCoord.out.bam
      echo "Performing RSeQC for $i"
      # RSeQC
      singularity exec \
        -B /bigdata/littlab/shared/Nobtusifolia/Genome_Files/:/bigdata/littlab/shared/Nobtusifolia/Genome_Files/ \
	SIFs/RSeQC_4.0.0.sif read_distribution.py  \
        -i $OUTDIR/${i}.Aligned.sortedByCoord.out.bam \
        -r $REFBED > $OUTDIR/QC/${i}/RSeQC_${i}.out
      fi
  fi
done


