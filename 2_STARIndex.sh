#!/bin/bash -l
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=20G
#SBATCH --nodes=1
#SBATCH --time=2:00:00
#SBATCH --mail-user=araje002@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH -o /bigdata/littlab/arajewski/FULTranscriptomes/logs/STARIndex-%A.out
set -e

TAIRDIR=/rhome/arajewski/bigdata/FULTranscriptomes/ExternalData/TAIR10
NobtDIR=/rhome/arajewski/shared/Nobtusifolia/Genome_Files
SlycDIR=/rhome/arajewski/bigdata/Datura/Alkaloids/ExternalData/Slyc

module load STAR/2.5.3a
#Make index Files for Arabidopsis
if [ ! -e $TAIRDIR/SAindex ]; then
    echo Making STAR index for Arabidopsis...
    #Change chromosome names for TAIR.fa to match TAIR.gff3
    sed 's/>\([1-5]\).*/>Chr\1/' $TAIRDIR/TAIR10.fa | sed 's/>mito.*/>ChrM/' | sed 's/>chloro.*/>ChrC/' > $TAIRDIR/TAIR10.fa2
    mv $TAIRDIR/TAIR10.fa2 $TAIRDIR/TAIR10.fa
    STAR \
	--runThreadN $SLURM_CPUS_PER_TASK \
	--runMode genomeGenerate \
	--genomeDir $TAIRDIR/ \
	--genomeFastaFiles $TAIRDIR/TAIR10.fa \
	--sjdbGTFfile $TAIRDIR/TAIR10.gff3 \
	--sjdbOverhang 100 \
	--sjdbGTFtagExonParentTranscript Parent
    echo Done.
else
    echo STAR index for Arabidopsis already present.
fi

#Make index files for Nobt
if [ ! -e $NobtDIR/SAindex ]; then
    echo Making STAR index for N. obtusifolia...
    STAR \
        --runThreadN $SLURM_CPUS_PER_TASK \
        --runMode genomeGenerate \
        --genomeDir $NobtDIR/ \
        --genomeFastaFiles $NobtDIR/NIOBT_r1.0.fasta \
        --sjdbGTFfile $NobtDIR/NIOBT_r1.0.update.gff \
        --sjdbOverhang 100 \
        --sjdbGTFtagExonParentTranscript Parent \
	--genomeChrBinNbits 16 \
	--limitSjdbInsertNsj 150000
    #Add limitSjdbInsertNsj because of error (https://groups.google.com/forum/#!msg/rna-star/ddhJDgvZfNA/ULUGGYb0BgAJ)
    #Change genomeChrBinNbits because the genome is highly fragmented
    echo Done.
else
    echo STAR index for N. obtusifolia already present.
fi

#Make index files for Slyc
if [ ! -e $SlycDIR/SAindex ]; then
    echo Making STAR index for S. lycopersicum...
    STAR \
        --runThreadN $SLURM_CPUS_PER_TASK \
        --runMode genomeGenerate \
        --genomeDir $SlycDIR/ \
        --genomeFastaFiles $SlycDIR/S_lycopersicum_chromosomes.4.00.fa \
        --sjdbGTFfile $SlycDIR/ITAG4.0_gene_models.gff \
        --sjdbOverhang 100 \
        --sjdbGTFtagExonParentTranscript Parent
    echo Done.
else
    echo STAR index for S. lycopersicum already present.
fi
