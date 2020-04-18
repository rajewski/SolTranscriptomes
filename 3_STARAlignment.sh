#!/bin/bash -l
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=8G
#SBATCH --nodes=1
#SBATCH --time=2:00:00
#SBATCH --mail-user=araje002@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH -o /bigdata/littlab/arajewski/FULTranscriptomes/logs/STARAlignment-%A.out
set -e

TAIRDIR=/rhome/arajewski/bigdata/FULTranscriptomes/ExternalData/TAIR10
SRADIR=/rhome/arajewski/bigdata/FULTranscriptomes/ExternalData/RNAseq
NobtDIR=/rhome/arajewski/shared/Nobtusifolia/Genome_Files
SlycDIR=/rhome/arajewski/bigdata/Datura/Alkaloids/ExternalData/Slyc
SpimpDIR=SpimpRNA

module load STAR/2.5.3a
#Map Arabidopsis RNAseq Reads
SRAList=( ERR2809815 ERR2809798 ERR2809807 ERR2809799 ERR2809816 ERR2809808 ERR2809817 ERR2809794 ERR2809796 ERR2809806 ERR2809797 ERR2809795 )

for SRA in ${SRAList[@]}; do
    if [ ! -s DEGAnalysis/STAR/TAIR10/${SRA}.Aligned.sortedByCoord.out.bam ]; then
	echo Mapping Arabidopsis reads for $SRA...
	STAR \
	    --runThreadN $SLURM_CPUS_PER_TASK \
	    --genomeDir $TAIRDIR/ \
	    --outFileNamePrefix DEGAnalysis/STAR/TAIR10/$SRA. \
	    --outSAMtype BAM SortedByCoordinate \
	    --readFilesIn $SRADIR/${SRA}_1_trimmed.fq.gz \
	    --readFilesCommand zcat
	echo Done.
    else
	echo Reads for $SRA already mapped.
    fi
done

#Map Arabidopsis ChIP reads
TAIRCHIPDIR=/rhome/arajewski/bigdata/FULTranscriptomes/ExternalData/ChIPseq
ChIPList=(SRR6412402 SRR6412403 SRR6412404 SRR3288009 SRR3288010 SRR3288011 SRR3288012)
mkdir -p /rhome/arajewski/bigdata/FULTranscriptomes/ChIPAnalysis/ChIPseq/STAR

for SRA in ${ChIPList[@]}; do
    if [ ! -s ChIPAnalysis/ChIPseq/STAR/${SRA}.Aligned.sortedByCoord.out.bam ]; then
        echo Mapping Arabidopsis reads for $SRA...
        STAR \
            --runThreadN $SLURM_CPUS_PER_TASK \
            --genomeDir $TAIRDIR/ \
            --outFileNamePrefix ChIPAnalysis/ChIPseq/STAR/${SRA}. \
            --outSAMtype BAM SortedByCoordinate \
            --readFilesIn $TAIRCHIPDIR/${SRA}_1_trimmed.fq.gz \
            --readFilesCommand zcat
        echo Done.
    else
        echo Reads for $SRA already mapped.
    fi
done


#Map Tomato SRA reads comment out to run in parallel with others
SRALISTTom=( SRR943813 SRR943814 SRR943815 SRR943816 SRR943817 SRR943818 SRR943825 SRR943826 SRR943827 SRR943828 SRR943829 SRR943830 )

for SRA in ${SRALISTTom[@]}; do
    if [ ! -s DEGAnalysis/STAR/Slyc/${SRA}.Aligned.sortedByCoord.out.bam ]; then
        echo Mapping Tomato reads for $SRA...
        STAR \
            --runThreadN $SLURM_CPUS_PER_TASK \
            --genomeDir $SlycDIR/ \
            --outFileNamePrefix DEGAnalysis/STAR/Slyc/$SRA. \
            --outSAMtype BAM SortedByCoordinate \
            --readFilesIn $SRADIR/${SRA}_1_trimmed.fq.gz \
	    --readFilesCommand zcat
        echo Done.
    else
        echo Reads for $SRA already mapped.
    fi
done

#Map Tomato in-house reads
IHList=(AC1DPA1 AC1DPA2 AC1DPA3 AC3DPA1 AC3DPA2 AC3DPA3 AC15DPA1 AC15DPA2 AC15DPA3 ACbreaker1 ACbreaker2 ACbreaker3 ACrr1 ACrr2 ACrr3 )

for IH in ${IHList[@]}; do
    if [ ! -s DEGAnalysis/STAR/Slyc/${IH}.Aligned.sortedByCoord.out.bam ]; then
        echo Mapping Tomato reads for $IH...
        STAR \
            --runThreadN $SLURM_CPUS_PER_TASK \
            --genomeDir $SlycDIR/ \
            --outFileNamePrefix DEGAnalysis/STAR/Slyc/$IH. \
            --outSAMtype BAM SortedByCoordinate \
            --readFilesIn SlycRNA/${IH}_val_*.fq.gz \
	    --readFilesCommand zcat
        echo Done.
    else
        echo Reads for $IH already mapped.
    fi
done

#Map Spimp in-house reads to Slyc genome
SpimpList=(PIMP1DPA1 PIMP1DPA2 PIMP1DPA3 PIMP3DPA1 PIMP3DPA2 PIMP3DPA3 PIMP15DPA1 PIMP15DPA2 PIMP15DPA3 PIMPbreaker1 PIMPbreaker2 PIMPbreaker3 PIMPrr1 PIMPrr2 PIMPrr3 )
for IH in ${SpimpList[@]}; do
    if [ ! -s DEGAnalysis/STAR/Spimp/${IH}.Aligned.sortedByCoord.out.bam ]; then
        echo Mapping Tomato reads for $IH...
        STAR \
            --runThreadN $SLURM_CPUS_PER_TASK \
            --genomeDir $SlycDIR/ \
            --outFileNamePrefix DEGAnalysis/STAR/Spimp/$IH. \
            --outSAMtype BAM SortedByCoordinate \
            --readFilesIn SpimpRNA/${IH}_val_*.fq.gz \
	    --readFilesCommand zcat
        echo Done.
    else
        echo Reads for $IH already mapped.
    fi
done


#Map Tobacco in house reads
IHListNobt=(NobtPre1 NobtPre2 NobtPre3 Nobt3DPA1 Nobt3DPA2 Nobt3DPA3 Nobt6DPA1 Nobt6DPA2 Nobt6DPA3 )

for IH in ${IHListNobt[@]}; do
    if [ ! -s DEGAnalysis/STAR/Nobt/${IH}.Aligned.sortedByCoord.out.bam ]; then
        echo Mapping Tobacco reads for $IH...
        STAR \
            --runThreadN $SLURM_CPUS_PER_TASK \
            --genomeDir $NobtDIR/ \
            --outFileNamePrefix DEGAnalysis/STAR/Nobt/$IH. \
            --outSAMtype BAM SortedByCoordinate \
            --readFilesIn NobtRNA/${IH}_*_val_*.fq.gz \
            --readFilesCommand zcat
        echo Done.
    else
        echo Reads for $IH already mapped.
    fi
done

