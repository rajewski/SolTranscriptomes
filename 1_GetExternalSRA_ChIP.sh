#!/bin/bash -l
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=7G
#SBATCH --nodes=1
#SBATCH --time=2:00:00
#SBATCH --mail-user=araje002@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH -o /bigdata/littlab/arajewski/FULTranscriptomes/logs/GetExternalData-%A.out
set -e

mkdir -p /rhome/arajewski/bigdata/FULTranscriptomes/ExternalData/ChIPseq
cd /rhome/arajewski/bigdata/FULTranscriptomes/ExternalData/ChIPseq

# For the Chip-seq in Arabidopsis, there are two experiments, but I don't want all of the data from each
ChIPseq=(SRR6412402 SRR6412403 SRR6412404 SRR3288009 SRR3288010 SRR3288011 SRR3288012)

#Prefetch Data
if [ ! -d ${ChIPseq[$SLURM_ARRAY_TASK_ID]} ]; then
    echo Downloading ${ChIPseq[$SLURM_ARRAY_TASK_ID]} data from SRA...
    module load sratoolkit/2.10.0
    prefetch ${ChIPseq[$SLURM_ARRAY_TASK_ID]}
    echo Done.
else
    echo ${ChIPseq[$SLURM_ARRAY_TASK_ID]} data already present.
fi

#Dump Fastq files 
if [ ! -e ${ChIPseq[$SLURM_ARRAY_TASK_ID]}_1.fastq.gz ]; then
    echo Dumping fastq data for ${ChIPseq[$SLURM_ARRAY_TASK_ID]}...
    module load sratoolkit/2.10.0
    fastq-dump --defline-seq '@$sn[_$rn]/$ri' --defline-qual '+$sn[_$rn]/$ri' --split-files --gzip -B ${ChIPseq[$SLURM_ARRAY_TASK_ID]}
    echo Done.
else
    echo Fastq data for ${ChIPseq[$SLURM_ARRAY_TASK_ID]} already present.
fi

#Trim Reads                                                                                                                                                                      
if [ ! -e ${ChIPseq[$SLURM_ARRAY_TASK_ID]}_1_trimmed.fq.gz ]; then
    echo Running Trim Galore on ${ChIPseq[$SLURM_ARRAY_TASK_ID]}...
    module load trim_galore/0.4.2
    trim_galore \
        --no_report_file \
        ${ChIPseq[$SLURM_ARRAY_TASK_ID]}_1.fastq.gz
    echo Done.
else
    echo ${ChIPseq[$SLURM_ARRAY_TASK_ID]} RNA seq already trimmed.
fi
