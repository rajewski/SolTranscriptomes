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

RNAseq=( SRR943813 SRR943814 SRR943815 SRR943816 SRR943817 SRR943818 SRR943825 SRR943826 SRR943827 SRR943828 SRR943829 SRR943830 ERR2809794 ERR2809795 ERR2809796 ERR2809797 ERR2809798 ERR2809799 ERR2809806 ERR2809807 ERR2809808 ERR2809815 ERR2809816 ERR2809817 )

cd /rhome/arajewski/bigdata/FULTranscriptomes/ExternalData/RNAseq

#Prefetch Data
if [ ! -d ${RNAseq[$SLURM_ARRAY_TASK_ID]} ]; then
    echo Downloading ${RNAseq[$SLURM_ARRAY_TASK_ID]} data from SRA...
    module load sratoolkit/2.10.0
    prefetch ${RNAseq[$SLURM_ARRAY_TASK_ID]}
    echo Done.
else
    echo ${RNAseq[$SLURM_ARRAY_TASK_ID]} data already present.
fi

#Dump Fastq files for data
if [ ! -e ${RNAseq[$SLURM_ARRAY_TASK_ID]}_1.fastq.gz ]; then
    echo Dumping fastq data for ${RNAseq[$SLURM_ARRAY_TASK_ID]}...
    module load sratoolkit/2.10.0
    fastq-dump --defline-seq '@$sn[_$rn]/$ri' --defline-qual '+$sn[_$rn]/$ri' --split-files --gzip -B ${RNAseq[$SLURM_ARRAY_TASK_ID]}
    echo Done.
else
    echo Fastq data for ${RNAseq[$SLURM_ARRAY_TASK_ID]} already present.
fi

#Trim Reads
if [ ! -e ${RNAseq[$SLURM_ARRAY_TASK_ID]}_1_trimmed.fq.gz ]; then
    echo Running Trim Galore on ${RNAseq[$SLURM_ARRAY_TASK_ID]}...
    module load trim_galore/0.4.2
    trim_galore \
        --no_report_file \
        ${RNAseq[$SLURM_ARRAY_TASK_ID]}_1.fastq.gz
    echo Done.
else
    echo ${RNAseq[$SLURM_ARRAY_TASK_ID]} RNA seq already trimmed.
fi
