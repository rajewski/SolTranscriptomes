#!/bin/bash -l
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --nodes=1
#SBATCH --mail-user=rajewski23@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -o ./logs/SRADownload-%A_%a.out
set -eu

# Load SRA Toolkit
module load sratoolkit/3.0.0

# Get vars based on array ID
accession=$(awk "NR==$SLURM_ARRAY_TASK_ID" SRA_IDs.tsv | cut -f1)
stem=$(awk "NR==$SLURM_ARRAY_TASK_ID" SRA_IDs.tsv | cut -f2)

# Get FASTQs
echo "Dumping FASTQ(s) for $accession with name $stem"
fastq-dump \
	--defline-seq '@$sn[_$rn]/$ri' \
	--defline-qual '+$sn[_$rn]/$ri' \
	--split-3 \
	--gzip \
	-B \
	-O ./SRA/ \
	$accession

# Rename the file
# Thanks to: https://stackoverflow.com/a/40029320/13954432
find . -type f | sed -n "s/$accession\(_[1|2]\).fastq.gz/& Test\1.fastq.gz/p" | xargs -n 2 mv
