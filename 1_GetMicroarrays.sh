#!/bin/bash -l
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=7G
#SBATCH --nodes=1
#SBATCH --time=2:00:00
#SBATCH --mail-user=araje002@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH -o /bigdata/littlab/arajewski/FULTranscriptomes/logs/GetMicroarrays.out
set -e

#For the two tomato microarray based experiments I am going to remap the probe sequences to the genome or transcriptome (depending) and either confirm that they bind where they bind or find out what gene they bind in. The ChIP-chip array from nimblegen already has gene names associated but they are based on ITAG2, which is a several years old by now.

#First for the Tomato Agilent Array
#remove the header, select probe name and sequence columns and output to a fasta
if [ ! -e ExternalData/Microarray/Agilent022270.fa ]; then
    tail -n +18 ExternalData/Microarray/Agilent022270.txt | cut -f4,16 | awk '{if ($2) print ">"$1"\n"$2}' > ExternalData/Microarray/Agilent022270.fa
fi

#map those probes to the transcriptome with exonerate tolerating no mismatches
if [ ! -e ExternalData/Microarray/Agilent022270.nodups.txt ]; then
    module load exonerate
    exonerate \
	--showalignment FALSE \
	--percent 100 \
	-S no \
	--bestn 10 \
	ExternalData/Microarray/Agilent022270.fa \
	SlycDNA/Slyc.transcripts.fa > ExternalData/Microarray/Agilent022270.exonerate.out
    #get a list of duplicates from the exonerate output
    tail -n +3 ExternalData/Microarray/Agilent022270.exonerate.out | cut -f2 -d " " |sort |uniq -d > ExternalData/Microarray/Agilent022270.exonerate.dups.txt
    grep -Fwf ExternalData/Microarray/Agilent022270.exonerate.dups.txt ExternalData/Microarray/Agilent022270.exonerate.out #Print multimappers with targets for personal entertainment
    tail -n +3 ExternalData/Microarray/Agilent022270.exonerate.out | cut -f2 -d " " |sort |uniq -dc #Show how many matches each multimapper has (between 2 and 80)
    #Remove multimappers from the output
    grep -vFwf ExternalData/Microarray/Agilent022270.exonerate.dups.txt ExternalData/Microarray/Agilent022270.exonerate.out > ExternalData/Microarray/Agilent022270.exonerate.nodups.txt
    #Create a final list of just the single mapping probes
    tail -n +3 ExternalData/Microarray/Agilent022270.exonerate.nodups.txt |cut -f2 -d " " | head -n -1 >ExternalData/Microarray/Agilent022270.nodups.txt
fi

#now for the Nimblegen Array, pray for me.
