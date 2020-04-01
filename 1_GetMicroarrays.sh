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

#For the tomato gene expression microarray, I cannot find the gene list with sequences publically available for direct download. Instad I have copies the data table from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL10570 to a new text file called Agilent022270.txt

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
    grep -Fwf ExternalData/Microarray/Agilent022270.exonerate.dups.txt ExternalData/Microarray/Agilent022270.exonerate.out | less #Print multimappers with targets for personal entertainment
    tail -n +3 ExternalData/Microarray/Agilent022270.exonerate.out | cut -f2 -d " " |sort |uniq -dc | less #Show how many matches each multimapper has (between 2 and 80)
    #Remove multimappers from the output
    grep -vFwf ExternalData/Microarray/Agilent022270.exonerate.dups.txt ExternalData/Microarray/Agilent022270.exonerate.out > ExternalData/Microarray/Agilent022270.exonerate.nodups.txt
    #Create a final list of just the single mapping probes
    tail -n +3 ExternalData/Microarray/Agilent022270.exonerate.nodups.txt |cut -f2 -d " " | head -n -1 >ExternalData/Microarray/Agilent022270.nodups.txt
fi

#now for the Nimblegen Array, pray for me.
#similar ot the Agilent array, I simply downloaded that datatabel with the probe names and sequences from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL15968 and named it GPL15968-24228.txt
if [ ! -e ExternalData/Microarray/GPL15968-24228.fa ]; then
    tail -n +8 ExternalData/Microarray/GPL15968-24228.txt | cut -f1,2 | awk '{if ($2) print ">"$1"\n"$2}' > ExternalData/Microarray/GPL15968-24228.fa
fi
#Map these targets to the genome tolerating no mismatches
#consider outputting with --showtargetgff TRUE as a help for mapping these loci back to the genome maybe?
if [ ! -e ExternalData/Microarray/GPL15968-24228.nodups.txt ]; then
    module load exonerate
    exonerate \
        --showalignment FALSE \
        --percent 100 \
        -S no \
        --bestn 10 \
	--showtargetgff FALSE \
        ExternalData/Microarray/GPL15968-24228.fa \
        SlycDNA/S_lycopersicum_chromosomes.4.00.fa > ExternalData/Microarray/GPL15968-24228.exonerate.out
    #get a list of duplicates from the exonerate output
    tail -n +3 ExternalData/Microarray/GPL15968-24228.exonerate.out | cut -f2 -d " " |sort |uniq -d > ExternalData/Microarray/GPL15968-24228.exonerate.dups.txt
    grep -Fwf ExternalData/Microarray/GPL15968-24228.exonerate.dups.txt ExternalData/Microarray/GPL15968-24228.exonerate.out | less #Print multimappers with targets for personal entertainment
    tail -n +3 ExternalData/Microarray/GPL15968-24228.exonerate.out | cut -f2 -d " " |sort |uniq -dc | less#Show how many matches each multimapper has (between 2 and 80)
    #Remove multimappers from the output
    grep -vFwf ExternalData/Microarray/GPL15968-24228.exonerate.dups.txt ExternalData/Microarray/GPL15968-24228.exonerate.out > ExternalData/Microarray/GPL15968-24228.exonerate.nodups.txt
    #Create a final list of just the single mapping
    tail -n +3 ExternalData/Microarray/GPL15968-24228.exonerate..nodups.txt |cut -f2 -d " " | head -n -1 >ExternalData/Microarray/GPL15968-24228.nodups.txt
fi



#The arabidopsis tiling array is mapped to the TAIR9 genome. I want to update this to be consistent, but I am having trouble finding the design...stand by


#Download the CEL files for the Arabidopsis Affy Expt
if [ ! -e ExternalData/Microarray/GSE79553/GSM2098198_GR0b-v4.CEL.gz ]; then
    echo Downlaoding Arabidopsis microarray GSE79553...
    mkdir -p ExternalData/Microarray/GSE79553
    curl ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE79nnn/GSE79553/suppl/GSE79553_RAW.tar > ExternalData/Microarray/GSE79553/GSE79553_RAW.tar
    tar -C ExternalData/Microarray/GSE79553/ -xvf ExternalData/Microarray/GSE79553/GSE79553_RAW.tar
    echo Done.
fi
