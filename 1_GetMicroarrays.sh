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

#For the tomato gene expression microarray, I cannot find the gene list with sequences publically available for direct download. Instead I have copied the data table from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL10570 to a new text file called Agilent022270.txt

cd ExternalData/Microarray

#First for the Tomato Agilent Array
#remove the header, select probe name and sequence columns and output to a fasta
if [ ! -e GSE41560/Agilent022270.fa ]; then
    mkdir -p GSE41560
    tail -n +18 GSE41560/Agilent022270.txt | cut -f4,16 | awk '{if ($2) print ">"$1"\n"$2}' > GSE41560/Agilent022270.fa
fi

#map those probes to the transcriptome with exonerate tolerating no mismatches
if [ ! -e GSE41560/Agilent022270.nodups.txt ]; then
    module load exonerate
    exonerate \
	--showalignment FALSE \
	--percent 100 \
	-S no \
	--bestn 10 \
	GSE41560/Agilent022270.fa \
	../../SlycDNA/Slyc.transcripts.fa > GSE41560/Agilent022270.exonerate.out
    #get a list of duplicates from the exonerate output
    tail -n +3 GSE41560/Agilent022270.exonerate.out | cut -f2 -d " " |sort |uniq -d > GSE41560/Agilent022270.exonerate.dups.txt
    grep -Fwf GSE41560/Agilent022270.exonerate.dups.txt GSE41560/Agilent022270.exonerate.out | less #Print multimappers with targets for personal entertainment
    tail -n +3 GSE41560/Agilent022270.exonerate.out | cut -f2 -d " " |sort |uniq -dc | less #Show how many matches each multimapper has (between 2 and 80)
    #Remove multimappers from the output
    grep -vFwf GSE41560/Agilent022270.exonerate.dups.txt GSE41560/Agilent022270.exonerate.out > GSE41560/Agilent022270.exonerate.nodups.txt
    #Create a final list of just the single mapping probes
    tail -n +3 GSE41560/Agilent022270.exonerate.nodups.txt |cut -f2 -d " " | head -n -1 > GSE41560/Agilent022270.nodups.txt
fi

#Download the raw data files for the Tomato Agilent Array
if [ ! -e GSE41560/GSM1019121_252227010044_A04.txt.gz ]; then
    mkdir -p GSE41560
    #I cant seem to download the data directly via FTP, but I've updated it myself and the command contine from here on how to process it.
    tar -xvf GSE41560_RAW.tar
fi

#now for the Nimblegen Array, pray for me.
#similar ot the Agilent array, I simply downloaded that datatabel with the probe names and sequences from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL15968 and named it GPL15968-24228.txt
if [ ! -e GPL15968-24228.fa ]; then
    tail -n +8 GPL15968-24228.txt | cut -f1,2 | awk '{if ($2) print ">"$1"\n"$2}' > GPL15968-24228.fa
fi
#Map these targets to the genome tolerating no mismatches
#consider outputting with --showtargetgff TRUE as a help for mapping these loci back to the genome maybe?
if [ ! -e GPL15968-24228.nodups.txt ]; then
    module load exonerate
    exonerate \
        --showalignment FALSE \
        --percent 100 \
        -S no \
        --bestn 10 \
	--showtargetgff FALSE \
        GPL15968-24228.fa \
        ../../SlycDNA/S_lycopersicum_chromosomes.4.00.fa > GPL15968-24228.exonerate.out
    #get a list of duplicates from the exonerate output
    tail -n +3 GPL15968-24228.exonerate.out | cut -f2 -d " " |sort |uniq -d > GPL15968-24228.exonerate.dups.txt
    grep -Fwf GPL15968-24228.exonerate.dups.txt GPL15968-24228.exonerate.out | less #Print multimappers with targets for personal entertainment
    tail -n +3 GPL15968-24228.exonerate.out | cut -f2 -d " " |sort |uniq -dc | less#Show how many matches each multimapper has (between 2 and 80)
    #Remove multimappers from the output
    grep -vFwf GPL15968-24228.exonerate.dups.txt GPL15968-24228.exonerate.out > GPL15968-24228.exonerate.nodups.txt
    #Create a final list of just the single mapping
    tail -n +3 GPL15968-24228.exonerate..nodups.txt |cut -f2 -d " " | head -n -1 > GPL15968-24228.nodups.txt
fi

### I'm not goin to mess wit hthe Affy array data. It's too complex and also the results didnt turn out to be compelling at all.
#Download the CEL files for the Arabidopsis Affy Expt
if [ ! -e GSM2098198_GR0b-v4.CEL.gz ]; then
    echo Downlaoding Arabidopsis microarray GSE79553...
    mkdir -p GSE79553
    curl ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE79nnn/GSE79553/suppl/GSE79553_RAW.tar > GSE79553/GSE79553_RAW.tar
    tar -C GSE79553/ -xvf GSE79553/GSE79553_RAW.tar
    echo Done.
fi

#Load in the normalized data from NCBI
if [ ! -e GSE79553/GSE79553_Normalized_data_with_all_controls.txt.gz ]; then
    echo Getting non-log2 normalized data from NCBI
    curl https://ftp.ncbi.nlm.nih.gov/geo/series/GSE79nnn/GSE79553/suppl/GSE79553%5FNormalized%5Fdata%5Fwith%5Fall%5Fcontrols%2Etxt%2Egz > GSE79553/GSE79553_Normalized_data_with_all_controls.txt.gz
    cd GSE79553
    gunzip GSE79553_Normalized_data_with_all_controls.txt.gz
    cd ../
fi

#Download the CDF for that Affy array
if [ ! -e GSE79553/GPL14926_agronomics1_TAIR9_gene.CDF.gz ]; then
    echo Downloading Affy CDF file
    curl https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL14nnn/GPL14926/suppl/GPL14926%5Fagronomics1%5FTAIR9%5Fgene%2ECDF%2Egz > GSE79553/GPL14926_agronomics1_TAIR9_gene.CDF.gz
    cd GSE79553
    gunzip GPL14926_agronomics1_TAIR9_gene.CDF.gz
    cd ../
    echo Done
fi

#try anoter CDF for the Affy array
if [ ! -e GSE79553/AGRONOMICS1_At_TAIRG.cdf ]; then
    cd GSE79553
    wget http://mbni.org/customcdf/18.0.0/tairg.download/AGRONOMICS1_At_TAIRG_18.0.0.zip
    unzip AGRONOMICS1_At_TAIRG_18.0.0.zip
    cd ../
fi

