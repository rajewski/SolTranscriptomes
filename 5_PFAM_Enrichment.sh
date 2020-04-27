#!/bin/bash -l
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=60
#SBATCH --mem-per-cpu=7G
#SBATCH --nodes=1
#SBATCH -p short
#SBATCH --time=02:00:00
#SBATCH --mail-user=araje002@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH -o /bigdata/littlab/arajewski/FULTranscriptomes/logs/PFAM-%A.out
set -e

# This analysis is adapted from the Supplemental Text S1 of https://doi.org/10.1104/pp.108.132985

# Make a variable to store all the protein fasta files
PEPs=( ExternalData/TAIR10/TAIR10.proteins.fa NobtDNA/NIOBT_r1.0.proteins.fa SlycDNA/Slyc.proteins.fa)
#Run an interpro scan on the protein fasta to get a list of PFAM domains for each protein in the genome.
####Add check to avoid repeating this step####
module load interproscan
for PEP in ${!PEPs[@]}
do
    if [ ! -e DEGAnalysis/Pfam/$(basename ${PEPs[PEP]}).tsv ]; then
	echo Running IPRScan on $(basename ${PEPs[PEP]})...
	interproscan.sh \
	    -i ${PEPs[PEP]} \
	    -f TSV \
	    -appl Pfam,TIGRFAM,PRINTS,ProSiteProfiles \
	    --goterms \
	    -dra \
	    -d DEGAnalysis/Pfam \
	    -cpu $SLURM_CPUS_PER_TASK
   fi
done

# Clean up the Interpro output files
cd DEGAnalysis/Pfam

# Get a list of all proteins' names
grep ">" ../../SlycDNA/Slyc.proteins.fa | sed 's/\S*\(Solyc\S*\)\s.*/\1/' | sort > Slyc.protein.names.txt
grep ">" ../../NobtDNA/NIOBT_r1.0.proteins.fa | sed 's/>\(\S*\)\s\S*\s\S*/\1/g' | sort > Nobt.protein.names.txt
grep ">" ../../ExternalData/TAIR10/TAIR10.proteins.fa | sed 's/>\(\S*\)\s.*/\1/' | sort > TAIR10.protein.names.txt

# Make an associative array to aid with renaming of the outputs
declare -A PFAMs
PFAMs=([TAIR10]=TAIR10.proteins.fa.tsv [Nobt]=NIOBT_r1.0.proteins.fa.tsv [Slyc]=Slyc.proteins.fa.tsv )
for TAB in ${!PFAMs[@]}
do
    cut -f1 ${PFAMs[$TAB]} | sort | uniq > $TAB.pfamhits.tsv # Names of all genes with a pfam hit
    #cut -f1,5 ${PFAMs[$TAB]} |sort | uniq > $TAB.gene2pfam.tsv # Association of all genes and their pfam domains
    cut -f1,12 ${PFAMs[$TAB]} |sort |uniq |grep "IPR" > $TAB.gene2ipr.tsv # Association of all genes and their IPR domains
    #cut -f5,6 ${PFAMs[$TAB]} > $TAB.pfam2desc.tsv # Descriptions of pfam domains
    cut -f12,13 ${PFAMs[$TAB]} > $TAB.ipr2desc.tsv # Descriptions of IPR domains
    cut -f1,14 ${PFAMs[$TAB]} | sort |uniq | grep "GO" > $TAB.genes2go.tsv # Get the genes to go mapping
    comm -1 -3 $TAB.pfamhits.tsv $TAB.protein.names.txt > $TAB.nopfam.tsv # Names of all genes without a pfam hit
done
cd ../../


#I want to create a bed file of the transcription start sites for each gene in tomato (arabidopsis to follow) so that I can track the distance between the TSS and the nearest FUL ChIP binding site
cd ChIPAnalysis
# Make a bed file of just the TSS
awk 'BEGIN{FS=OFS="\t"}($7=="+" && $3=="mRNA"){print $1,$2,"TSS",$4,$4+2,$6,$7,$8,$9}' ../SlycDNA/ITAG4.0_gene_models.gff > Slyc.TSS.gff
awk 'BEGIN{FS=OFS="\t"}($7=="-" && $3=="mRNA"){print $1,$2,"TSS",$5-2,$5,$6,$7,$8,$9}' ../SlycDNA/ITAG4.0_gene_models.gff >> Slyc.TSS.gff
sort -k1,1 -k4,4n Slyc.TSS.gff > Slyc.TSS.sort.gff
module load bedops/2.4.24
convert2bed --input=GFF  < Slyc.TSS.sort.gff > Slyc.TSS.sort.bed
# Get distance from every TSS to nearest upstream FUL binding site
module load bedtools/2.28.0
bedtools closest \
    -a Slyc.TSS.sort.bed \
    -b ChIP-chip/chers.bed \
    -D a


cd ../../
