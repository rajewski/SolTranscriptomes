#!/bin/bash -l
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --nodes=1
#SBATCH --time=1:00:00
#SBATCH --mail-user=araje002@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH -o ./logs/GetExternalData-%A_%a.out
set -e

# Define a list of necessary external data sources
# This  should be run as an array job with indices specified below
# Index 0 for Arabidopsis data
# Index 1-3 for microarray data
# Index 4-51 are RNAseq 
# Index  52-58 are for ChIP data
# Index 59 is for the melon data
AllExternal=(TAIR10 GSE41560 GSE49125 GSE79553 SRR943813 SRR943814 SRR943815 SRR943816 SRR943817 SRR943818 SRR943825 SRR943826 SRR943827 SRR943828 SRR943829 SRR943830 ERR2809794 ERR2809795 ERR2809796 ERR2809797 ERR2809798 ERR2809799 ERR2809800 ERR2809801 ERR2809802 ERR2809803 ERR2809804 ERR2809805 ERR2809806 ERR2809807 ERR2809808 ERR2809809 ERR2809810 ERR2809811 ERR2809812 ERR2809813 ERR2809814 ERR2809815 ERR2809816 ERR2809817 SRR3199616 SRR3199617 SRR3199618 SRR3199656 SRR3199657 SRR3199622 SRR3199623 SRR3199624 SRR3199628 SRR3199629 SRR3199630 SRR3199635 SRR6412402 SRR6412403 SRR6412404 SRR3288009 SRR3288010 SRR3288011 SRR3288012 Melon )

####### Download Arabidopsis Data
if [ "$SLURM_ARRAY_TASK_ID" == 0 ]; then
  cd  ExternalData/${AllExternal[$SLURM_ARRAY_TASK_ID]}
  # Get genome fasta
  if [ ! -e ${AllExternal[$SLURM_ARRAY_TASK_ID]}.fa ]; then
      echo Downloading ${AllExternal[$SLURM_ARRAY_TASK_ID]} genome...
      curl https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas > ${AllExternal[$SLURM_ARRAY_TASK_ID]}.fa
      echo Done.
  else
      echo ${AllExternal[$SLURM_ARRAY_TASK_ID]} genome already present.
  fi
  # Get annotation gff
  if [ ! -e ${AllExternal[$SLURM_ARRAY_TASK_ID]}.gff3 ]; then
      echo Downloading ${AllExternal[$SLURM_ARRAY_TASK_ID]} genome annotation...
      curl https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff > ${AllExternal[$SLURM_ARRAY_TASK_ID]}.gff3
      echo Done.
  else
      echo ${AllExternal[$SLURM_ARRAY_TASK_ID]} genome annotation already present.
  fi
  # Get protein fasta
  if [ ! -e ${AllExternal[$SLURM_ARRAY_TASK_ID]}.proteins.fa ]; then
      echo Downloading ${AllExternal[$SLURM_ARRAY_TASK_ID]} proteins...
      curl https://www.arabidopsis.org/download_files/Proteins/TAIR10_protein_lists/TAIR10_pep_20110103_representative_gene_model > ${AllExternal[$SLURM_ARRAY_TASK_ID]}.proteins.fa
      echo Done.
  else
      echo ${AllExternal[$SLURM_ARRAY_TASK_ID]} proteins already present.
  fi
  cd ../../
fi

####### Tomato Agilent Array
if [ "$SLURM_ARRAY_TASK_ID" == 1 ]; then
  mkdir -p ExternalData/Microarray/${AllExternal[$SLURM_ARRAY_TASK_ID]}
  cd ExternalData/Microarray/${AllExternal[$SLURM_ARRAY_TASK_ID]}
  # remove the header, select probe name and sequence columns and output to a fasta
  if [ ! -e Agilent022270.fa ]; then
      tail -n +18 Agilent022270.txt | cut -f4,16 | awk '{if ($2) print ">"$1"\n"$2}' > Agilent022270.fa
  fi
  # map those probes to the transcriptome with exonerate tolerating no mismatches
  if [ ! -e Agilent022270.final.txt ]; then
      module load exonerate
      exonerate \
      --showalignment FALSE \
  	  --percent 100 \
  	  -S no \
  	  --bestn 10 \
  	  Agilent022270.fa \
  	  ../../../SlycDNA/Slyc.transcripts.fa > Agilent022270.exonerate.out
      #get a list of duplicates from the exonerate output
      tail -n +3 Agilent022270.exonerate.out | cut -f2 -d " " |sort |uniq -d > Agilent022270.exonerate.dups.txt
      grep -Fwf Agilent022270.exonerate.dups.txt Agilent022270.exonerate.out | less #Print multimappers with targets for personal entertainment
      tail -n +3 Agilent022270.exonerate.out | cut -f2 -d " " |sort |uniq -dc | less #Show how many matches each multimapper has (between 2 and 80)
      #Remove multimappers from the output
      grep -vFwf Agilent022270.exonerate.dups.txt Agilent022270.exonerate.out > Agilent022270.exonerate.nodups.txt
      #Create a final list of just the single mapping probes and one with targets
      tail -n +3 Agilent022270.exonerate.nodups.txt | cut -f2 -d " " | head -n -1 > Agilent022270.nodups.txt
      tail -n +3 Agilent022270.exonerate.nodups.txt | cut -f2,6 -d " "  | head -n -1 |  sed 's/mRNA://' > Agilent022270.final.txt
  fi
  # Download the raw data files for the Tomato Agilent Array
  if [ ! -e GSM1019121_252227010044_A04.txt.gz ]; then
      wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE41nnn/GSE41560/suppl/GSE41560_RAW.tar
      tar -xvf GSE41560_RAW.tar
  fi
  cd ../../../
fi

####### Tomato Nimblegen Array
if [ "$SLURM_ARRAY_TASK_ID" == 2 ]; then
  mkdir -p ExternalData/Microarray/${AllExternal[$SLURM_ARRAY_TASK_ID]}
  cd ExternalData/Microarray/${AllExternal[$SLURM_ARRAY_TASK_ID]}
  # I downloaded that data table with the probe names and sequences from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL15968 and named it GPL15968-24228.txt
  # clean this and convert to fasta
  if [ ! -e GPL15968-24228.fa ]; then
      tail -n +8 GPL15968-24228.txt | cut -f1,2 | awk '{if ($2) print ">"$1"\n"$2}' > GPL15968-24228.fa
      cd ../
  fi
  # map these probes to the genome tolerating no mismatches
  if [ ! -e GPL15968-24228.nodups.txt ]; then
      module load exonerate
      exonerate \
          --showalignment FALSE \
          --percent 100 \
          -S no \
          --bestn 10 \
        	--showtargetgff FALSE \
          GPL15968-24228.fa \
          ../../../SlycDNA/S_lycopersicum_chromosomes.4.00.fa > GPL15968-24228.exonerate.out
      #get a list of duplicates from the exonerate output
      tail -n +3 GPL15968-24228.exonerate.out | cut -f2 -d " " |sort |uniq -d > GPL15968-24228.exonerate.dups.txt
      grep -Fwf GPL15968-24228.exonerate.dups.txt GPL15968-24228.exonerate.out | less #Print multimappers with targets for personal entertainment
      tail -n +3 GPL15968-24228.exonerate.out | cut -f2 -d " " |sort |uniq -dc | less #Show how many matches each multimapper has (between 2 and 80)
      #Remove multimappers from the output
      grep -vFwf GPL15968-24228.exonerate.dups.txt GPL15968-24228.exonerate.out > GPL15968-24228.exonerate.nodups.txt
      #Create a final list of just the single mapping and another with data suitable for import into Ringo
      tail -n +3 GPL15968-24228.exonerate.nodups.txt |cut -f2 -d " " | head -n -1 > GPL15968-24228.nodups.txt
      tail -n +3 GPL15968-24228.exonerate.nodups.txt |head -n -1 | cut -f2,6,7,8,9 -d " " > GPL15968-24228.final.txt
      # get a list of the random probes
      zgrep "RANDOM" GPL15968_110707_RDKK310_Slyc_ChIP.ndf.gz | cut -f13,5 > GPL15968-24228.randomprobes.txt
      cd ../
  fi
  # download the raw data files for the Tomato Nimblegen array
  if [ ! -e GSM1194061_71972605_ratio_peaks_mapToFeatures_All_Peaks.txt.gz ]; then
      wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE49nnn/GSE49125/suppl/GSE49125_RAW.tar
      tar -xvf GSE49125_RAW.tar
  fi
  cd ../../../
fi

####### Tomato Nimblegen Array
if [ "$SLURM_ARRAY_TASK_ID" == 3 ]; then
  mkdir -p ExternalData/Microarray/${AllExternal[$SLURM_ARRAY_TASK_ID]}
  cd ExternalData/Microarray/${AllExternal[$SLURM_ARRAY_TASK_ID]}
  # Having done this once, I need to say that this data is bad. Consider ignoring it totally
  # Affy arrays are too complex for me to remap them so I will skip that steps
  # Download the CEL files for the Arabidopsis Affy Expt
  if [ ! -e GSM2098198_GR0b-v4.CEL.gz ]; then
      echo Downlaoding Arabidopsis microarray ${AllExternal[$SLURM_ARRAY_TASK_ID]}...
      curl ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE79nnn/GSE79553/suppl/GSE79553_RAW.tar > ${AllExternal[$SLURM_ARRAY_TASK_ID]}_RAW.tar
      tar -xvf GSE79553_RAW.tar
      echo Done.
  fi
  # Load in the normalized data from NCBI just in case
  if [ ! -e GSE79553_Normalized_data_with_all_controls.txt.gz ]; then
      echo Getting non-log2 normalized data from NCBI
      curl https://ftp.ncbi.nlm.nih.gov/geo/series/GSE79nnn/GSE79553/suppl/GSE79553%5FNormalized%5Fdata%5Fwith%5Fall%5Fcontrols%2Etxt%2Egz > GSE79553_Normalized_data_with_all_controls.txt.gz
      gunzip GSE79553_Normalized_data_with_all_controls.txt.gz
  fi
  # Download the CDF for that Affy array
  if [ ! -e GPL14926_agronomics1_TAIR9_gene.CDF.gz ]; then
      echo Downloading Affy CDF file
      curl https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL14nnn/GPL14926/suppl/GPL14926%5Fagronomics1%5FTAIR9%5Fgene%2ECDF%2Egz > GPL14926_agronomics1_TAIR9_gene.CDF.gz
      gunzip GPL14926_agronomics1_TAIR9_gene.CDF.gz
      echo Done
  fi
  # try another CDF for the Affy array
  if [ ! -e AGRONOMICS1_At_TAIRG.cdf ]; then
      wget http://mbni.org/customcdf/18.0.0/tairg.download/AGRONOMICS1_At_TAIRG_18.0.0.zip
      unzip AGRONOMICS1_At_TAIRG_18.0.0.zip
  fi
  cd ../../../
fi

####### Download SRA RNAseq Data
if [ "$SLURM_ARRAY_TASK_ID" -ge 4 ] && [ "$SLURM_ARRAY_TASK_ID" -le 51 ]; then
  mkdir -p ExternalData/RNAseq
  cd ExternalData/RNAseq
  # Prefetch Data
  if [ ! -d ${AllExternal[$SLURM_ARRAY_TASK_ID]} ]; then
      echo Downloading ${AllExternal[$SLURM_ARRAY_TASK_ID]} data from SRA...
      module load sratoolkit/2.10.0
      prefetch ${AllExternal[$SLURM_ARRAY_TASK_ID]}
      echo Done.
  else
      echo ${AllExternal[$SLURM_ARRAY_TASK_ID]} data already present.
  fi
  # Dump Fastq files for data
  if [ ! -e ${AllExternal[$SLURM_ARRAY_TASK_ID]}_1.fastq.gz ]; then
      echo Dumping fastq data for ${AllExternal[$SLURM_ARRAY_TASK_ID]}...
      module load sratoolkit/2.10.0
      fastq-dump --defline-seq '@$sn[_$rn]/$ri' --defline-qual '+$sn[_$rn]/$ri' --split-files --gzip -B ${AllExternal[$SLURM_ARRAY_TASK_ID]}
      echo Done.
  else
      echo Fastq data for ${AllExternal[$SLURM_ARRAY_TASK_ID]} already present.
  fi
  #Trim Reads
  if [ ! -e ${AllExternal[$SLURM_ARRAY_TASK_ID]}_1_trimmed.fq.gz ]; then
      echo Running Trim Galore on ${AllExternal[$SLURM_ARRAY_TASK_ID]}...
      module unload perl
      module swap miniconda2 miniconda3
      module swap python/2.7.5 python/3.6.0
      module load fastqc
      conda activate cutadaptenv
      export PATH=/bigdata/littlab/arajewski/Datura/software/TrimGalore-0.6.5:$PATH
      trim_galore \
	  --no_report_file \
	  -j $SLURM_CPUS_PER_TASK \
          ${AllExternal[$SLURM_ARRAY_TASK_ID]}_1.fastq.gz
      echo Done.
  else
      echo ${AllExternal[$SLURM_ARRAY_TASK_ID]} RNA seq already trimmed.
  fi
  cd ../../
fi

####### Download SRA ChIP-seq Data
if [ "$SLURM_ARRAY_TASK_ID" -ge 52 ] && [ "$SLURM_ARRAY_TASK_ID" -le 58 ]; then
  mkdir -p ExternalData/ChIPseq
  cd ExternalData/ChIPseq
  # Prefetch Data
  if [ ! -d ${AllExternal[$SLURM_ARRAY_TASK_ID]} ]; then
      echo Downloading ${AllExternal[$SLURM_ARRAY_TASK_ID]} data from SRA...
      module load sratoolkit/2.10.0
      prefetch ${AllExternal[$SLURM_ARRAY_TASK_ID]}
      echo Done.
  else
      echo ${AllExternal[$SLURM_ARRAY_TASK_ID]} data already present.
  fi
  #Dump Fastq files 
  if [ ! -e ${AllExternal[$SLURM_ARRAY_TASK_ID]}_1.fastq.gz ]; then
      echo Dumping fastq data for ${AllExternal[$SLURM_ARRAY_TASK_ID]}...
      module load sratoolkit/2.10.0
      fastq-dump --defline-seq '@$sn[_$rn]/$ri' --defline-qual '+$sn[_$rn]/$ri' --split-files --gzip -B ${AllExternal[$SLURM_ARRAY_TASK_ID]}
      echo Done.
  else
      echo Fastq data for ${AllExternal[$SLURM_ARRAY_TASK_ID]} already present.
  fi
  #Trim Reads                                                                                                                                                                      
  if [ ! -e ${AllExternal[$SLURM_ARRAY_TASK_ID]}_1_trimmed.fq.gz ]; then
      echo Running Trim Galore on ${AllExternal[$SLURM_ARRAY_TASK_ID]}...
      module load trim_galore/0.4.2
      trim_galore \
          --no_report_file \
          ${AllExternal[$SLURM_ARRAY_TASK_ID]}_1.fastq.gz
      echo Done.
  else
      echo ${AllExternal[$SLURM_ARRAY_TASK_ID]} RNA seq already trimmed.
  fi
  cd ../../
fi

# Melon Genome and RNA seq data
if [ "$SLURM_ARRAY_TASK_ID" == 59 ]; then
    mkdir -p ExternalData/C_melo
    cd  ExternalData/C_melo
  # Get genome fasta
  if [ ! -e CM3.5.1_genome.fa ]; then
      echo Downloading C. melo genome fasta genome...
      curl ftp://cucurbitgenomics.org/pub/cucurbit/genome/melon/v3.5.1/CM3.5.1_genome.fa.gz > CM3.5.1_genome.fa.gz
      echo Done downloading, now unzipping...
      gunzip CM3.5.1_genome.fa.gz
      echo Done.
  else
      echo C. melo genome already present.
  fi
  # Get annotation gff                                                                                                                                                           
  if [ ! -e CM3.5.1_gene.gff ]; then
      echo Downloading C. melo genome annotation...
      curl ftp://cucurbitgenomics.org/pub/cucurbit/genome/melon/v3.5.1/CM3.5.1_gene.gff.gz > CM3.5.1_gene.gff.gz
      echo Done downloading, now unzipping...
      gunzip CM3.5.1_gene.gff.gz
      echo Done.
  else
      echo C. melo genome annotation already present.
  fi
  # Get protein fasta
  if [ ! -e CM3.5.1_protein.fa ]; then
      echo Downloading C. melo proteins...
      curl ftp://cucurbitgenomics.org/pub/cucurbit/genome/melon/v3.5.1/CM3.5.1_protein.fa.gz > CM3.5.1_protein.fa.gz
      echo Done downlaoding, now unzipping...
      gunzip CM3.5.1_protein.fa.gz
      echo Done.
  else
      echo C. melo proteins already present.
  fi
  # Get GO Terms
  if [ ! -e CM3.5.1_GO_anno.txt ]; then
      echo Downloading C. melo GO Terms...
      curl ftp://cucurbitgenomics.org/pub/cucurbit/genome/melon/v3.5.1/CM3.5.1_GO_anno.txt.gz > CM3.5.1_GO_anno.txt.gz
      echo Done downlaoding, now unzipping...
      gunzip CM3.5.1_GO_anno.txt
      echo Done.
  else
      echo C. melo proteins already present.
  fi
  cd ../../
fi
