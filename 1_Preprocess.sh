#!/bin/bash -l
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=7G
#SBATCH --nodes=1
#SBATCH --time=04:00:00
#SBATCH --mail-user=araje002@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH -o /bigdata/littlab/arajewski/FULTranscriptomes/logs/Rcorrector-%A.out
set -e

# Do a kmer based read correction with RCorrector
# Reads are going to be corrected by species batch
BATCHES=(Nobt Slyc Spimp TAIR10)
module load jellyfish/2.3.0

if [ "${BATCHES[$SLURM_ARRAY_TASK_ID]}" == "Nobt" ]; then
  cd NobtRNA
  mkdir -p Rcorrector
  perl /bigdata/littlab/arajewski/Datura/software/rcorrector/run_rcorrector.pl \
    -1 NobtPre1_1_val_1.fq.gz,Nobt6DPA3_1_val_1.fq.gz,Nobt6DPA2_1_val_1.fq.gz,Nobt3DPA2_1_val_1.fq.gz,NobtPre3_1_val_1.fq.gz,NobtPre2_1_val_1.fq.gz,Nobt3DPA3_1_val_1.fq.gz,Nobt3DPA1_1_val_1.fq.gz,Nobt6DPA1_1_val_1.fq.gz\
    -2 NobtPre1_2_val_2.fq.gz,Nobt6DPA3_2_val_2.fq.gz,Nobt6DPA2_2_val_2.fq.gz,Nobt3DPA2_2_val_2.fq.gz,NobtPre3_2_val_2.fq.gz,NobtPre2_2_val_2.fq.gz,Nobt3DPA3_2_val_2.fq.gz,Nobt3DPA1_2_val_2.fq.gz,Nobt6DPA1_2_val_2.fq.gz\
    -t $SLURM_CPUS_PER_TASK \
    -od Rcorrector
  # remove reads where one of the pair is unfixable
  /bigdata/littlab/arajewski/Datura/software/TranscriptomeAssemblyTools/FilterUncorrectabledPEfastq.py \
  -1 Rcorrector/NobtPre1_1_val_1.fq.gz,Rcorrector/Nobt6DPA3_1_val_1.fq.gz,Rcorrector/Nobt6DPA2_1_val_1.fq.gz,Rcorrector/Nobt3DPA2_1_val_1.fq.gz,Rcorrector/NobtPre3_1_val_1.fq.gz,Rcorrector/NobtPre2_1_val_1.fq.gz,Rcorrector/Nobt3DPA3_1_val_1.fq.gz,Rcorrector/Nobt3DPA1_1_val_1.fq.gz,Rcorrector/Nobt6DPA1_1_val_1.fq.gz\
  -2 Rcorrector/NobtPre1_2_val_2.fq.gz,Rcorrector/Nobt6DPA3_2_val_2.fq.gz,Rcorrector/Nobt6DPA2_2_val_2.fq.gz,Rcorrector/Nobt3DPA2_2_val_2.fq.gz,Rcorrector/NobtPre3_2_val_2.fq.gz,Rcorrector/NobtPre2_2_val_2.fq.gz,Rcorrector/Nobt3DPA3_2_val_2.fq.gz,Rcorrector/Nobt3DPA1_2_val_2.fq.gz,Rcorrector/Nobt6DPA1_2_val_2.fq.gz
  cd ../
fi

if [ "${BATCHES[$SLURM_ARRAY_TASK_ID]}" == "Slyc" ]; then
  cd SlycRNA
  SRADIR=../ExternalData/RNAseq
  mkdir -p $SRADIR/Rcorrector
  mkdir -p Rcorrector
  mkdir 
  perl /bigdata/littlab/arajewski/Datura/software/rcorrector/run_rcorrector.pl \
    -s $SRADIR/SRR943813_1_trimmed.fq.gz,$SRADIR/SRR943814_1_trimmed.fq.gz,$SRADIR/SRR943815_1_trimmed.fq.gz,$SRADIR/SRR943816_1_trimmed.fq.gz,$SRADIR/SRR943817_1_trimmed.fq.gz,$SRADIR/SRR943818_1_trimmed.fq.gz,$SRADIR/SRR943825_1_trimmed.fq.gz,$SRADIR/SRR943826_1_trimmed.fq.gz,$SRADIR/SRR943827_1_trimmed.fq.gz,$SRADIR/SRR943828_1_trimmed.fq.gz,$SRADIR/SRR943829_1_trimmed.fq.gz,$SRADIR/SRR943830_1_trimmed.fq.gz \
    -1 ACrr2_val_1.fq.gz,ACrr3_val_1.fq.gz,ACrr1_val_1.fq.gz,AC1DPA1_val_1.fq.gz,AC1DPA3_val_1.fq.gz,AC15DPA3_val_1.fq.gz,ACbreaker1_val_1.fq.gz,ACbreaker2_val_1.fq.gz,AC15DPA2_val_1.fq.gz,AC15DPA1_val_1.fq.gz,AC3DPA3_val_1.fq.gz,ACbreaker3_val_1.fq.gz,AC1DPA2_val_1.fq.gz,AC3DPA1_val_1.fq.gz,AC3DPA2_val_1.fq.gz \
    -2 ACrr2_val_2.fq.gz,ACrr3_val_2.fq.gz,ACrr1_val_2.fq.gz,AC1DPA1_val_2.fq.gz,AC1DPA3_val_2.fq.gz,AC15DPA3_val_2.fq.gz,ACbreaker1_val_2.fq.gz,ACbreaker2_val_2.fq.gz,AC15DPA2_val_2.fq.gz,AC15DPA1_val_2.fq.gz,AC3DPA3_val_2.fq.gz,ACbreaker3_val_2.fq.gz,AC1DPA2_val_2.fq.gz,AC3DPA1_val_2.fq.gz,AC3DPA2_val_2.fq.gz \
    -t $SLURM_CPUS_PER_TASK \
    -od Rcorrector
  # Remove the unfixable paired reads
  /bigdata/littlab/arajewski/Datura/software/TranscriptomeAssemblyTools/FilterUncorrectabledPEfastq.py \
  -1 Rcorrector/ACrr2_val_1.cor.fq.gz,Rcorrector/ACrr3_val_1.cor.fq.gz,Rcorrector/ACrr1_val_1.cor.fq.gz,Rcorrector/AC1DPA1_val_1.cor.fq.gz,Rcorrector/AC1DPA3_val_1.cor.fq.gz,Rcorrector/AC15DPA3_val_1.cor.fq.gz,Rcorrector/ACbreaker1_val_1.cor.fq.gz,Rcorrector/ACbreaker2_val_1.cor.fq.gz,Rcorrector/AC15DPA2_val_1.cor.fq.gz,Rcorrector/AC15DPA1_val_1.cor.fq.gz,Rcorrector/AC3DPA3_val_1.cor.fq.gz,Rcorrector/ACbreaker3_val_1.cor.fq.gz,Rcorrector/AC1DPA2_val_1.cor.fq.gz,Rcorrector/AC3DPA1_val_1.cor.fq.gz,Rcorrector/AC3DPA2_val_1.cor.fq.gz \
  -2 Rcorrector/ACrr2_val_2.cor.fq.gz,Rcorrector/ACrr3_val_2.cor.fq.gz,Rcorrector/ACrr1_val_2.cor.fq.gz,Rcorrector/AC1DPA1_val_2.cor.fq.gz,Rcorrector/AC1DPA3_val_2.cor.fq.gz,Rcorrector/AC15DPA3_val_2.cor.fq.gz,Rcorrector/ACbreaker1_val_2.cor.fq.gz,Rcorrector/ACbreaker2_val_2.cor.fq.gz,Rcorrector/AC15DPA2_val_2.cor.fq.gz,Rcorrector/AC15DPA1_val_2.cor.fq.gz,Rcorrector/AC3DPA3_val_2.cor.fq.gz,Rcorrector/ACbreaker3_val_2.cor.fq.gz,Rcorrector/AC1DPA2_val_2.cor.fq.gz,Rcorrector/AC3DPA1_val_2.cor.fq.gz,Rcorrector/AC3DPA2_val_2.cor.fq.gz 
  # removed the unfixable reads
  /bigdata/littlab/arajewski/Datura/software/TranscriptomeAssemblyTools/FilterUncorrectabledSEfastq.py \
  -1 Rcorrector/SRR943813_1_trimmed.cor.fq.gz,Rcorrector/SRR943814_1_trimmed.cor.fq.gz,Rcorrector/SRR943815_1_trimmed.cor.fq.gz,Rcorrector/SRR943816_1_trimmed.cor.fq.gz,Rcorrector/SRR943817_1_trimmed.cor.fq.gz,Rcorrector/SRR943818_1_trimmed.cor.fq.gz,Rcorrector/SRR943825_1_trimmed.cor.fq.gz,Rcorrector/SRR943826_1_trimmed.cor.fq.gz,Rcorrector/SRR943827_1_trimmed.cor.fq.gz,Rcorrector/SRR943828_1_trimmed.cor.fq.gz,Rcorrector/SRR943829_1_trimmed.cor.fq.gz,Rcorrector/SRR943830_1_trimmed.cor.fq.gz
  mv Rcorrector/*.trimmed.cor.fq.gz $SRADIR/Rcorrector

  cd ../
fi

if [ "${BATCHES[$SLURM_ARRAY_TASK_ID]}" == "Spimp" ]; then
  cd SpimpRNA
  mkdir -p Rcorrector
  perl /bigdata/littlab/arajewski/Datura/software/rcorrector/run_rcorrector.pl \
    -1 PIMPrr2_val_1.fq.gz,PIMPrr3_val_1.fq.gz,PIMPrr1_val_1.fq.gz,PIMP1DPA1_val_1.fq.gz,PIMP1DPA3_val_1.fq.gz,PIMP15DPA3_val_1.fq.gz,PIMPbreaker1_val_1.fq.gz,PIMPbreaker2_val_1.fq.gz,PIMP15DPA2_val_1.fq.gz,PIMP15DPA1_val_1.fq.gz,PIMP3DPA3_val_1.fq.gz,PIMPbreaker3_val_1.fq.gz,PIMP1DPA2_val_1.fq.gz,PIMP3DPA1_val_1.fq.gz,PIMP3DPA2_val_1.fq.gz \
    -2 PIMPrr2_val_2.fq.gz,PIMPrr3_val_2.fq.gz,PIMPrr1_val_2.fq.gz,PIMP1DPA1_val_2.fq.gz,PIMP1DPA3_val_2.fq.gz,PIMP15DPA3_val_2.fq.gz,PIMPbreaker1_val_2.fq.gz,PIMPbreaker2_val_2.fq.gz,PIMP15DPA2_val_2.fq.gz,PIMP15DPA1_val_2.fq.gz,PIMP3DPA3_val_2.fq.gz,PIMPbreaker3_val_2.fq.gz,PIMP1DPA2_val_2.fq.gz,PIMP3DPA1_val_2.fq.gz,PIMP3DPA2_val_2.fq.gz \
    -t $SLURM_CPUS_PER_TASK \
    -od Rcorrector
  cd ../
fi

if [ "${BATCHES[$SLURM_ARRAY_TASK_ID]}" == "TAIR10" ]; then
  cd ExternalData/RNAseq
  mkdir -p Rcorrector
  perl /bigdata/littlab/arajewski/Datura/software/rcorrector/run_rcorrector.pl \
    -s ERR2809815_1_trimmed.fq.gz,ERR2809798_1_trimmed.fq.gz,ERR2809807_1_trimmed.fq.gz,ERR2809799_1_trimmed.fq.gz,ERR2809816_1_trimmed.fq.gz,ERR2809808_1_trimmed.fq.gz,ERR2809817_1_trimmed.fq.gz,ERR2809794_1_trimmed.fq.gz,ERR2809796_1_trimmed.fq.gz,ERR2809806_1_trimmed.fq.gz,ERR2809797_1_trimmed.fq.gz,ERR2809795_1_trimmed.fq.gz \
    -t $SLURM_CPUS_PER_TASK \
    -od Rcorrector
  cd ../
fi