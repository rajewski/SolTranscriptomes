#!/bin/bash -l
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --nodes=1
#SBATCH --time=04-00:00:00
#SBATCH --mail-user=araje002@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH -o /bigdata/littlab/arajewski/FULTranscriptomes/logs/NoninteractiveClustering-%A_%a.out

if [ $SLURM_ARRAY_TASK_ID == 1 ]; then
    Rscript 4_NoninteractiveClustering.R TAIR
fi

if [ $SLURM_ARRAY_TASK_ID == 2 ]; then
    Rscript 4_NoninteractiveClustering.R SlycIH
fi

if [ $SLURM_ARRAY_TASK_ID == 3 ]; then
    Rscript 4_NoninteractiveClustering.R SlycSRA
fi

if [ $SLURM_ARRAY_TASK_ID == 4 ]; then
    Rscript 4_NoninteractiveClustering.R Nobt
fi

if [ $SLURM_ARRAY_TASK_ID == 5 ]; then
    Rscript 4_NoninteractiveClustering.R SlycIH3Stage
fi

if [ $SLURM_ARRAY_TASK_ID == 6 ]; then
    Rscript 4_NoninteractiveClustering.R TAIR3stage
fi
