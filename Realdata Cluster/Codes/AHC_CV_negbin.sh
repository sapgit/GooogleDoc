#!/bin/bash
#
#SBATCH --job-name=AHC_CV_negbin
#SBATCH --ntasks=1
#SBATCH --partition=medium
#
# Time format = HH:MM:SS, DD-HH:MM:SS
#
#SBATCH --time=50:00:00
#
# Minimum memory required per allocated  CPU  in  MegaBytes. 
#
#SBATCH --mem-per-cpu=2000
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=himel@uab.edu
#
#SBATCH --array=1-10

module load R/3.3.1-foss-2016b-no-X

Rscript --no-save --no-restore --verbose main.R 10 $SLURM_ARRAY_TASK_ID > AHC_CV_negbin_$SGE_TASK_ID.Rout 