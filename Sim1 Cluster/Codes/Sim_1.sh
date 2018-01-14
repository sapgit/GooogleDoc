#!/bin/bash
#
#SBATCH --job-name=Sim_1_negbin
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

Rscript --no-save --no-restore --verbose Sim_1.R ${RHO} 'negbin' ${PHI} ${SAMPLE_SIZE} 20 $SLURM_ARRAY_TASK_ID > Sim_1_${RHO}_negbin_${PHI}_${SAMPLE_SIZE}_$SGE_TASK_ID.Rout
