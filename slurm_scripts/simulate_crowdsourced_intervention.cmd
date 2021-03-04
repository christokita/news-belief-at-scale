#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=60gb
#SBATCH --array=0-126
#SBATCH --time=23:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ctokita@princeton.edu
#SBATCH --output=slurm_outfiles/slurm-%A_%a.out

## NOTE: We have 127 unique news stories stories as of now, hence the job array of 0 to 126 (zero indexing in python).

##Load anaconda python packages
module load anaconda3 
##Run script
srun python3 scripts/09_crowdsourced_interventions.py $SLURM_ARRAY_TASK_ID 
