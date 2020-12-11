#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=48gb
#SBATCH --array=0-28
#SBATCH --time=71:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ctokita@princeton.edu
#SBATCH --output=slurm_outfiles/slurm-%A_%a.out

## NOTE: We have 29 unique fake news stories stories as of now, hence the job array of 0 to 28 (zero indexing in python).

##Load anaconda python packages
module load anaconda3 
##Run script
srun python3 scripts/08a_model_interventions.py $SLURM_ARRAY_TASK_ID 
