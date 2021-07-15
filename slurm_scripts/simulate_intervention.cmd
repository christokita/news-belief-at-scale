#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=60gb
#SBATCH --array=0-27 #we have 28 fake news articles
##SBATCH --array=17 #only the largest fake news article
#SBATCH --time=71:00:00
##SBATCH --time=120:00:00 #only largest fake news article
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ctokita@princeton.edu
#SBATCH --output=slurm_outfiles/slurm-%A_%a.out

## NOTE: We have 29 unique fake news stories stories as of now, hence the job array of 0 to 28 (zero indexing in python).

##Load anaconda python packages
module load anaconda3 
##Run script
srun python3 scripts/08_model_interventions.py $SLURM_ARRAY_TASK_ID 
