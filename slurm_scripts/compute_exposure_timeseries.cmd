#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=12gb
#SBATCH --array=0-135
#SBATCH --time=36:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ctokita@princeton.edu
#SBATCH --output=slurm_outfiles/slurm-%A_%a.out

## NOTE: We have 136 unique stories as of now, hence the job array of 0 to 135 (zero indexing in python).

##Load anaconda python packages
module load anaconda3 
##Run script
srun python3 scripts/05a_exposure_over_time.py $SLURM_ARRAY_TASK_ID 
