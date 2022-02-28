#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=24gb
#SBATCH --array=0-135
#SBATCH --time=23:59:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ctokita@princeton.edu
#SBATCH --output=slurm_outfiles/slurm-%A_%a.out

## NOTE: We have 136 unique stories as of now, hence the job array of 0 to 135 (zero indexing in python).

##Load anaconda python packages
module load anaconda3/2021.11 
conda activate my_conda
##Run script
srun python3 scripts/05a_calculate_exposure.py $SLURM_ARRAY_TASK_ID 
