#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=0-99
#SBATCH --time=11:59:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ctokita@princeton.edu
#SBATCH --output=slurm_outfiles/slurm-%A_%a.out

##Load anaconda python packages
module load anaconda3 
##Run script passing the array job number to the script
srun python3 scripts/03a_parse_fakenews_spreaders.py $SLURM_ARRAY_TASK_ID 
