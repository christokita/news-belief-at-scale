#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=18gb
#SBATCH --array=0-299
#SBATCH --time=11:59:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ctokita@princeton.edu
#SBATCH --output=slurm_outfiles/slurm-%A_%a.out

##Load anaconda python packages
module load anaconda3/2021.11 
conda activate my_conda 
##Run script passing the array job number to the script
srun python3 scripts/02a_parse_users.py $SLURM_ARRAY_TASK_ID 
