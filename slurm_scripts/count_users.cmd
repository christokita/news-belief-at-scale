#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=23:59:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ctokita@princeton.edu
#SBATCH --output=slurm_outfiles/slurm-%A.out

##Load anaconda python packages
module load anaconda3/2021.11 
conda activate my_conda 
##Run script
srun python3 scripts/02b_count_users.py
