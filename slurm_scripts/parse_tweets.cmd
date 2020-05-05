#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=11:59:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ctokita@princeton.edu
#SBATCH --output=slurm_outfiles/slurm-%A.out

##Load anaconda python packages
module load anaconda3 
##Run script
srun python3 scripts/01_parse_tweets.py
