#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=40gb
#SBATCH --array=0-299
#SBATCH --time=23:59:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ctokita@princeton.edu
#SBATCH --output=slurm_outfiles/slurm-%A_%a.out

##Load anaconda python packages
module purge
module load anaconda3/2021.11 
conda activate my_conda
##Run script passing the array job number to the script
srun python3 scripts/04b_follower_ideologies_inference.py $SLURM_ARRAY_TASK_ID 
