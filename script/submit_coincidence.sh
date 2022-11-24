#!/bin/bash
#SBATCH --job-name=coincidence

#SBATCH --output=log/coincidence_%A_%a.out

#SBATCH --error=log/coincidence_%A_%a.err

#SBATCH --nodes=1

#SBATCH --ntasks=1

#SBATCH --cpus-per-task=1

#SBATCH --mem-per-cpu 16G

#SBATCH --time=01:30:00

#SBATCH --array=0-5

module load anaconda3/2021.11
conda activate ds-pu

python root_coincidence.py $SLURM_ARRAY_TASK_ID