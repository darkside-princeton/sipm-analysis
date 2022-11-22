#!/bin/bash
#SBATCH --job-name=spewf

#SBATCH --output=log/spewf_%A_%a.out

#SBATCH --error=log/spewf_%A_%a.err

#SBATCH --nodes=1

#SBATCH --ntasks=1

#SBATCH --cpus-per-task=1

#SBATCH --mem-per-cpu 16G

#SBATCH --time=01:30:00

#SBATCH --array=0-9

module load anaconda3/2021.11
conda activate ds-pu

python root_spewf.py $SLURM_ARRAY_TASK_ID