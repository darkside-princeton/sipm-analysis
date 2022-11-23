#!/bin/bash -l
#SBATCH --partition physics
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem-per-cpu 32G
#SBATCH --time 3:00:00
#SBATCH --job-name jupyter-notebook
#SBATCH --output jupyter-notebook-%J.log

# get tunneling info
XDG_RUNTIME_DIR=""
port=$(shuf -i8000-9999 -n1)
node=$(hostname -s)
user=$(whoami)
cluster=$(hostname -f | awk -F"." '{print $2}')

# load modules or conda environments here
module load anaconda3/2021.11
conda activate myenv


jupyter-notebook --no-browser --port=8888 --ip=${node}