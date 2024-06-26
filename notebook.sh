#!/bin/bash -l
#SBATCH --partition physics
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem-per-cpu 16G
#SBATCH --time 4:00:00
#SBATCH --job-name jupyter-notebook
#SBATCH --output jupyter-notebook-%J.log

# delete all previous log files except the newly generated one for this job
find . -type f  -name "jupyter-notebook-*.log" ! -name "jupyter-notebook-${SLURM_JOB_ID}.log" -delete

# get tunneling info
XDG_RUNTIME_DIR=""
port=$(shuf -i8000-9999 -n1)
node=$(hostname -s)
user=$(whoami)
cluster=$(hostname -f | awk -F"." '{print $2}')

# load modules or conda environments here
conda activate sipm

jupyter-notebook --no-browser --port=8888 --ip=${node}