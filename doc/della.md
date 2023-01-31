# Getting an Account

You can request an account on Della by sending an email to [cses@princeton.edu](mailto:cses@princeton.edu) providing the following information:

- Department: Physics
- Advisor: Cristiano Galbiati

# Logging in

```
ssh your-netid@della.princeton.edu
```

Your private data directory is ```/scratch/gpfs/your-netid/```, and the shared data directory of the group is ```/scratch/gpfs/GALBIATI/```.

# Running Jupyter Notebooks

On the Della cluster the login node is not meant for running extensive scripts or notebooks.

First you need to run `notebook.sh`. Make sure you adjust the parameters according to your needs and use your own conda environment. 

```bash
#!/bin/bash -l
#SBATCH --partition physics
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem-per-cpu 8G
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
conda activate ds-pu

jupyter-notebook --no-browser --port=8888 --ip=${node}
```

Once logged into Della you can submit the above job script via

```bash
sbatch notebook.sh
```

You can check if your job has been successfully submitted and is running via

```bash
squeue -u <netid>
```

If done correctly you should be seeing something like this: 

```bash
(ds-pu) [aj9512@della8 sipm-analysis]$ squeue -u aj9512
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
          42920715   physics jupyter-   aj9512  R       0:02      1 della-i13n24
```

In addition there should be a new `.log` file that was created for your job that contains information about the address of the Jupyter server that we started with our job script. It might look something like this

```bash
[I 2022-09-26 00:27:52.501 LabApp] JupyterLab extension loaded from /home/aj9512/.conda/envs/ds-pu/lib/python3.8/site-packages/jupyterlab
[I 2022-09-26 00:27:52.501 LabApp] JupyterLab application directory is /home/aj9512/.conda/envs/ds-pu/share/jupyter/lab
[I 00:27:52.504 NotebookApp] Serving notebooks from local directory: /home/aj9512/sipm-analysis
[I 00:27:52.505 NotebookApp] Jupyter Notebook 6.4.12 is running at:
[I 00:27:52.505 NotebookApp] http://della-i13n24:8888/?token=71acbd6139acb6da1a1ef02b731f44be051b689e149c8cd0
[I 00:27:52.505 NotebookApp]  or http://127.0.0.1:8888/?token=71acbd6139acb6da1a1ef02b731f44be051b689e149c8cd0
[I 00:27:52.505 NotebookApp] Use Control-C to stop this server and shut down all kernels (twice to skip confirmation).
[C 00:27:52.508 NotebookApp] 
    
    To access the notebook, open this file in a browser:
        file:///home/aj9512/.local/share/jupyter/runtime/nbserver-2175509-open.html
    Or copy and paste one of these URLs:
        http://della-i13n24:8888/?token=71acbd6139acb6da1a1ef02b731f44be051b689e149c8cd0
     or http://127.0.0.1:8888/?token=71acbd6139acb6da1a1ef02b731f44be051b689e149c8cd0
```

Let's copy the Jupyter address that starts with the name of the node where the job is running, i.e.

```bash
http://della-i13n24:8888/?token=71acbd6139acb6da1a1ef02b731f44be051b689e149c8cd0
```

Afterwards open any Juptyer notebook you would like to run and click on the bottom right on **Jupyter Server**. Then click on **Existing** from the drop-down menu, type in the server address we just copied, and press Enter.

You can double-check that you are actually connected to that Jupyter server by editing and running some cell and then saving the notebook and double-checking the content of the `.log` file from earlier. If you see new entries, it means that you are now correctly using Jupyter on Della via a batch job.