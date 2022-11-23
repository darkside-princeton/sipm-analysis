import subprocess

class Scheduler():
    def __init__(self, script="sipm/recon/analysis.py", dirs="", num=1000000000):
        self.num = num
        self.dirs = dirs
        self.script = script
        self.scratch = "/scratch/gpfs/GALBIATI/"
        self.partition = "physics"
        self.nodes = 1 
        self.tasks_per_node = 1 
        self.cpu_memory = 4 # in GB
        self.wall_time = "1:10:00"
        self.output = "test3.log"
        self.conda_module = "module load anaconda3/2021.11"
        self.conda_activate = "conda activate ds-pu"

    def submit(self):
        for index,directory in enumerate(self.dirs):
            self.batch_script(index,directory)
            cmd = f"sbatch {self.scratch}/job_{index}.sh"
            p = subprocess.Popen(cmd.split(' '), stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    def batch_script(self, index, directory):
        with open(f"{self.scratch}/job_{index}.sh", "w") as f:
            f.write("#!/bin/bash -l\n")
            f.write(f"#SBATCH --partition {self.partition}\n")
            f.write(f"#SBATCH --nodes {self.nodes}\n")
            f.write(f"#SBATCH --ntasks-per-node {self.tasks_per_node}\n")
            f.write(f"#SBATCH --mem-per-cpu {self.cpu_memory}G\n")
            f.write(f"#SBATCH --time {self.wall_time}\n")
            f.write(f"#SBATCH --job-name job_{index}\n")
            f.write(f"#SBATCH --output {self.scratch}/log_{index}.log\n\n")
            f.write("now=$(date +"%Y%m%d")\n\n")
            f.write(f"module load anaconda3/2021.11\n")
            f.write(f"conda activate ds-pu\n\n")
            f.write(f"python {self.script} -f {directory} -n {self.num}\n")









