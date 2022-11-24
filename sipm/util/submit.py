import subprocess

class Scheduler():
    """Class for managing and submitting jobs for a large number of analyses to the Della cluster at Princeton.

    Atrributes
    ----------
    dirs : list of string
        List of directories containging the wavedump files to be analyzed.
    script : string, optional
        The analysis python script that will run the default analysis on all files. By default this is located under `sipm/recon/analysis.py`.
    num : string, optional
        Number of waveforms for each file to analyze. Default is set to 1e9, which effectively means that the entire file will be read in and analyzed. 
    """
    def __init__(self, dirs, script="sipm/recon/analysis.py",  num=1000000000):
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
        """For each subdirectory a new shell script is generated and then executed by submitting it to the cluster.
        """
        for index,directory in enumerate(self.dirs):
            self.batch_script(index,directory)
            cmd = f"sbatch {self.scratch}/job_{index}.sh"
            p = subprocess.Popen(cmd.split(' '), stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    def batch_script(self, index, directory):
        """Generates shell script for submitting slurm jobs to the Della cluster at Princeton and saves them under the common directory defined by `self.scratch`.

        Parameters
        ----------
        index : int
            The job index for naming the batch script.
        directory : string
            Directory containg the wavedump files for this analysis.
        """
        with open(f"{self.scratch}/job_{index}.sh", "w") as f:
            f.write("#!/bin/bash -l\n")
            f.write(f"#SBATCH --partition {self.partition}\n")
            f.write(f"#SBATCH --nodes {self.nodes}\n")
            f.write(f"#SBATCH --ntasks-per-node {self.tasks_per_node}\n")
            f.write(f"#SBATCH --mem-per-cpu {self.cpu_memory}G\n")
            f.write(f"#SBATCH --time {self.wall_time}\n")
            f.write(f"#SBATCH --job-name job_{index}\n")
            f.write(f"#SBATCH --output {self.scratch}/log_{index}.log\n\n")
            f.write(f"module load anaconda3/2021.11\n")
            f.write(f"conda activate ds-pu\n\n")
            f.write(f"python {self.script} -f {directory} -n {self.num}\n")