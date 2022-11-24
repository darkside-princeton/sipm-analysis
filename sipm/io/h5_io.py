import h5py
import pandas as pd
import re
import pwd
import os

class IO():
    """Class for handling the saving of high level variables form the SiPM analysis.

    Atrributes
    ----------
    filename : string
        List of directories containging the wavedump files to be analyzed.
    script : string, optional
        The analysis python script that will run the default analysis on all files. By default this is located under `sipm/recon/analysis.py`.
    num : string, optional
        Number of waveforms for each file to analyze. Default is set to 1e9, which effectively means that the entire file will be read in and analyzed. 
    """
    def __init__(self, dataset):
        self.filename = dataset.path
        self.date = self.get_date(self.filename)
        self.username = pwd.getpwuid(os.getuid())[0]
        self.scratch = f"/scratch/gpfs/{self.username}/results/{self.date}"

    def get_date(self, filename):
        """Extracts the date of data taking from the filepath.

        Parameters
        ----------
        filename : string
            The parent directory for all measurements.
        """
        indices_object = re.finditer(pattern='-', string=filename)
        indices = [index.start() for index in indices_object]
        return filename[indices[0]-4:indices[1]+3]

    def get_tag(self,filename):
        tags = ['volt', 'pos', 'light', 'coinc', 'cond']
        for tag in tags: 
            pass
    
    # def save(self):
    #     pd.




