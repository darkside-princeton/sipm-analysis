import h5py
import pandas as pd
import re
import pwd
import os
import time
import numpy as np

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
        self.d = dataset
        self.filename = self.d.path
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

    def get_metadata(self,path):
        self.tag_list = ['volt', 'pos', 'light', 'coinc', 'cond']
        tag_dict = {}
        for tag in self.tag_list:
            index = path.find(tag)
            filename = path[index:].split('/')[0]
            indices_object = re.finditer(pattern='_', string=filename)
            indices = [index.start() for index in indices_object]
            if len(indices) < 2:
                tag_dict[tag] = filename[indices[0]+1:]
            else:
                tag_dict[tag] = filename[indices[0]+1:indices[1]]
        return tag_dict
    
    def save(self):

        if not os.path.isdir(self.scratch):
            os.makedirs(f"{self.scratch}")

        for i in self.d.channels:
            data = {} 
            data['amplitude'] = self.d.ch[i].peak
            data['integral'] = self.d.ch[i].integral
            data['avg_waveform'] = np.mean(self.d.ch[i].traces, axis=0)
            data['time'] = self.d.ch[i].time
            df = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in data.items() ]))

            metadata = self.get_metadata(self.filename)

            tag = ""
            for x,y in metadata.items():
                tag += f"_{x}_{y}"

            counter = 0
            while True:
                if counter > 10:
                    break
                try:
                    store = pd.HDFStore(f"{self.scratch}/{self.date}{tag}.h5")
                    store.put(f"{metadata['volt']}/{i}", df, format='t', append=True, data_columns=True)
                    store.get_storer(f"{metadata['volt']}/{i}").attrs.metadata = metadata
                    store.close()
                    break
                except:
                    time.sleep(np.random.randint(0,10))
                    counter += 1




