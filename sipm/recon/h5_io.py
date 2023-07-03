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
        The analysis python script that will run the default analysis on all files. 
        By default this is located under `sipm/exe/analysis.py`.
    num : string, optional
        Number of waveforms for each file to analyze. 
        Default is set to 1e9, which effectively means that the entire file will be read in and analyzed. 
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

    def get_metadata(self, path):
        """Extracts the metadata about the measurement parameters from the file path.
        They are being saved to the HDF5 file for future reference. 

        Parameters
        ----------
        path : string
            The directory containing the waveform files.
        """
        self.tag_list = ['volt', 'pos', 'light', 'coinc', 'cond', 'config', 'intensity']
        tag_dict = {}
        for tag in self.tag_list:
            index = path.find(tag)
            if index!=-1:
                filename = path[index:].split('/')[0]
                indices_object = re.finditer(pattern='_', string=filename)
                indices = [index.start() for index in indices_object]
                if len(indices) < 2 or tag=='cond':#condition description can contain multiple words separated by _
                    tag_dict[tag] = filename[indices[0]+1:]
                else:
                    tag_dict[tag] = filename[indices[0]+1:indices[1]]
        return tag_dict

    def get_run_number(self, path):
        """Extracts the run number by looking for the two closest '/' in the file path.

        Parameters
        ----------
        filename : string
            The directory containing the waveform files.
        """
        indices_object = re.finditer(pattern='/', string=path)
        indices = [index.start() for index in indices_object]
        arr = np.where((np.diff(indices)<4) & (np.diff(indices)>1))[0]
        if arr.shape[0]>0:
            run_number = arr[0]
            return path[indices[run_number]+1:indices[run_number+1]]
        else:
            return '0'
        

    def set_h5_filename(self,script):
        self.metadata = self.get_metadata(self.filename)
        self.run_number = self.get_run_number(self.filename)
        print("run number is: ", self.run_number)

        tag = ""
        for x,y in self.metadata.items():
            tag += f"_{x}_{y}"
        self.h5_filename = f"{self.scratch}/{self.date}{tag}_run{self.run_number}_{script}.h5"
    
    def save(self, script):
        """Saving the relevant parameters to a HDF5 file for further analysis. 
        Parameters
        ----------
        TODO: add possibility to decide which variables to save. 
        filename : string
            The directory containing the waveform files.
        """

        if not os.path.isdir(self.scratch):
            os.makedirs(f"{self.scratch}")
        self.set_h5_filename(script)
        # data output for individual channels
        for i in self.d.channels:
            data = self.d.ch[i].output
            print('keys: ',data.keys())
            
            df = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in data.items() ]))
            print('df=',df)
            store = pd.HDFStore(self.h5_filename)
            store.put(f"{self.metadata['volt']}/{i}", df, format='t', append=False, data_columns=True)
            store.get_storer(f"{self.metadata['volt']}/{i}").attrs.metadata = self.metadata
            store.close()

        # data output for all channels (e.g. total pe, fprompt)
        if self.d.output!={}:
            data = self.d.output
            print('keys: ',data.keys())
            
            df = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in data.items() ]))
            print('df=',df)

            store = pd.HDFStore(self.h5_filename)
            store.put(f"{self.metadata['volt']}/-1", df, format='t', append=False, data_columns=True)
            store.get_storer(f"{self.metadata['volt']}/-1").attrs.metadata = self.metadata
            store.close()
        
