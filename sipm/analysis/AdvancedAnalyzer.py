from typing import Dict
import glob
import pandas as pd
import yaml
import numpy as np

class AdvancedAnalyzer:
    """The parent class for the post-analysis classes.
    """
    def __init__(self, directory:str, metadata_dict:Dict, script:str, merge:bool, verbose:bool):
        """AdvancedAnalyzer constructor.

        Args:
            directory (str): Directory containing processed HDF5 files (e.g. '/scratch/gpfs/as111/results/')
            metadata_dict (Dict): Metadata arranged into a nested dictionary. Need to match the file names. (See jupyter/calibration_liq2.ipynb for example)
            script (str): The name of the pre-processing script (without .py)
            merge (bool): Whether to merge different runs.
            verbose (bool): Whether to print more information.
        """
        self.directory = directory
        self.metadata = metadata_dict
        self.data = {}
        self.proto_data = {}
        self.script = script
        self.merge = merge
        self.verbose = verbose
        self.baseline = {}

    def print_info(self):
        """Print nested structure of self.data with every lowest-level entry replaced by '###'.
        """
        print('\nData structure:')
        print(yaml.dump(self.proto_data,default_flow_style=False))

    def load_file(self, date:str, volt:float, pos:str, light:str, coinc:str, cond:str, run:int, ch:int):
        """Load the file specified by the metadata in the arguments. Can use '*' for pattern matching.

        Args:
            date (str): Date (e.g. '2022-11-22')
            volt (float): Voltage (e.g. 67)
            pos (str): Position (e.g. 'top')
            light (str): Light (e.g. 'scintillation')
            coinc (str): Coincidence (e.g. '111')
            cond (str): Condition (e.g. 'calibration')
            run (int): Run number (e.g. 0)
            ch (int): Channel number (e.g. 0)

        Returns:
            pandas.DataFrame: data
        """
        data = []
        files = glob.glob(f"{self.directory}{date}/{date}_volt_{volt}_pos_{pos}_light_{light}_coinc_{coinc}_cond_{cond}_run{run}_{self.script}.h5")
        for f in files:
            ind = f.find('run')+3
            run_number = int(f[ind:f.find('_',ind)])
            df = pd.read_hdf(f, key=f'{volt}/{ch}')
            df['event'] = df.index
            df['run'] = np.array([run_number]*df.shape[0])
            data.append(df)
        if len(files)>0 and self.merge:
            data = pd.concat(data, ignore_index=True).sort_values(by=['run','event']).reset_index(drop=True)
        if self.verbose:
            if len(files)==0:
                events = 0
            else:
                keys = data.keys()
                events = len(data[keys[0]])
            print(f'{date} {volt}V {pos} {light} coinc={coinc} cond={cond} run{run} ch{ch} script={self.script}.py - {len(files)} files {events} events')
        return data

    def load_files(self):
        """Load processed HDF5 files specified by self.metadata.
        """
        self.metadata_dfs(self.metadata, self.data, self.proto_data)
        if self.verbose:
            self.print_info()

    def metadata_dfs(self, mydict:Dict, mydata:Dict, myproto:Dict):
        """Loop over and load all files recursively.

        Args:
            mydict (Dict): metadata dictionary
            mydata (Dict): data dictionary
            myproto (Dict): proto-data dictionary (same structure as data but with data content replaced by '###')
        """
        if 'metadata' not in mydict:
            for k,v in mydict.items():
                mydata[k] = {}
                myproto[k] = {}
                self.metadata_dfs(v,mydata[k],myproto[k])
        else:
            m = mydict['metadata']
            mydata['metadata'] = mydict['metadata']
            mydata['data'] = self.load_file(m['date'],m['volt'],m['pos'],m['light'],m['coinc'],m['cond'],m['run'],m['ch'])
            for k in mydata['data'].keys():
                myproto[k] = '###'
                
    def baseline_analysis(self, threshold_dict:Dict):
        """Apply baseline quality cuts and calculate cut fraction.

        Args:
            threshold_dict (Dict): Baseline rms thresholds arranged into a nested dictionary with the same structure as self.metadata.
        """
        self.bsl_rms_thre = threshold_dict
        self.baseline_cut_dfs(threshold_dict, self.data, self.baseline)
    
    def baseline_cut_dfs(self, mythre:Dict, mydata:Dict, mybsl:Dict):
        """Loop over all data entry and apply baseline rms cut.

        Args:
            mythre (Dict): Baseline rms thresholds arranged into a nested dictionary with the same structure as self.metadata.
            mydata (Dict): Data dictionary
            mybsl (Dict): Baseline dictionary. Contains baseline rms and baseline mean distributions as well as cut fractions.
        """
        if 'data' not in mydata:
            for k,v in mydata.items():
                mybsl[k] = {}
                self.baseline_cut_dfs(mythre[k],v,mybsl[k])
        else:
            data = mydata['data']
            data['bsl_cut'] = data['baseline_rms']<mythre
            mybsl['rms_threshold'] = mythre
            mybsl['rms_hist'], mybsl['rms_bins'] = np.histogram(data['baseline_rms'], bins=500, range=(0, 5))
            mybsl['mean_hist'], mybsl['mean_bins'] = np.histogram(data['baseline_mean'], bins=1500, range=(3720, 3850))
            mybsl['mean_hist_cut'], mybsl['mean_bins_cut'] = np.histogram(data['baseline_mean'].loc[data['bsl_cut']], bins=1500, range=(3720, 3850))
            mybsl['cut_fraction'] = 1-np.sum(data['bsl_cut'])/data.shape[0]
            if self.verbose:
                print(f'{mydata["metadata"]["pos"]} ch{mydata["metadata"]["ch"]} {mydata["metadata"]["volt"]}V cut fraction = {mybsl["cut_fraction"]*100:.5f}%')