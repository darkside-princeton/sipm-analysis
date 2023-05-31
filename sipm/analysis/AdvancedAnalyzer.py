from typing import Dict
import glob
import pandas as pd
import yaml
import numpy as np

class AdvancedAnalyzer:
    """The parent class for the post-analysis code.
    """
    def __init__(self, directory:str, metadata_dict:Dict, wf:bool, merge:bool, verbose:bool):
        """Constructor.

        Args:
            dir (str): Directory containing the hdf5 files (e.g. '/scratch/gpfs/as111/results/') \n
            volts (List[float]): SiPM bias voltage (e.g. [63,65,67,69,71]) \n
            channels (List[int]): Channel number (e.g. [0,1,2,3]) \n
            positions (List[str]): Position (e.g. ['top','bottom'])
        """
        self.directory = directory
        self.metadata = metadata_dict
        self.data = {}
        self.proto_data = {}
        self.results = None
        self.wf = wf
        self.merge = merge
        self.verbose = verbose
        self.baseline = {}

    def print_info(self):
        print('\nData structure:')
        print(yaml.dump(self.proto_data,default_flow_style=False))

    def load_file(self, date:str, volt:float, pos:str, light:str, coinc:str, cond:str, run:int, ch:int):
        data = []
        if self.wf:
            files = glob.glob(f"{self.directory}{date}/{date}_volt_{volt}_pos_{pos}_light_{light}_coinc_{coinc}_cond_{cond}_run{run}_wf.h5")
        else:
            files = glob.glob(f"{self.directory}{date}/{date}_volt_{volt}_pos_{pos}_light_{light}_coinc_{coinc}_cond_{cond}_run{run}[!_wf].h5")
        for f in files:
            ind = f.find('run')+3
            if self.wf:
                run_number = int(f[ind:f.find('_',ind)])
            else:
                run_number = int(f[ind:f.find('.',ind)])
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
            if self.wf:
                print(f'{date} {volt}V {pos} {light} coinc={coinc} cond={cond} run{run} ch{ch} wf - {len(files)} files {events} events')
            else:
                print(f'{date} {volt}V {pos} {light} coinc={coinc} cond={cond} run{run} ch{ch} - {len(files)} files {events} events')
        return data

    def load_files(self):
        """Load pre-processed hdf5 files. Arguments need to match the file name substrings. Can use '*' for pattern matching except for the 'wf' argument.

        Args:
            date (str): Date (e.g. '2023-05-10') \n
            light (str): Light source (e.g. 'scintillation') \n
            coinc (str): Coincidence setting (e.g. '000') \n
            cond (str): Condition (e.g. 'gamma') \n
            run (str): Run number (e.g. '0') \n
            wf (bool): Processed by laser_pulse.py or scintillation_pulse.py = True; processed by laser_waveform.py or scintillation_waveform.py = False. \n
            verbose (bool): Print messages \n
        """
        self.metadata_dfs(self.metadata, self.data, self.proto_data)
        if self.verbose:
            self.print_info()

    def metadata_dfs(self, mydict, mydata, myproto):
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
                
    def baseline_analysis(self, threshold_dict):
        self.bsl_rms_thre = threshold_dict
        self.baseline_cut_dfs(threshold_dict, self.data, self.baseline)
    
    def baseline_cut_dfs(self, mythre, mydata, mybsl):
        if 'data' not in mydata:
            for k,v in mydata.items():
                mybsl[k] = {}
                self.baseline_cut_dfs(mythre[k],v,mybsl[k])
        else:
            data = mydata['data']
            data['bsl_cut'] = data['baseline_rms']<mythre
            mybsl['rms_threshold'] = mythre
            nbins = 500
            range_min = 0
            range_max = 5
            mybsl['rms_hist'], mybsl['rms_bins'] = np.histogram(data['baseline_rms'], bins=nbins, range=(range_min, range_max))
            nbins = 1500
            range_min = 3720
            range_max = 3850
            mybsl['mean_hist'], mybsl['mean_bins'] = np.histogram(data['baseline_mean'], bins=nbins, range=(range_min, range_max))
            mybsl['mean_hist_cut'], mybsl['mean_bins_cut'] = np.histogram(data['baseline_mean'].loc[data['bsl_cut']], bins=nbins, range=(range_min, range_max))
            mybsl['cut_fraction'] = 1-np.sum(data['bsl_cut'])/data.shape[0]
            if self.verbose:
                print(f'{mydata["metadata"]["pos"]} ch{mydata["metadata"]["ch"]} {mydata["metadata"]["volt"]}V cut fraction = {mybsl["cut_fraction"]*100:.5f}%')