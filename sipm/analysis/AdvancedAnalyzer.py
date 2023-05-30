from typing import List
import glob
import pandas as pd
import yaml

class AdvancedAnalyzer:
    """The parent class for the post-analysis code.
    """
    def __init__(self, dir:str,voltages:List[float], channels:List[int], positions:List[str]):
        """Constructor.

        Args:
            dir (str): Directory containing the hdf5 files (e.g. '/scratch/gpfs/as111/results/') \n
            volts (List[float]): SiPM bias voltage (e.g. [63,65,67,69,71]) \n
            channels (List[int]): Channel number (e.g. [0,1,2,3]) \n
            positions (List[str]): Position (e.g. ['top','bottom'])
        """
        self.directory = dir
        self.data = {}
        self.proto_data = {}
        self.results = {}
        self.voltages = voltages
        self.channels = channels
        self.positions = positions

    def print_info(self):
        print('Data structure:\n','level1=position','level2=channel','level3=voltage')
        print(yaml.dump(self.proto_data,default_flow_style=False))

    def load_files(self, date:str, light:str, coinc:str, cond:str, run:str, wf:bool, verbose:bool):
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
        for pos in self.positions:
            self.data[pos] = {}
            self.proto_data[pos] = {}
            for ch in self.channels:
                self.data[pos][ch] = {}
                self.proto_data[pos][ch] = {}
                for volt in self.voltages:
                    self.data[pos][ch][volt] = []
                    self.proto_data[pos][ch][volt] = {}
                    if wf:
                        files = glob.glob(f"{self.directory}{date}/{date}_volt_{volt}_pos_{pos}_light_{light}_coinc_{coinc}_cond_{cond}_run{run}_wf.h5")
                    else:
                        files = glob.glob(f"{self.directory}{date}/{date}_volt_{volt}_pos_{pos}_light_{light}_coinc_{coinc}_cond_{cond}_run{run}[!_wf].h5")
                    for f in files:
                        df = pd.read_hdf(f, key=f'{volt}/{ch}')
                        self.data[pos][ch][volt].append(df)
                    if len(files)>0:
                        self.data[pos][ch][volt] = pd.concat(
                            self.data[pos][ch][volt], ignore_index=True)
                    for k in self.data[pos][ch][volt].keys():
                        self.proto_data[pos][ch][volt][k] = []
                    if verbose and ch==self.channels[0]:
                        if len(files)==0:
                            events = 0
                        else:
                            keys = self.data[pos][ch][volt].keys()
                            events = len(self.data[pos][ch][volt][keys[0]])
                        if wf:
                            print(f'{date} {volt}V {pos} {light} coinc={coinc} cond={cond} run{run} wf - {len(files)} files {events} events')
                        else:
                            print(f'{date} {volt} {pos} {light} coinc={coinc} cond={cond} run{run} - {len(files)} files {events} events')