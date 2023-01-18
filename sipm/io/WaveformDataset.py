import numpy as np
import sipm.io.WaveformAnalyzer as wfa

class Dataset: 
    def __init__(self,  path, pol=-1, channels=[0,1,2,3], samples=3000):
        self.datadir = "/scratch/gpfs/GALBIATI/data/sipm/reflector_studies/"
        self.path = f"{self.datadir}/{path.replace(self.datadir, '')}"
        self.channels = channels
        self.samples = samples
        self.pol = pol
        self.ch = self.InitializeChannels()
        
    def InitializeChannels(self):
        channels = []
        for i in self.channels:
            new_channel = wfa.WaveformAnalyzer(id=i, pol=self.pol, path=self.path, samples=self.samples)
            channels.append(new_channel)
        return np.array(channels)
    
    def analyze(self, header=True, num_events=1e9, clear=True, sum=False):
        for i in self.channels:
            self.ch[i].read_data(header=header, num_events=num_events)
            self.ch[i].bandpass_filter(low=1, high=8e6, order=1)
            self.ch[i].baseline_subtraction()
        if sum:
            self.sum_traces(clear=clear)
            self.sum.get_max()
            self.sum.get_integral()
        for i in self.channels:
            self.ch[i].get_max()
            self.ch[i].get_integral()
        if clear:
            self.clear()
            self.sum.clear()

    def sum_traces(self, clear=True):
        self.sum = wfa.WaveformAnalyzer(id=-1, pol=self.pol, path=self.path, samples=self.samples)
        self.sum.traces = np.array([self.ch[i].traces for i in self.channels])
        self.sum.traces = np.sum(self.sum.traces, axis=0)
    
    def calibrate(self):
        for i in self.channels:
            self.ch[i].calibrate(self.ch[i].peak, width=10, verbose=1)
            self.ch[i].plot_calibration(self.ch[i].peak)
    
    def clear(self):
        for i in self.channels:
            self.ch[i].clear()
    
        