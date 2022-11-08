import numpy as np
import yaml 

import sipm.sipm as sipm

class Dataset: 
    def __init__(self,  path, pol=-1, channels=4, samples=3000):
        self.datadir = "/scratch/gpfs/GALBIATI/data/sipm/reflector_studies/"
        self.path = f"{self.datadir}/{path}"
        self.channels = channels
        self.samples = samples
        self.pol = pol
        self.ch = self.InitializeChannels()
        
    def InitializeChannels(self):
        channels = []
        for i in self.channels:
            new_channel = sipm.SiPM(id=i, pol=self.pol, path=self.path, samples=self.samples)
            channels.append(new_channel)
        return channels
    
    def analyze(self, header=True, num_events=1e9):
        for i in self.channels:
            self.ch[i].read_data(header=header, num_events=num_events)
            self.ch[i].bandpass_filter(low=1, high=8e6, order=1)
            self.ch[i].baseline_subtraction()
            self.ch[i].get_max()
            self.ch[i].get_integral()
    
    def calibrate(self):
        for i in self.channels:
            self.ch[i].calibrate(self.ch[i].peak, width=10, verbose=1)
            self.ch[i].plot_calibration(self.ch[i].peak)
    
    def clear(self):
        for i in self.channels:
            self.ch[i].clear()
        