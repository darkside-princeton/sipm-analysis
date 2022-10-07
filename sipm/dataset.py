import numpy as np
import yaml 

import sipm.sipm as sipm

class Dataset: 
    def __init__(self,  path, pol=1, channels=4, samples=3000):
        self.path = path
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