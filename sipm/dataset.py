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
        self.gain = [263.83304269154604, 225.54550142289833, 246.43480002250982, 216.58835157109758] #Gain of each channel of top tile. Hard-coded for now
        
    def InitializeChannels(self):
        channels = []
        for i in self.channels:
            new_channel = sipm.SiPM(id=i, pol=self.pol, path=self.path, samples=self.samples)
            channels.append(new_channel)
        return channels

    def get_summed_integral_pe(self):
        self.summed_integral_pe = np.zeros(self.ch[0].nevents)
        for i in self.channels:
            self.summed_integral_pe += np.array(self.ch[i].integral)/self.gain[i]