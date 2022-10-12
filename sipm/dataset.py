import numpy as np
import yaml 

import sipm.sipm as sipm

class Dataset: 
    def __init__(self,  path, pol=1, channels=4, samples=3000, bias=65, pos='top'):
        self.path = path
        self.channels = channels
        self.samples = samples
        self.pol = pol
        self.ch = self.InitializeChannels()
        self.gain = []
        self.summed_integral_pe = []
        if pos=='top':
            slope = np.array([43.893, 39.983, 41.399, 41.417])
            vbd = np.array([54.378, 54.880, 54.436, 55.156])
            self.gain = slope*(bias-vbd)
        if pos=='bottom':
            slope = np.array([39.327, 41.039, 37.207, 43.592])
            vbd = np.array([54.766, 55.099, 55.224, 54.624])
            self.gain = slope*(bias-vbd)
        
    def InitializeChannels(self):
        channels = []
        for i in self.channels:
            new_channel = sipm.SiPM(id=i, pol=self.pol, path=self.path, samples=self.samples)
            channels.append(new_channel)
        return channels

    def get_summed_integral_pe(self):
        self.summed_integral_pe = np.zeros(self.ch[0].cumulative_nevents)
        for i in self.channels:
            self.summed_integral_pe += np.array(self.ch[i].integral)/self.gain[i]

    def get_waveforms_id(self, count=10, integral_range=(0,0)):
        count_ = 0
        event_id = []
        ev = 0
        while count_ < count and ev < self.ch[0].nevents:
            if self.summed_integral_pe[ev]<integral_range[1] and self.summed_integral_pe[ev]>integral_range[0]:
                event_id.append(ev)
                count_ += 1
            ev += 1
        return event_id