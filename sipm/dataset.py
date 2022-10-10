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
        self.gain = [325.283, 306.669, 324.455, 330.207] #Gain of each channel of top tile @65V. Hard-coded for now
        
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