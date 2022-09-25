import numpy as np

import sipm.sipm as sipm

class Dataset: 
    def __init__(self,  path, pol=1, channels=4):
        self.path = path
        self.channels = channels
        self.pol = pol
        self.ch = self.InitializeChannels(self.channels, self.pol)

    def InitializeChannels(self, num_channels=2, pol=1):
        return [sipm.SiPM(id=ii, pol=pol, path=self.path) for ii in self.channels]