import datetime
import glob
import h5py
import numpy as np
import matplotlib.pyplot as plt
from source import waveform as wvf

class Dataset: 
    def __init__(self,  Path, ShowPlots=True, Selection='*', Pol=1, NumChannels=2):
        self.Path = Path
        self.NumChannels = NumChannels
        self.ShowPlots = ShowPlots
        self.Selection = Selection
        self.Pol = Pol 
        self.Ch = self.InitializeChannels(self.NumChannels, self.Pol)
        self.Files = glob.glob(self.Path+self.Selection)

    def RunStandardAnalysis(self, NoiseDataset=None): 
        for File in self.Files: 
            self.ImportDataFromHDF5(File, self.Ch)
        self.DoAnalysis(self.Ch, NoiseDataset=NoiseDataset)
        self.ChargeCollection = self.Ch[0].Max / self.Ch[1].Max
        self.DiffMinute = int((np.max(self.Ch[0].TimeStamp) - np.min(self.Ch[0].TimeStamp)).seconds/60.0 + 0.5)
        self.XTicks = int((self.DiffMinute/12.0 + 0.5))+1
        self.NoiseCut = 1000
        self.Cut = np.where(self.Ch[0].BaseStd < self.NoiseCut)
        self.InverseCut = np.where(self.Ch[0].BaseStd > self.NoiseCut)

    def InitializeChannels(self, NumChannels=2, Pol=1):
        return [wvf.Waveform(ID=ii, Pol=(-1)**ii*-1*Pol) for ii in range(1,NumChannels+1)]

    def ImportDataFromHDF5(self, File, channels, var=['trig','timestamp']):
        f = h5py.File(File, 'r')  
        for ch in channels:
            ch.Time = np.array(f.get('Time')).flatten() * ch.TScale
            if 'trig' in var:
                ch.Trigger = np.array(f.get('Trigger')).flatten() * ch.VScale
            Group = f.get(ch.ChName)
            GroupKeys = Group.keys()
            ch.Files.append(len(GroupKeys))
            for key in GroupKeys:
                ch.Amp.append(np.array(Group.get(key)).flatten() * ch.VScale * ch.Pol)
                if "timestamp" in var:
                    ch.TimeStamp.append(datetime.datetime.strptime(Group.get(key).attrs["TimeStamp"].decode('utf-8'), '%Y%m%d%H%M%S'))
        f.close()
            

    def DoAnalysis(self, channels, NoiseDataset=None):
    ###### Basic analysis: baseline subtraction, waveform averaging, obtaining fourier spectra, frequency bandpass filter and finding extrema.
        Print = False 
        for ii, ch in enumerate(channels):
            # print(" | Processing data in channel %d..." % (ch.ID))
            ch.GetSampling()
            ch.Amp = [x for _, x in sorted(zip(ch.TimeStamp, ch.Amp))]
            ch.Amp = np.array(ch.Amp)
            ch.TimeStamp = np.array(sorted(ch.TimeStamp))

            # ch.TimeStamp = np.array(ch.TimeStamp)
            ch.Amp = ch.SubtractBaseline(Data=ch.Amp, state=Print)
            ch.Amp = ch.RemoveNoise(Data=ch.Amp, HighPass=80000, state=Print)
            if NoiseDataset is not None:
                for jj,amp in enumerate(ch.Amp):
                    ch.Amp[jj] =  ch.Amp[jj]-np.mean(NoiseDataset.Ch[ii].Amp,axis=0)

            # ch.RunFit(Data=ch.Amp)
            ch.GetAllMaxima(Data=ch.Amp, state=Print)
            # ch.FindMaxGradient(Data=ch.Amp ,state=Print)
            ch.GetDriftTime(Data=ch.Amp, Threshold=0.1)
            ch.GetIntegral(Data=ch.Amp, state=Print)
            # ch.GetBaselineNoise(Data=ch.Amp)

    def RoundUpToNext(self, Num, Ceil): 
        return int(np.ceil(Num / float(Ceil))) * float(Ceil)

    def RoundDownToNext(self, Num, Floor): 
        return int(np.floor(Num / float(Floor))) * float(Floor)