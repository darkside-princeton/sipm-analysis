import numpy as np
import sipm.recon.WaveformAnalyzer as wfa

class WaveformDataset: 
    def __init__(self,  path, pol=-1, channels=[0,1,2,3], samples=3000):
        """Class that analyzes the waveform data. Contains a list of WaveformAnalyzer objects for all the channels in a SiPM array.

        Args:
            path (str): Data folder path
            pol (int, optional): Polarity. Defaults to -1.
            channels (list, optional): A list of channel numbers. Defaults to [0,1,2,3].
            samples (int, optional): Length of acquisition window. Defaults to 3000.
        """
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

    def process_laser_pulses(self, header=True, num_events=1e9):
        """Obtain pulse information from laser data.

        Args:
            header (bool, optional): Whether wavedump file contains header. Defaults to True.
            num_events (int, optional): Number of events. Defaults to 1e9.
        """
        for i in self.channels:
            self.ch[i].read_data(header=header, num_events=num_events)
            self.ch[i].baseline_subtraction(samples=self.ch[i].trigger_position-int(0.5/self.ch[i].sample_step)) # average maximum position - 0.5us
            self.ch[i].ar_filter(tau=20) # 20 samples = 80us = fast component
            self.ch[i].get_max(ar=True, trig=True) # AR matched filter, maximum near trigger position
            self.ch[i].get_integral(length_us=[5]) # 5us integral
        self.clear()

    def process_laser_waveforms(self, header=True, num_events=1e9, calib=""):
        """Obtain single PE average waveform from laser data. Only a placeholder for now.

        Args:
            header (bool, optional): Whether wavedump file contains header. Defaults to True.
            num_events (int, optional): Number of events. Defaults to 1e9.
            calib (str, optional): Directory that contains SiPM calibration results. Defaults to "".
        """
        pass

    def process_scintillation_pulses(self, header=True, num_events=1e9, calib=""):
        """Obtain pulse information from scintillation data. Only a placeholder for now.

        Args:
            header (bool, optional): Whether wavedump file contains header. Defaults to True.
            num_events (int, optional): Number of events. Defaults to 1e9.
            calib (str, optional): Directory that contains SiPM calibration results. Defaults to "".
        """
        pass

    def process_scintillation_waveforms(self, header=True, num_events=1e9, calib="", fprompt=[0.1,0.6], pe=[300,700]):
        """Obtain average waveform from scintillation data. Only a placeholder for now.

        Args:
            header (bool, optional): Whether wavedump file contains header. Defaults to True.
            num_events (_type_, optional): Number of events. Defaults to 1e9.
            calib (str, optional): Directory that contains SiPM calibration results. Defaults to "".
            fprompt (list, optional): Range of f_prompt. Defaults to [0.1,0.6].
            pe (list, optional): Range of PEs. Defaults to [300,700].
        """
        pass

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
    
        