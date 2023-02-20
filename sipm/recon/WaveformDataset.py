import numpy as np
import sipm.recon.WaveformAnalyzer as wfa
import glob
import csv

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
        self.pos = ''
        self.volt = 0
        self.ch = self.InitializeChannels()
        self.output = {}
        
    def InitializeChannels(self):
        channels = []
        for i in self.channels:
            new_channel = wfa.WaveformAnalyzer(id=i, pol=self.pol, path=self.path, samples=self.samples)
            channels.append(new_channel)
        # Get position and voltage
        id_pos = [self.path.find('pos_')+4, self.path.find('_light')]
        self.pos = self.path[id_pos[0]:id_pos[1]]
        id_volt = [self.path.find('volt_')+5, self.path.find('_pos')]
        self.volt = int(self.path[id_volt[0]:id_volt[1]])
        return np.array(channels)

    def read_calibration(self,calib_dir):
        """Read calibration file to get SPE gain and SPE filtered amplitude range. The effective gain, taking into account DiCT and AP, is Q_peak*(1+Q_ap)/(1-p).

        Args:
            calib_dir (string): Directory of calibration csv files
        """
        self.gain = []
        self.a1min = []
        self.a1max = []
        print(f'Calibration csv file: {calib_dir}calibration_1122_{self.volt}V.csv')
        try:
            with open(f'{calib_dir}calibration_1122_{self.volt}V.csv') as f:
                r = csv.reader(f)
                line_count = 0
                for row in r:
                    if line_count>0:
                        if self.pos=='top':
                            if line_count<=4:
                                self.gain.append(float(row[7])*(1+float(row[9]))/(1-float(row[3])))
                                self.a1min.append(float(row[1]))
                                self.a1max.append(float(row[2]))
                        elif self.pos=='bottom':
                            if line_count>4:
                                self.gain.append(float(row[7])*(1+float(row[9]))/(1-float(row[3])))
                                self.a1min.append(float(row[1]))
                                self.a1max.append(float(row[2]))
                    line_count += 1
            print(f'Gain of {self.pos} SiPMs @{self.volt}V: {self.gain[0]:.2f} {self.gain[1]:.2f} {self.gain[2]:.2f} {self.gain[3]:.2f}')
            print(f'A1min of {self.pos} SiPMs @{self.volt}V: {self.a1min[0]:.2f} {self.a1min[1]:.2f} {self.a1min[2]:.2f} {self.a1min[3]:.2f}')
            print(f'A1max of {self.pos} SiPMs @{self.volt}V: {self.a1max[0]:.2f} {self.a1max[1]:.2f} {self.a1max[2]:.2f} {self.a1max[3]:.2f}')
        except:
            print('No calibration csv file. Use default gain of 500.')
            self.gain = [500]*4
            pass
    
    def get_total_pe(self):
        self.output['total_pe'] = np.zeros(self.ch[0].nevents)
        for ch in self.channels:
            self.output['total_pe'] += np.array(self.ch[ch].output['integral_5p00us'])/self.gain[ch]

    def get_fprompt(self):
        integral_long = np.zeros(self.ch[0].nevents)
        integral_short = np.zeros(self.ch[0].nevents)
        for ch in self.channels:
            integral_long += np.array(self.ch[ch].output['integral_5p00us'])
            integral_short += np.array(self.ch[ch].output['integral_0p50us'])
        self.output['fprompt'] = integral_short/integral_long
        integral_long = None
        integral_short = None

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
        self.read_calibration(calib)
        for i in self.channels:
            self.ch[i].read_data(header=header, num_events=num_events)
            self.ch[i].baseline_subtraction(samples=self.ch[i].trigger_position-int(0.5/self.ch[i].sample_step))
            self.ch[i].ar_filter(tau=20) # 20 samples = 80us = fast component
            self.ch[i].get_max(ar=True, trig=True) # AR matched filter, maximum near trigger position
            self.ch[i].get_integral() # full integral (from trigger-10 samples to end)
            # Make cut on filtered amplitude->SPE, baseline rms->pre-trigger pulses, and total integral->post-trigger scintillation pulses
            cut = (np.array(self.ch[i].output['baseline_rms'])<2.5) & \
                (np.array(self.ch[i].output['amplitude_trig'])<self.a1max[i]) & \
                (np.array(self.ch[i].output['amplitude_trig'])>self.a1min[i]) & \
                (np.array(self.ch[i].output['integral'])<6*self.gain[i])
            # Store SPE average waveform and number of selected waveforms
            self.ch[i].output['n_spe_wfs'] = np.sum(cut)
            self.ch[i].output['avg_spe_wf'] = np.dot(self.ch[i].traces.T,cut)
            self.ch[i].output['time'] = self.ch[i].time
        self.clear()

    def process_scintillation_pulses(self, header=True, num_events=1e9, calib=""):
        """Obtain pulse information from scintillation data.

        Args:
            header (bool, optional): Whether wavedump file contains header. Defaults to True.
            num_events (int, optional): Number of events. Defaults to 1e9.
            calib (str, optional): Directory that contains SiPM calibration results. Defaults to "".
        """
        self.read_calibration(calib)
        for i in self.channels:
            self.ch[i].read_data(header=header, num_events=num_events)
            self.ch[i].baseline_subtraction(samples=self.ch[i].trigger_position-int(0.5/self.ch[i].sample_step))
            self.ch[i].get_integral(length_us=[0.5,5]) # 0.5us for Fprompt analysis
        self.get_total_pe()
        self.get_fprompt()
        self.clear()

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
    
        