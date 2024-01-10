import numpy as np
import sipm.recon.WaveformAnalyzer as wfa
import glob
import csv
from datetime import datetime
import pandas as pd

class WaveformDataset: 
    def __init__(self,  path, pol=-1, channels=[0,1,2,3], samples=3000):
        """Class that analyzes the waveform data. Contains a list of WaveformAnalyzer objects for all the channels in a SiPM array.

        Args:
            path (str): Data folder path
            pol (int, optional): Polarity. Defaults to -1.
            channels (list, optional): A list of channel numbers. Defaults to [0,1,2,3].
            samples (int, optional): Length of acquisition window. Defaults to 3000.
        """
        self.datadir = "/scratch/gpfs/GALBIATI/data/sipm/"
        self.path = f"{self.datadir}/{path.replace(self.datadir, '')}"
        self.channels = channels
        self.samples = samples
        self.pol = pol
        self.pos = ''
        self.volt = 0
        self.output = {}
        self.read_timestamp()
        self.ch = self.InitializeChannels()
        
    def read_timestamp(self):
        """Read the time information for the data file.
        """
        format_str = '%Y/%m/%d %H:%M:%S'
        try:
            with open(f'{self.path}time.txt', 'r') as f:
                lines = f.readlines()
            self.start_datetime = datetime.strptime(lines[0][:-1],format_str)
            self.output['start_datetime'] = np.array([self.start_datetime.year,
                                                    self.start_datetime.month, 
                                                    self.start_datetime.day,
                                                    self.start_datetime.hour,
                                                    self.start_datetime.minute,
                                                    self.start_datetime.second]).astype(float)
            self.end_datetime = datetime.strptime(lines[1][:-1],format_str)
            self.output['end_datetime'] = np.array([self.end_datetime.year,
                                                    self.end_datetime.month,
                                                    self.end_datetime.day,
                                                    self.end_datetime.hour,
                                                    self.end_datetime.minute,
                                                    self.end_datetime.second]).astype(float)
            self.duration_seconds = float(lines[2][:-2])
            self.output['duration_seconds'] = np.array([self.duration_seconds]).astype(float)
        except:
            print('No time.txt file')
 
    def InitializeChannels(self):
        channels = []
        for i in self.channels:
            new_channel = wfa.WaveformAnalyzer(id=i, pol=self.pol, path=self.path, samples=self.samples)
            channels.append(new_channel)
        # Get position and voltage
        if self.path.find('pos_')!=-1:
            id_pos = [self.path.find('pos_')+4, self.path.find('_',self.path.find('pos_')+4)]
            self.pos = self.path[id_pos[0]:id_pos[1]]
        if  self.path.find('intensity_')!=-1:
            id_intn = [self.path.find('intensity_')+10, self.path.find('/',self.path.find('intensity_')+10)]
            self.intensity = self.path[id_intn[0]:id_intn[1]]
        id_volt = [self.path.find('volt_')+5, self.path.find('_',self.path.find('volt_')+5)]
        self.volt = int(self.path[id_volt[0]:id_volt[1]])
        return np.array(channels)

    def read_calibration_h5(self, filename):
        """Read calibration result HDF5 file. See the member function 'SipmCalibration::write_to_h5()' in SipmCalibration.py for an example to generate such a file.

        Args:
            filename (_type_): Full path of the file.
        """
        self.calib_df = pd.read_hdf(filename, key=f'/{self.pos}/{self.volt}V')
        self.calib_df['cn_corrected_gain'] = self.calib_df['Qpeak']*(1+self.calib_df['Qap'])/(1-self.calib_df['DiCT']) # effective SPE gain corrected for correlated noises (DiCT and afterpulsing)

    def get_total_pe(self):
        self.output['total_pe'] = np.zeros(self.ch[0].nevents)
        for ch in self.channels:
            self.output['total_pe'] += np.array(self.ch[ch].output['integral_5p00us'])/self.calib_df['cn_corrected_gain'][ch]

    def get_fprompt(self, tprompt=[0.5], channels=np.arange(4)):
        integral_long = np.zeros(self.ch[0].nevents)
        integral_short = np.zeros(self.ch[0].nevents)
        channels_str = ''.join(channels.astype(str))
        for tp in tprompt:
            length_digits = str(int(tp*100))
            while len(length_digits)<3:
                length_digits = '0'+length_digits
            name = f'{length_digits[:-2]}p{length_digits[-2:]}'
            for ch in channels:
                integral_long += np.array(self.ch[ch].output['integral_5p00us'])
                integral_short += np.array(self.ch[ch].output[f'integral_{name}us'])
            self.output[f'fprompt_{name}us_{channels_str}'] = integral_short/integral_long
            integral_long = np.zeros(self.ch[0].nevents)
            integral_short = np.zeros(self.ch[0].nevents)

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
        self.read_calibration_h5(calib)
        for i in self.channels:
            self.ch[i].read_data(header=header, num_events=num_events)
            self.ch[i].baseline_subtraction(samples=self.ch[i].trigger_position-int(0.5/self.ch[i].sample_step))
            self.ch[i].ar_filter(tau=20) # 20 samples = 80us = fast component
            self.ch[i].get_max(ar=True, trig=True) # AR matched filter, maximum near trigger position
            self.ch[i].get_integral() # full integral (from trigger-10 samples to end)
            # Make cut on filtered amplitude->SPE, baseline rms->pre-trigger pulses, and total integral->post-trigger scintillation pulses
            cut = (np.array(self.ch[i].output['baseline_rms'])<self.calib_df['bsl_rms'][i]) & \
                (np.array(self.ch[i].output['amplitude_trig'])<self.calib_df['A1max'][i]) & \
                (np.array(self.ch[i].output['amplitude_trig'])>self.calib_df['A1min'][i]) & \
                (np.array(self.ch[i].output['integral'])<6*self.calib_df['cn_corrected_gain'][i])
            # Store SPE average waveform and number of selected waveforms
            self.ch[i].output['n_spe_wfs'] = np.sum(cut)
            self.ch[i].output['avg_spe_wf'] = np.dot(self.ch[i].traces.T,cut)/self.ch[i].output['n_spe_wfs']
            self.ch[i].output['time'] = self.ch[i].time
            # Clean up unnecessary variables
            self.ch[i].output.pop('baseline_mean')
            self.ch[i].output.pop('baseline_rms')
            self.ch[i].output.pop('amplitude_trig')
            self.ch[i].output.pop('peakpos_trig')
            self.ch[i].output.pop('integral')
        self.clear()

    def process_scintillation_pulses(self, header=True, num_events=1e9, calib=""):
        """Obtain pulse information from scintillation data.

        Args:
            header (bool, optional): Whether wavedump file contains header. Defaults to True.
            num_events (int, optional): Number of events. Defaults to 1e9.
            calib (str, optional): Directory that contains SiPM calibration results. Defaults to "".
        """
        self.read_calibration_h5(calib)
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
        self.read_calibration_h5(calib)
        for i in self.channels:
            self.ch[i].read_data(header=header, num_events=num_events)
            self.ch[i].baseline_subtraction(samples=self.ch[i].trigger_position-int(0.5/self.ch[i].sample_step))
            self.ch[i].get_integral(length_us=[0.5,5]) # 0.5us for Fprompt analysis
        self.get_total_pe()
        self.get_fprompt()
        # Make cut on total pe, fprompt, baseline rms of all the channels
        cut = (np.array(self.output['total_pe'])<pe[1]) & \
            (np.array(self.output['total_pe'])>pe[0]) & \
            (np.array(self.output['fprompt'])<fprompt[1]) & \
            (np.array(self.output['fprompt'])>fprompt[0])
        for i in self.channels:
            cut = cut & (np.array(self.ch[i].output['baseline_rms'])<self.calib_df['bsl_rms'][i])
        # Store average LAr scintillation waveform and number of selected waveforms
        for i in self.channels:
            self.ch[i].output['n_scint_wfs'] = np.sum(cut)
            self.ch[i].output['avg_scint_wf'] = np.dot(self.ch[i].traces.T,cut)/self.ch[i].output['n_scint_wfs']
            self.ch[i].output['time'] = self.ch[i].time
        # Clean up unnecessary variables
        self.output.pop('total_pe')
        self.output.pop('fprompt')
        for i in self.channels:
            self.ch[i].output.pop('baseline_mean')
            self.ch[i].output.pop('baseline_rms')
            self.ch[i].output.pop('integral_0p50us')
            self.ch[i].output.pop('integral_5p00us')
        self.clear()

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
    
        