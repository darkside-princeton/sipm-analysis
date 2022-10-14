import numpy as np
import glob
from scipy import signal

class SiPM():
    def __init__(self, id, pol, path, samples):
        # basic information
        self.path = path
        self.file = None
        self.id = id # channel number
        self.pol = pol # polarization
        self.sampling = 250000000 # in Hz
        self.sample_step = 1./float(self.sampling)*1e6 # in us
        self.header = [0]*6
        self.timestamp = []
        self.samples = samples #waveform length
        self.nevents = 0
        self.cumulative_nevents = 0 # in case it needs to read multiple files
        # waveforms
        self.traces = [] #raw
        self.filtered_traces = [] # band pass filtered
        self.ar_filtered_traces = []# ar filtered
        self.time = [] #time array
        self.avgwf = np.zeros(0) # overall average
        self.spe_avgwf = None # spe waveform
        self.baseline_samples = 0 # will be modified to trigger_position - 10 samples
        self.trigger_position = 0   
        # band pass filter
        self.filt_pars = None
        # histograms
        self.famp = []
        self.famp_hist = None
        self.famp_hist_bin = None
        self.integral_prompt = []
        self.integral_prompt_hist = None
        self.integral_prompt_hist_bin = None
        self.integral_short = []
        self.integral_short_hist = None
        self.integral_short_hist_bin = None
        self.integral_long = []
        self.integral_long_hist = None
        self.integral_long_hist_bin = None
        # things from histograms
        self.gain_integral = 0 # spe integral
        self.integral_peaks = None # integral peak range (mu+/-2sigma)
        self.gain_famp = 0 # spe filtered amp
        self.famp_peaks = None # filtered amp peak range (mu+/-2sigma)
        # fitted pulse shape parameters and errors
        self.a1 = None
        self.tau1 = None
        self.a2 = None
        self.tau2 = None
        self.tau_singlet = None
        self.tau_triplet = None
        # after-pulse
        self.ap_charge = [] # list of list of long integral for each famp pe peak
        self.ap_charge_hist = [] # list of long integral histogram for each famp pe peak
        self.ap_charge_hist_bin = []
        
    
    def read_data(self, header=True, simple=False, verbose=False):
        self.file = glob.glob(self.path+"wave{}.dat".format(self.id))[0]
        file = open(self.file, 'rb')
        if header:
            for i in range(1000000):
                self.header = np.fromfile(file, dtype=np.dtype('I'), count=6)
                if len(self.header) == 0:
                    break
                self.samples = (self.header[0] - 24) // 2
                self.timestamp.append(self.header[-1])
                trace = np.fromfile(file, dtype=np.dtype('<H'), count=self.samples)
                self.traces.append(trace)
        if not header:
            self.traces = np.fromfile(file, dtype=np.dtype('<H'), count=-1)
        file.close()
        self.traces = np.array(self.traces)
        self.traces = self.traces.reshape((-1,self.samples)).astype(float)
        self.time = np.arange(0,self.sample_step*self.samples,self.sample_step)
        if self.avgwf.shape[0]==0:
            self.avgwf = np.zeros(self.samples)
        if not simple:
            self.filtered_traces = np.zeros(np.shape(self.traces))
            self.ar_filtered_traces = np.zeros(np.shape(self.traces))
            self.trigger_position = np.argmax(self.pol*np.mean(self.traces,axis=0))
            self.baseline_samples = self.trigger_position-50
            self.nevents = np.shape(self.traces)[0]
            self.cumulative_nevents += self.nevents
            if verbose:
                print('WAVEFORM LENGTH = {} SAMPLES'.format(self.samples))
                print('TRIGGER POSITION = SAMPLE {}'.format(self.trigger_position))
                print('CUMULATIVE WAVEFORMS = {}'.format(self.cumulative_nevents))

    def get_waveforms(self, event_id=[], header=True):
        if self.traces==[]:
            self.read_data(header=header, spe=True)
            self.baseline_subtraction()
        waveforms = []
        for event_id_ in event_id:
            waveforms.append(self.traces[event_id_,:])
        self.clear()
        return waveforms

    def baseline_subtraction(self):
        for ii,x in enumerate(self.traces):
            baseline = np.mean(self.traces[ii][:self.baseline_samples])
            self.traces[ii] -= baseline
            self.traces[ii] *= self.pol

    def bandpass_filter(self, low, high, order=3):
        if not self.filt_pars:
            b, a = signal.butter(order, [low,high], analog=False, fs=self.sampling, btype='band')
            self.filt_pars = [b,a]
        for ii,x in enumerate(self.traces):
            self.filtered_traces[ii] = signal.filtfilt(*self.filt_pars, x)

    def ar_filter(self, tau):
        wf_filt = np.zeros(self.traces.transpose().shape)
        for i,raw in enumerate(list(reversed(self.traces.transpose()))):
            if i>0:
                wf_filt[i] = raw + wf_filt[i-1]*np.exp(-1/tau)
            else:
                wf_filt[i] = raw
        self.ar_filtered_traces = np.array(list(reversed(wf_filt))).transpose()
            
    def get_max(self):
        self.peak = []
        self.peak_pos = []
        for ii,x in enumerate(self.filtered_traces):
            self.peak.append(np.max(x))
            self.peak_pos.append(np.argmax(x))

    def get_famp(self):
        self.famp = np.max(self.ar_filtered_traces, axis=1)

    def get_famp_hist(self, bin=[]):
        '''
        bin=[min, max, nbins]
        '''
        self.famp_hist, self.famp_hist_bin = np.histogram(self.famp, bins=bin[2], range=(bin[0],bin[1]))

    def get_integral(self, prompt=None, short=None, long=None):
        t0 = self.trigger_position
        if prompt!=None:
            tp = int(prompt/self.sample_step)
            for i,wf in enumerate(self.traces):
                self.integral_prompt.append(np.sum(wf[t0:t0+tp]))
        if short!=None:
            ts = int(short/self.sample_step)
            for i,wf in enumerate(self.traces):
                self.integral_short.append(np.sum(wf[t0:t0+ts]))
        if long!=None:
            tl = int(long/self.sample_step)
            for i,wf in enumerate(self.traces):
                self.integral_long.append(np.sum(wf[t0:t0+tl]))

    def get_integral_hist(self, prompt=None, short=None, long=None):
        '''
        prompt = [min_prompt, max_prompt, nbins_prompt]
        '''
        if prompt!=None:
            self.integral_prompt_hist, self.integral_prompt_hist_bin = np.histogram(self.integral_prompt, bins=prompt[2], range=(prompt[0],prompt[1]))
        if short!=None:
            self.integral_short_hist, self.integral_short_hist_bin = np.histogram(self.integral_short, bins=short[2], range=(short[0],short[1]))
        if long!=None:
            self.integral_long_hist, self.integral_long_hist_bin = np.histogram(self.integral_long, bins=long[2], range=(long[0],long[1]))

    def get_avgwf(self):
        self.avgwf = self.avgwf*(1-self.nevents/self.cumulative_nevents) + np.mean(self.traces,axis=0)*self.nevents/self.cumulative_nevents
        while self.avgwf[self.trigger_position]>1:
            self.trigger_position -= 1

    def get_spe_avgwf(self):
        if self.traces==[]:
            self.read_data(header=False,simple=True)
            self.baseline_subtraction()
        self.spe_avgwf = np.zeros(self.samples)
        count = 0
        for i,fa in enumerate(self.famp):
            if fa>self.famp_peaks[0][0] and fa<self.famp_peaks[0][1]:
                self.spe_avgwf *= count
                self.spe_avgwf += self.traces[i,:]
                count += 1
                self.spe_avgwf /= count

    def get_afterpulse_charge(self, xmin=-300, xmax=3000, nbins=500):
        self.ap_charge = []
        self.ap_charge_hist = []
        self.ap_charge_hist_bin = []
        for pe in range(len(self.famp_peaks)):
            ap_charge = []
            for i,fa in enumerate(self.famp):
                if fa>self.famp_peaks[pe][0] and fa<self.famp_peaks[pe][1]:
                    ap_charge.append(self.integral_long[i])
            ap_charge_hist, ap_charge_hist_bin = np.histogram(ap_charge, bins=nbins, range=(xmin,xmax))
            self.ap_charge.append(ap_charge)
            self.ap_charge_hist.append(ap_charge_hist)
            self.ap_charge_hist_bin.append(ap_charge_hist_bin)

    def set_calibration(self, gain_integral=0, integral_peaks=None, gain_famp=0, famp_peaks=None):
        if gain_integral!=0:
            self.gain_integral = gain_integral
        if integral_peaks!=None:
            self.integral_peaks = integral_peaks
        if gain_famp!=0:
            self.gain_famp = gain_famp
        if famp_peaks!=None:
            self.famp_peaks = famp_peaks

    def set_pulse_pars(self, a1=0, tau1=0, a2=0, tau2=0):
        self.a1 = a1
        self.tau1 = tau1
        self.a2 = a2
        self.tau2 = tau2
    
    def get_pulse_pars(self):
        return self.a1, self.tau1, self.a2, self.tau2

    def get_pulse_shape(self, t, a1=[0,0], tau1=[0,0], a2=[0,0], tau2=[0,0]):
        '''
        a1=[a1, a1_err]
        '''
        if a1==0:
            a1 = self.a1
        if tau1==0:
            tau1 = self.tau1
        if a2==0:
            a2 = self.a2
        if tau2==0:
            tau2 = self.tau2        
        return a1*np.exp(-(t-self.trigger_position*self.sample_step)/tau1) + a2*np.exp(-(t-self.trigger_position*self.sample_step)/tau2)        

    def get_scintillation(self, t, a_s, tau_s, a_t, tau_t):
        t_trg = self.trigger_position*self.sample_step
        s1 = self.a1*a_t*tau_t*self.tau1*(np.exp(-(t-t_trg)/tau_t)-np.exp(-(t-t_trg)/self.tau1))/(tau_t-self.tau1)
        s2 = self.a2*a_t*tau_t*self.tau2*(np.exp(-(t-t_trg)/tau_t)-np.exp(-(t-t_trg)/self.tau2))/(tau_t-self.tau2)
        s3 = self.a1*a_s*tau_s*self.tau1*(np.exp(-(t-t_trg)/tau_s)-np.exp(-(t-t_trg)/self.tau1))/(tau_s-self.tau1)
        s4 = self.a2*a_s*tau_s*self.tau2*(np.exp(-(t-t_trg)/tau_s)-np.exp(-(t-t_trg)/self.tau2))/(tau_s-self.tau2)
        return s1+s2+s3+s4

    def clear(self):
        self.traces = []
        self.filtered_traces = []
        self.ar_filtered_traces = []