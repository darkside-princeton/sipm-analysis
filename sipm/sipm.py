import numpy as np
import glob
from scipy import signal

class SiPM():
    def __init__(self, id, pol, path, samples):
        self.path = path
        self.id = id
        self.pol = pol
        self.file = glob.glob(self.path+"wave{}.dat".format(self.id))[0]
        self.sampling = 250000000 # in MHz
        self.sample_step = 1./float(self.sampling)*1e6 # in us
        self.traces = []
        self.filtered_traces = []
        self.time = []
        self.baseline_samples = 100 # will be modified to trigger_position - 10 samples
        self.filt_pars = None
        self.samples = samples
        self.header = [0]*6
        self.traces = []
        self.timestamp = []
        self.trigger_position = 0
        self.nevents = 0
        self.avgwf = []
        self.gain_integral = 0
        self.spe_integral = [0,0]
        self.gain_famp = 0
        self.spe_famp = [0,0]
    
    def read_data(self, header=True, spe=False):
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
        if not spe:
            self.filtered_traces = np.zeros(np.shape(self.traces))
            self.ar_filtered_traces = np.zeros(np.shape(self.traces))
            self.time = np.arange(0,self.sample_step*self.samples,self.sample_step)
            self.trigger_position = np.argmax(self.pol*np.mean(self.traces,axis=0))
            self.baseline_samples = self.trigger_position-100
            self.nevents = np.shape(self.traces)[0]
            print('WAVEFORM LENGTH = {} SAMPLES'.format(self.samples))
            print('TRIGGER POSITION = SAMPLE {}'.format(self.trigger_position))
            print('NUMBER OF WAVEFORMS = {}'.format(self.nevents))

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

    def get_ar_filtered_amp(self):
        self.famp = np.max(self.ar_filtered_traces, axis=1)

    def get_famp_hist(self, min=0, max=1e3, nbins=1000):
        self.famp_hist, self.famp_hist_bin = np.histogram(self.famp, bins=nbins, range=(min,max))

    def get_integral(self, length=5):
        '''
        default integral length = 5 us
        '''
        self.integral = []
        for ii,x in enumerate(self.traces):
            tmax = int(min(self.baseline_samples+length/self.sample_step, self.samples))
            self.integral.append(np.sum(x[self.baseline_samples:tmax]))

    def get_integral_hist(self, min=0, max=5e3, nbins=1000):
        self.integral_hist, self.integral_hist_bin = np.histogram(self.integral, bins=nbins, range=(min,max))

    def get_avgwf(self):
        self.avgwf = np.mean(self.traces,axis=0)

    def get_spe_avgwf(self):
        self.read_data(header=False,spe=True)
        self.baseline_subtraction()
        self.spe_avgwf = np.zeros(self.samples)
        count = 0
        for i,intg in enumerate(self.integral):
            if intg<self.spe_integral[1] and intg>self.spe_integral[0]:
                self.spe_avgwf *= count
                self.spe_avgwf += self.traces[i,:]
                count += 1
                self.spe_avgwf /= count
        self.clear()

    def get_afterpulse_charge(self, xmin=0, xmax=1e3, nbins=100):
        self.read_data(header=False,spe=True)
        self.baseline_subtraction()
        self.afterpulse_charge = []
        for i,fa in enumerate(self.famp):
            if fa<self.spe_famp[1] and fa>self.spe_famp[0]:
                self.afterpulse_charge.append(np.sum(self.traces[i,self.baseline_samples:self.baseline_samples+1250]))
        self.ap_charge_hist, self.ap_charge_hist_bin = np.histogram(self.afterpulse_charge, bins=nbins, range=(xmin,xmax))
        self.clear()

    def set_calibration(self, gain_integral=0, spe_integral=[0,0], gain_famp=0, spe_famp=[0,0]):
        if gain_integral!=0:
            self.gain_integral = gain_integral
        if spe_integral!=[0,0]:
            self.spe_integral = spe_integral
        if gain_famp!=0:
            self.gain_famp = gain_famp
        if spe_famp!=[0,0]:
            self.spe_famp = spe_famp

    def set_pulse_shape(self, a1=0, tau1=0, a2=0, tau2=0):
        self.a1 = a1
        self.tau1 = tau1
        self.a2 = a2
        self.tau2 = tau2

    def get_pulse_shape(self, t, a1=0, tau1=0, a2=0, tau2=0):
        if a1==0:
            a1 = self.a1
        if tau1==0:
            tau1 = self.tau1
        if a2==0:
            a2 = self.a2
        if tau2==0:
            tau2 = self.tau2        
        return a1*np.exp(-(t-self.trigger_position*self.sample_step)/tau1) + a2*np.exp(-(t-self.trigger_position*self.sample_step)/tau2)        
    
    def clear(self):
        self.traces = []
        self.filtered_traces = []
        self.ar_filtered_traces = []