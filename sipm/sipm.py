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
        self.time = []
        self.baseline_samples = 100
        self.filt_pars = None
        self.samples = samples
        self.header = [0]*6
        self.traces = []
        self.timestamp = []
    
    def read_data(self, header=True):
        file = open(self.file, 'rb')
        for i in range(1000000):
            if header:
                self.header = np.fromfile(file, dtype=np.dtype('I'), count=6)
                if len(self.header) == 0:
                    break
                self.samples = (self.header[0] - 24) // 2
                self.timestamp.append(self.header[-1])
                trace = np.fromfile(file, dtype=np.dtype('<H'), count=self.samples)
            else:
                trace = np.fromfile(file, dtype=np.dtype('<H'), count=-1)
            self.traces.append(trace)
        file.close()
        self.traces = np.array(self.traces)
        self.traces = self.traces.reshape((-1,self.samples)).astype(float)
        self.time = np.arange(0,self.sample_step*self.samples,self.sample_step)

    def baseline_subtraction(self):
        for ii,x in enumerate(self.traces):
            baseline = np.mean(self.traces[ii][:500])
            self.traces[ii] -= baseline
            self.traces[ii] *= self.pol

    def bandpass_filter(self, low, high, order=3):
        if not self.filt_pars:
            b, a = signal.butter(order, [low,high], analog=False, fs=self.sampling, btype='band')
            self.filt_pars = [b,a]
        for ii,x in enumerate(self.traces):
            self.traces[ii] = signal.filtfilt(*self.filt_pars, x)
    
    def get_max(self):
        self.peak = []
        self.peak_pos = []
        for ii,x in enumerate(self.traces):
            self.peak.append(np.max(x))
            self.peak_pos.append(np.argmax(x))

    def get_integral(self):
        self.integral = []
        for ii,x in enumerate(self.traces):
            self.integral.append(np.sum(x[self.peak_pos[ii]-50:self.peak_pos[ii]+500]))
            # self.integral.append(np.sum(x[self.peak_pos[ii]-50:]))
            # self.integral.append(np.sum(x))
    
    def clear(self):
        self.traces = []