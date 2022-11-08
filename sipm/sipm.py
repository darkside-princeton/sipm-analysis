import numpy as np
import glob
import matplotlib.pyplot as plt

from scipy import signal
from scipy.fft import fft
from scipy.optimize import curve_fit

import sipm.functions as func

from BaselineRemoval import BaselineRemoval


class SiPM():
    def __init__(self, id, pol, path, samples):
        self.path = path
        self.id = id
        self.pol = pol
        self.files = glob.glob(self.path+"wave{}.dat".format(self.id))
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
    
    def read_data(self, header=True, num_events=1e9):
        """Reads data from the binary wavedump file storing the waveforms.

        Parameters
        ----------
        header : bool, optional
            Whether or not the wavedump file includes an header for each waveform which includes information about the length of the waveform and the unique timestamp (default is True)
        num_events : int, optional
            Number of waveforms to read from each wavedump file (default is all (1e9))
        """
        for f in sorted(self.files):
            file = open(f, 'rb')
            if header:
                for i in range(1000000):
                    self.header = np.fromfile(file, dtype=np.dtype('I'), count=6)
                    if len(self.header) == 0 or i>num_events:
                        break
                    self.samples = (self.header[0] - 24) // 2
                    self.timestamp.append(self.header[-1])
                    trace = np.fromfile(file, dtype=np.dtype('<H'), count=self.samples)
                    self.traces.append(trace)
            else:
                self.traces = np.fromfile(file, dtype=np.dtype('<H'), count=-1)
            file.close()

        self.traces = np.array(self.traces)
        self.traces = self.traces.reshape((-1,self.samples)).astype(float)
        self.time = np.arange(0,self.sample_step*self.samples,self.sample_step)

    def baseline_subtraction(self, samples=500):
        for ii,x in enumerate(self.traces):
            baseline = np.mean(self.traces[ii][:samples])
            self.traces[ii] -= baseline
            self.traces[ii] *= self.pol

    def bandpass_filter(self, low, high, order=3, keep=False):
        if not self.filt_pars:
            b, a = signal.butter(order, [low,high], analog=False, fs=self.sampling, btype='band')
            self.filt_pars = [b,a]
        if keep:
            self.traces_orig = self.traces.copy()
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
            self.integral.append(np.sum(x[1500:]))
            # self.integral.append(np.sum(x[self.peak_pos[ii]-50:]))
            # self.integral.append(np.sum(x))

    def rolling_baseline(self):
        return 0
    
    def get_rolling_integral(self):
        self.rolling_integral = []
        for x in self.traces:
            trace_corr = BaselineRemoval(x)
            self.rolling_integral.append(np.sum(trace_corr.IModPoly(10)))
        self.rolling_integral = np.array(self.rolling_integral)

    def get_fft(self):
        self.fft = []
        for ii,x in enumerate(self.traces):
            y = fft(x)
            y = 2.0/self.samples * np.abs(y[1:self.samples//2])
            self.fft.append(y)
            
    def clear(self):
        self.traces = []

    def calibrate(self,vals,width=6,prominence=2,verbose=-1,fitrange=10):
        self.fit_p = []
        self.fit_c = []

        h,hx = np.histogram(vals, bins=np.arange(0,np.max(vals),.1))
        self.pp = []
        while len(self.pp)<=3:
            print('Peak search with width={} ...'.format(width))
            self.pp,self.pdict = signal.find_peaks(h, prominence=prominence, width=width)
            width -= 1

        if verbose > 0:
            print("Found {} peaks".format(len(self.pp)))

        for x in self.pp: 
            fit_x = hx[:-1][x-fitrange:x+fitrange]
            fit_y = h[x-fitrange:x+fitrange]
            popt,pcov = curve_fit(func.gauss, fit_x, fit_y, p0=[h[x], hx[:-1][x], 2], maxfev=10000)
            self.fit_p.append(popt)
            self.fit_c.append(pcov)

        # quick estimate of gain
        gain = np.median(np.diff(np.array(self.fit_p)[:,1]))
        # get peak number
        peak_num = np.round(np.array(self.fit_p)[:,1]/gain)

        self.calib,self.calib_err = curve_fit(func.line, peak_num, np.array(self.fit_p)[:,1])
        self.calib_err = np.sqrt(np.diag(self.calib_err))

    def get_breakdown(self):
        return 0 

        
    def plot_calibration(self,vals):
        fig, ax = plt.subplots(figsize=(9,3), ncols=2,nrows=1)
        h,hx = np.histogram(vals, bins=np.arange(0,np.max(vals),.1))
        ax[0].step(hx[:-1],h, where='post')
        ax[0].scatter(hx[self.pp], h[self.pp], color='r', s=3, zorder=10)
        ax[0].set_xlabel('Amplitude [mV]')
        ax[0].set_ylabel('Counts')
        ax[0].set_yscale('log')
        ax[0].set_xlim(0,np.max(hx[self.pp])+20)
        ax[0].set_ylim(1e0,1e4)

        # quick estimate of gain
        gain = np.median(np.diff(np.array(self.fit_p)[:,1]))
        # get peak number
        peak_num = np.round(np.array(self.fit_p)[:,1]/gain)

        xfit = np.linspace(0,20,100)

        ax[1].grid()
        ax[1].set_ylim(0,50)
        ax[1].scatter(peak_num, np.array(self.fit_p)[:,1], marker='o', s=20, zorder=10)
        ax[1].plot(xfit, func.line(xfit, *self.calib), ls='--', label=r'$m=({:.3f}\pm {:.3f})$/p.e.'.format(self.calib[0], self.calib_err[0]))
        ax[1].legend(loc='lower right')
        plt.show()
