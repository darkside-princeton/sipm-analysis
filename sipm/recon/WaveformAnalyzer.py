import matplotlib.pyplot as plt
import numpy as np
import glob
import scipy
from BaselineRemoval import BaselineRemoval
import sipm.util.functions as func

class WaveformAnalyzer():
    def __init__(self, id, pol, path, samples):
        """Class that analyzes the waveform data to obtain higher-level information.

        Args:
            id (int): Channel number
            pol (int): Polarity
            path (str): Data folder path
            samples (int): Length of acquisition window
        """
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
        self.ar_filtered_traces = []
        self.timestamp = []
        self.trigger_position = 0
        self.nevents = 0
        self.output = {}
    
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
            print(f)
            file = open(f, 'rb')
            if header:
                for i in range(int(num_events)):
                    self.header = np.fromfile(file, dtype=np.dtype('I'), count=6)
                    if len(self.header) == 0 or i==num_events:
                        break
                    self.samples = (self.header[0] - 24) // 2
                    self.timestamp.append(self.header[-1])
                    trace = np.fromfile(file, dtype=np.dtype('<H'), count=self.samples)
                    self.traces.append(trace)
            else:
                self.traces = np.fromfile(file, dtype=np.dtype('<H'), count=-1)
                cutoff = -len(self.traces)%self.samples
                self.traces = self.traces[:cutoff].reshape((-1,self.samples)).astype(float)
            file.close()
        self.traces = np.array(self.traces).astype(float)
        self.nevents = self.traces.shape[0]
        print(f'{self.nevents} events')     
        self.time = np.arange(0,self.sample_step*self.samples,self.sample_step)
        self.trigger_position = np.argmax(self.pol*np.mean(self.traces, axis=0))

    def baseline_subtraction(self, samples=500):
        """Computes baseline mean and rms and subtracts baseline mean from raw waveform

        Args:
            samples (int, optional): Baseline evaluation window starting from the beginning of acquisition. Defaults to 500.
        """
        self.output['baseline_mean'] = []
        self.output['baseline_rms'] = []
        for ii,x in enumerate(self.traces):
            baseline = np.mean(self.traces[ii][:samples])
            self.output['baseline_mean'].append(baseline)
            self.output['baseline_rms'].append(np.std(self.traces[ii][:samples]))
            self.traces[ii] -= baseline
            self.traces[ii] *= self.pol

    def bandpass_filter(self, low, high, order=3, keep=False):
        if not self.filt_pars:
            b, a = scipy.signal.butter(order, [low,high], analog=False, fs=self.sampling, btype='band')
            self.filt_pars = [b,a]
        if keep:
            self.traces_orig = self.traces.copy()
        for ii,x in enumerate(self.traces):
            self.traces[ii] = scipy.signal.filtfilt(*self.filt_pars, x)
    
    def ar_filter(self, tau):
        """Matched filter using auto-regressive algorithm

        Args:
            tau (float): Time constant in units of sample.
        """
        wf_filt = np.zeros(self.traces.transpose().shape)
        for i,raw in enumerate(list(reversed(self.traces.transpose()))):
            if i>0:
                wf_filt[i] = raw + wf_filt[i-1]*np.exp(-1/tau)
            else:
                wf_filt[i] = raw
        self.ar_filtered_traces = np.array(list(reversed(wf_filt))).transpose()
    
    def get_waveforms(self, ev=[], ar_filter=True):
        """Return several baseline-subtracted waveforms

        Args:
            ev (list, optional): A list of event ids. Defaults to [].
            ar_filter (bool, optional): Whether to return ar-filtered waveforms. Defaults to True.

        Returns:
            list, list: baseline-subtracted waveforms, ar-filtered waveforms
        """
        if self.traces==[]:
            self.read_data()
            self.baseline_subtraction()
            if ar_filter:
                self.ar_filter(tau=20)
        wfs = []
        arwfs = []
        for ev_ in ev:
            wfs.append(self.traces[ev_,:])
            arwfs.append(self.ar_filtered_traces[ev_,:])
        self.clear()
        return wfs, arwfs

    def get_famp_hist(self, bin):
        self.famp_hist, self.famp_hist_bin = np.histogram(self.famp, bins=bin[2], range=(bin[0],bin[1]))

    def find_histo_peaks(self, thre, prom, wid, dist):
        self.peaks, self.pdict = scipy.signal.find_peaks(self.famp_hist[thre:], prominence=prom, width=wid, distance=dist)
        
    def get_charge_histo_pe(self, bin):
        self.ap_charge = []
        self.ap_charge_hist = []
        self.ap_charge_hist_bin = []
        for pe in range(len(self.pe_cuts)-1):
            ap_charge = []
            if pe==0:
                min_fa = self.pe_cuts[0]
            else:
                min_fa = self.pe_cuts[pe]
            if pe==len(self.pe_cuts)-2:
                max_fa = self.pe_cuts[-1]
            else:
                max_fa = self.pe_cuts[pe+1]

            for i,fa in enumerate(self.famp):
                if fa<max_fa and fa>min_fa:
                    ap_charge.append(self.integral_long[i])
            ap_charge_hist, ap_charge_hist_bin = np.histogram(ap_charge, bins=bin[2], range=(bin[0],bin[1]))
            self.ap_charge.append(ap_charge)
            self.ap_charge_hist.append(ap_charge_hist)
            self.ap_charge_hist_bin.append(ap_charge_hist_bin)

    def get_max(self, traces=None, ar=False, trig=False):
        """Evaluate amplitude.

        Args:
            traces (_type_, optional): _description_. Defaults to None.
            ar (bool, optional): Whether to use AR filter. Defaults to False.
            trig (bool, optional): Whether to confine the evaluation within 20 samples from the trigger position (average maximum position). Defaults to False.
        """
        if traces == None:
            if ar:
                traces = self.ar_filtered_traces
            else:
                traces = self.traces
        if trig:
            self.output['amplitude_trig'] = []
            self.output['peakpos_trig'] = []
            for ii,x in enumerate(traces):
                self.output['amplitude_trig'].append(np.max(x[self.trigger_position-20:self.trigger_position+20]))
                self.output['peakpos_trig'].append(self.trigger_position-20+np.argmax(x[self.trigger_position-20:self.trigger_position+20]))
        else:
            self.output['amplitude'] = []
            self.output['peakpos'] = []
            for ii,x in enumerate(traces):
                self.output['amplitude'].append(np.max(x))
                self.output['peakpos'].append(np.argmax(x))

    def get_integral(self, traces=None, length_us=[]):
        """Evaluate charge integral.

        Args:
            traces (_type_, optional): _description_. Defaults to None.
            length_us (list, optional): A list of integration time windows. Defaults to [].
        """
        if traces == None:
            traces = self.traces
        if len(length_us)==0:
            self.output['integral'] = []
            for ii,x in enumerate(traces):
                self.output['integral'].append(np.sum(x[self.trigger_position-10:]))
                # self.integral.append(np.sum(x[self.peak_pos[ii]-50:]))
                # self.integral.append(np.sum(x))
        else:
            for length in length_us:
                length_digits = str(int(length*100))
                while len(length_digits)<3:
                    length_digits = '0'+length_digits
                name = f'integral_{length_digits[:-2]}p{length_digits[-2:]}us'
                self.output[name] = []
                for ii,x in enumerate(traces):
                    self.output[name].append(np.sum(x[self.trigger_position-10:self.trigger_position+int(length/self.sample_step)]))

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
            y = scipy.fft.fft(x)
            y = 2.0/self.samples * np.abs(y[1:self.samples//2])
            self.fft.append(y)
            
    def clear(self):
        self.traces = []
        self.ar_filtered_traces = []

    def calibrate(self,vals,width=6,prominence=2,verbose=-1,fitrange=10):
        self.fit_p = []
        self.fit_c = []

        h,hx = np.histogram(vals, bins=np.arange(0,np.max(vals),.1))
        self.pp = []
        while len(self.pp)<=3:
            print('Peak search with width={} ...'.format(width))
            self.pp,self.pdict = scipy.signal.find_peaks(h, prominence=prominence, width=width)
            width -= 1

        if verbose > 0:
            print("Found {} peaks".format(len(self.pp)))

        for x in self.pp: 
            fit_x = hx[:-1][x-fitrange:x+fitrange]
            fit_y = h[x-fitrange:x+fitrange]
            popt,pcov = scipy.optimize.curve_fit(func.gauss, fit_x, fit_y, p0=[h[x], hx[:-1][x], 2], maxfev=10000)
            self.fit_p.append(popt)
            self.fit_c.append(pcov)

        # quick estimate of gain
        gain = np.median(np.diff(np.array(self.fit_p)[:,1]))
        # get peak number
        peak_num = np.round(np.array(self.fit_p)[:,1]/gain)

        self.calib,self.calib_err = scipy.optimize.curve_fit(func.line, peak_num, np.array(self.fit_p)[:,1])
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
