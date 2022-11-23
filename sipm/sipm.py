import numpy as np
import glob
from scipy import signal
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from scipy.special import gamma, erf
from scipy.fft import fft, ifft
import ROOT
from array import array

def gauss(x,N,mu,sigma):
    return N*np.exp(-(x-mu)**2/(2*sigma**2))/np.sqrt(2*np.pi)/sigma

def compound_poisson(x,n,mu,p):
    k = [int(x_+0.5) for x_ in x]
    ans = []
    for k_ in k:
        if k_==0:
            ans.append(n*np.exp(-mu))
        else:
            ans_ = 0
            for i in range(1,k_+1):
                ans_ += gamma(k_+1)*gamma(k_)/gamma(i+1)/gamma(i)/gamma(k_-i+1)*(mu*(1-p))**i*p**(k_-i)
            ans.append(n*ans_*np.exp(-mu)/gamma(k_+1))
    return ans

def pulse_jitter(t, a, tau, sigma, t0):
    return a*np.exp(sigma**2/2/tau**2)*np.exp(-(t-t0)/tau)*(1+erf((t-t0-sigma**2/tau)/sigma/np.sqrt(2)))/2

class SiPM():
    def __init__(self, id, pol, path, root_file_name=None):
        # basic information
        self.path = path
        self.file = None
        self.id = id # channel number
        self.pol = pol # polarization
        self.sampling = 250000000 # in Hz
        self.sample_step = 1./float(self.sampling)*1e6 # in us
        self.header = [0]*6
        self.timestamp = []
        self.samples = 0 #waveform length
        self.nevents = 0
        self.acquisition_time = 0 # in seconds
        self.cumulative_nevents = 0 # in case it needs to read multiple files
        self.cumulative_time = 0 # in seconds
        # waveforms
        self.traces = [] #raw
        self.filtered_traces = [] # band pass filtered
        self.ar_filtered_traces = []# ar filtered
        self.deconv = [] #deconvolution
        self.time = [] #time array
        self.scint_sumwf = np.array([]) # scintillation summed waveform
        self.scint_avgwf = np.array([]) # scintillation averaged waveform
        self.scint_sumwf_count = 0
        self.spe_avgwf = np.array([]) # spe averaged waveform
        self.spe_sumwf = np.array([]) # spe summed waveform
        self.spe_sumwf_count = 0
        self.baseline_samples = 0 # will be modified to trigger_position - 10 samples
        self.trigger_position = 0
        # baseline
        self.baseline_avg = []
        self.baseline_med = []
        self.baseline_std = []
        self.acquisition_max = []
        self.acquisition_min = []
        # band pass filter
        self.filt_pars = None
        # histograms
        self.famp = []
        self.famp_hist = None
        self.famp_hist_bin = None
        self.famp_hist_fit = None
        self.integral_prompt = []
        self.integral_prompt_hist = None
        self.integral_prompt_hist_bin = None
        self.integral_short = []
        self.integral_short_hist = None
        self.integral_short_hist_bin = None
        self.integral_short_hist_fit = None
        self.integral_long = []
        self.integral_long_hist = None
        self.integral_long_hist_bin = None
        self.integral_long_hist_fit = None
        # spe gain
        self.q_a = None # SPE filtered amplitude [Q_A, Q_A err]
        self.q_peak = None # SPE 5us integral [Q_peak, Q_peak err]
        self.q_avg = None # SPE 5us integral + AP [Q_avg, Q_avg err]
        # fitted pulse shape parameters and errors
        self.a1 = []
        self.tau1 = []
        self.a2 = []
        self.tau2 = []
        self.tau_singlet = []
        self.tau_triplet = []
        # after-pulse
        self.ap_charge = [] # list of list of long integral for each famp pe peak
        self.ap_charge_hist = [] # list of long integral histogram for each famp pe peak
        self.ap_charge_hist_bin = []
    
    def read_data(self, header=True, simple=False, verbose=False):
        self.file = glob.glob(self.path+"wave{}.dat".format(self.id))[0]
        file = open(self.file, 'rb')
        self.acquisition_time = 0
        if header:
            for i in range(100000):
                self.header = np.fromfile(file, dtype=np.dtype('I'), count=6)
                if len(self.header) == 0:
                    break
                self.samples = (self.header[0] - 24) // 2
                self.timestamp.append(self.header[-1]%(2**31))
                if i>0:
                    self.acquisition_time += (self.timestamp[-1]-self.timestamp[-2])*8e-9
                    if self.timestamp[-1]<self.timestamp[-2]:
                        self.acquisition_time += 8e-9*2**31
                else:
                    self.acquisition_time += self.timestamp[-1]*8e-9
                trace = np.fromfile(file, dtype=np.dtype('<H'), count=self.samples)
                if np.shape(trace)[0]==self.samples:
                    self.traces.append(trace)
        if not header:
            self.traces = np.fromfile(file, dtype=np.dtype('<H'), count=-1)
        file.close()
        self.traces = np.array(self.traces)
        self.traces = self.traces.reshape((-1,self.samples)).astype(float)
        self.time = np.arange(0,self.sample_step*self.samples,self.sample_step)
        if self.scint_sumwf.shape[0]==0:
            self.scint_sumwf = np.zeros(self.samples)
        if self.spe_sumwf.shape[0]==0:
            self.spe_sumwf = np.zeros(self.samples)
        if not simple:
            self.filtered_traces = np.zeros(np.shape(self.traces))
            self.ar_filtered_traces = np.zeros(np.shape(self.traces))
            self.trigger_position = np.argmax(self.pol*np.mean(self.traces,axis=0))
            self.baseline_samples = self.trigger_position-int(0.5/self.sample_step)
            self.nevents = np.shape(self.traces)[0]
            self.cumulative_nevents += self.nevents
            self.cumulative_time += self.acquisition_time
            if verbose:
                print('WAVEFORM LENGTH = {} SAMPLES'.format(self.samples))
                print('TRIGGER POSITION = SAMPLE {}'.format(self.trigger_position))
                print('CUMULATIVE WAVEFORMS = {}'.format(self.cumulative_nevents))
                for k in range(10):
                    print(self.traces[k,:])

    def get_waveforms(self, event_id=[], header=True):
        if self.traces==[]:
            self.read_data(header=header, simple=True)
            self.baseline_subtraction()
        waveforms = []
        for event_id_ in event_id:
            waveforms.append(self.traces[event_id_,:])
        self.clear()
        return waveforms

    def baseline_subtraction(self, analysis=True):
        for ii,x in enumerate(self.traces):
            if analysis:
                self.baseline_avg.append(np.mean(self.traces[ii][:self.baseline_samples]))
                self.baseline_med.append(np.median(self.traces[ii][:self.baseline_samples]))
                self.baseline_std.append(np.std(self.traces[ii][:self.baseline_samples]))
                self.acquisition_max.append(np.max(self.traces[ii][:]))
                self.acquisition_min.append(np.min(self.traces[ii][:]))
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
            
    def deconvolution(self):
        h = lambda x: self.a1*np.exp(-x/self.tau1)+self.a2*np.exp(-x/self.tau2)
        for i,wf in enumerate(self.traces):
            ftilde = fft(wf)
            htilde = fft(h(self.time))
            self.deconv.append(ifft(ftilde/htilde).real)

    def get_max(self):
        self.peak = []
        self.peak_pos = []
        for ii,x in enumerate(self.filtered_traces):
            self.peak.append(np.max(x))
            self.peak_pos.append(np.argmax(x))

    def get_famp(self):
        self.famp = np.concatenate((self.famp, np.max(self.ar_filtered_traces[:,self.trigger_position-15:self.trigger_position+20], axis=1)))

    def get_famp_hist(self, bin=[]):
        '''
        bin=[min, max, nbins]
        '''
        self.famp_hist, self.famp_hist_bin = np.histogram(self.famp, bins=bin[2], range=(bin[0],bin[1]))

    def get_integral(self, prompt=None, short=None, long=None, deconv=False):
        t0 = self.baseline_samples
        if deconv:
            if prompt!=None:
                tp = self.trigger_position+int(prompt/self.sample_step)
                for i,wf in enumerate(self.deconv):
                    self.integral_prompt.append(np.sum(wf[t0:tp]))
            if short!=None:
                ts = self.trigger_position+int(short/self.sample_step)
                for i,wf in enumerate(self.deconv):
                    self.integral_short.append(np.sum(wf[t0:ts]))
            if long!=None:
                tl = self.trigger_position+int(long/self.sample_step)
                for i,wf in enumerate(self.deconv):
                    self.integral_long.append(np.sum(wf[t0:tl]))
        else:
            if prompt!=None:
                tp = self.trigger_position+int(prompt/self.sample_step)
                for i,wf in enumerate(self.traces):
                    self.integral_prompt.append(np.sum(wf[t0:tp]))
            if short!=None:
                ts = self.trigger_position+int(short/self.sample_step)
                for i,wf in enumerate(self.traces):
                    self.integral_short.append(np.sum(wf[t0:ts]))
            if long!=None:
                tl = self.trigger_position+int(long/self.sample_step)
                for i,wf in enumerate(self.traces):
                    self.integral_long.append(np.sum(wf[t0:tl]))

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

    def find_histo_peaks(self, hist, thre, prom, wid, dist):
        self.thre = thre
        if hist=='integral_short':
            self.peaks, self.pdict = find_peaks(self.integral_short_hist[thre:], prominence=prom, width=wid, distance=dist)
        elif hist=='integral_long':
            self.peaks, self.pdict = find_peaks(self.integral_long_hist[thre:], prominence=prom, width=wid, distance=dist)
        elif hist=='famp':
            self.peaks, self.pdict = find_peaks(self.famp_hist[thre:], prominence=prom, width=wid, distance=dist)
        else:
            print('no option called {}'.format(hist))
            print('available options: integral_short/integral_long/famp')

    def fit_histo_peaks(self, hist, bad=False):
        histo = None
        histo_bins = None
        if hist=='integral_short':
            histo = self.integral_short_hist
            histo_bins = self.integral_short_hist_bin
        elif hist=='integral_long':
            histo = self.integral_long_hist
            histo_bins = self.integral_long_hist_bin
        elif hist=='famp':
            histo = self.famp_hist
            histo_bins = self.famp_hist_bin
        else:
            print('no option called {}'.format(hist))
            print('available options: integral_short/integral_long/famp')
            return None
        bin_width = histo_bins[1]-histo_bins[0]
        gauss_fit = []
        min_bins = []
        max_bins = []
        for ip,peak_bin in enumerate(self.peaks):
            peak_bin = peak_bin + self.thre
            pe_width_bin = int(self.pdict['widths'][ip])
            pe_width_x = bin_width*pe_width_bin
            if not bad:
                min_bins.append(peak_bin-pe_width_bin)
            else:
                print('bad channel - change fit range')
                min_bins.append(peak_bin-int(1.5*pe_width_bin))
            max_bins.append(peak_bin+pe_width_bin)
            peak_x = histo_bins[0] + bin_width*peak_bin
            popt,pcov = curve_fit(gauss, histo_bins[min_bins[ip]:max_bins[ip]], histo[min_bins[ip]:max_bins[ip]], p0=[histo[peak_bin], peak_x, pe_width_x], sigma=np.sqrt(histo[min_bins[ip]:max_bins[ip]]), maxfev=10000)
            gauss_fit.append([ [popt[ipar], np.sqrt(pcov[ipar,ipar])] for ipar in range(3) ]) #[ [N,N_err], [mu,mu_err], [sigma,sigma_err] ]
        if hist=='integral_short':
            self.integral_short_hist_fit = gauss_fit
        elif hist=='integral_long':
            self.integral_long_hist_fit = gauss_fit
        elif hist=='famp':
            self.famp_hist_fit = gauss_fit
        return min_bins, max_bins

    def get_scint_sumwf(self, indices=[]):
        for i in indices:
            self.scint_sumwf += self.traces[i]
            self.scint_sumwf_count += 1

    def get_scint_avgwf(self):
        self.scint_avgwf = self.scint_sumwf/self.scint_sumwf_count

    def get_spe_sumwf(self, famp_range=(0,0)):
        for i,wf_filt in enumerate(self.ar_filtered_traces):
            fa = np.max(wf_filt[self.trigger_position-15:self.trigger_position+20])
            if fa>famp_range[0] and fa<famp_range[1] and  self.baseline_std[i]<2.5:
                self.spe_sumwf += self.traces[i,:]
                self.spe_sumwf_count += 1
                
    def get_spe_avgwf(self):
        self.spe_avgwf = self.spe_sumwf/self.spe_sumwf_count

    def get_afterpulse_charge(self, nsigma=3, bin=[-300, 6000, 500]):
        self.ap_charge = []
        self.ap_charge_hist = []
        self.ap_charge_hist_bin = []
        for pe in range(len(self.famp_hist_fit)):
            ap_charge = []
            if pe==0:
                min_fa = 0.5*self.famp_hist_fit[pe][1][0]
            else:
                min_fa = 0.5*(self.famp_hist_fit[pe][1][0]+self.famp_hist_fit[pe][1][0])
            if pe==len(self.famp_hist_fit)-1:
                max_fa = 1.5*self.famp_hist_fit[pe][1][0]-0.5*self.famp_hist_fit[pe-1][1][0]
            else:
                max_fa = 0.5*(self.famp_hist_fit[pe+1][1][0]+self.famp_hist_fit[pe][1][0])

            for i,fa in enumerate(self.famp):
                # if abs(fa-self.famp_hist_fit[pe][1][0]) < nsigma*self.famp_hist_fit[pe][2][0]:
                if fa<max_fa and fa>min_fa:
                    ap_charge.append(self.integral_long[i])
            ap_charge_hist, ap_charge_hist_bin = np.histogram(ap_charge, bins=bin[2], range=(bin[0],bin[1]))
            self.ap_charge.append(ap_charge)
            self.ap_charge_hist.append(ap_charge_hist)
            self.ap_charge_hist_bin.append(ap_charge_hist_bin)

    def set_spe_gain(self, q_a=None, q_peak=None, q_avg=None):
        if q_a!=None:
            self.q_a = q_a
        if q_peak!=None:
            self.q_peak = q_peak
        if q_avg!=None:
            self.q_avg = q_avg

    def set_correlated_noise(self, ap=0, ct=0):
        if ap!=0:
            self.ap = ap
        if ct!=0:
            self.ct = ct

    def set_pulse_pars(self, a1=0, tau1=0, a2=0, tau2=0):
        self.a1 = a1
        self.tau1 = tau1
        self.a2 = a2
        self.tau2 = tau2
    
    def get_pulse_pars(self):
        return self.a1, self.tau1, self.a2, self.tau2

    def get_pulse_shape(self, t, a1, tau1, a2, tau2, sigma, t0):
        return pulse_jitter(t, a1, tau1, sigma, t0) + pulse_jitter(t, a2, tau2, sigma, t0)
        
    def get_scintillation(self, t, a_s, tau_s, a_t, tau_t, sigma, t0):
        return pulse_jitter(t, a_s, tau_s, sigma, t0) + pulse_jitter(t, a_t, tau_t, sigma, t0)

    def clear(self):
        self.traces = []
        self.filtered_traces = []
        self.ar_filtered_traces = []
        self.deconv = []
    
    def clear_all(self):
        self.clear()
        # basic information
        self.timestamp = []
        self.nevents = 0
        self.acquisition_time = 0 # in seconds
        # waveforms
        self.traces = [] #raw
        self.filtered_traces = [] # band pass filtered
        self.ar_filtered_traces = []# ar filtered
        self.deconv = [] #deconvolution
        # baseline
        self.baseline_avg = []
        self.baseline_med = []
        self.baseline_std = []
        self.acquisition_max = []
        self.acquisition_min = []
        # histograms
        self.famp = []
        self.integral_prompt = []
        self.integral_short = []
        self.integral_long = []