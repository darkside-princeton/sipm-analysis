import numpy as np
import glob

from scipy.signal import butter, filtfilt
from scipy.special import erfc
from scipy.optimize import curve_fit
from scipy.signal import find_peaks

from WaveformAnalysis import Dataset

class SiPM(Dataset.Dataset):
    def __init__(self, Path, Selection='*'):
        self.Path = Path
        self.Selection = Selection
        self.Files = glob.glob(self.Path+self.Selection)
        self.fit_parameters = []
        self.fit_covariance = []
        self.peak_pos = []
        self.peak_height = []
        self.integral = {}
        self.max = []
        self.deconv_filter = None
        self.shaping_time=[]
        self.filter_coefficients=[]
        self.sampling_freq=None #Sanpling rate for the waveforms in Hz

    def get_filtered_waveform(self,amp):
        filtered_waveform=[]
        for filt_ceoffs in self.filter_coefficients:
            b,a=filt_ceoffs[0],filt_ceoffs[1]
            filtered_waveform.append(filtfilt(b,a,amp))
        return filtered_waveform



    def get_averaged_waveform(self, time, amp, avg=4):
        avg_time = time[:-avg+1]
        avg_amp = np.convolve(amp,np.ones(avg),'valid')/avg
        return avg_time, avg_amp

    def run_fit(self, time, amp):
        max_val = np.max(Amp[cut])
        max_pos_cut = np.where(Amp[cut] == np.max(Amp[cut]))[0][0]
        max_pos = D0.Ch[1].Time[cut][max_pos_cut]



    def setup_butter_filter(self, order=3):
        '''
        calculates filter coefficients for different shaping times
        '''
    
        nyq = 0.5 * self.sampling_freq #this is the Nyquist frequency
        
        self.filter_coefficients=[]
        for s_time in self.shaping_time:

            b, a = butter(order, [(1/s_time)/nyq], btype='low', analog=False)
            self.filter_coefficients.append([b,a])

            
    
    def fit_peaks(self, time, data):
        peaks,pdict = find_peaks(data, height=35, width=20, distance=50)
        self.peak_pos.append(peaks)
        self.peak_height.append(pdict['peak_heights'])
        
        for ii,p in enumerate(peaks): 
            pp = time[p]
            ph = pdict['peak_heights'][ii]
            fit = np.where((time>time[p]-50) & (time<time[p]+50))
            try:
                popt, pcov = curve_fit(self.func, time[fit], data[fit], p0=[0.0, ph, 1.0, 40.0, pp], maxfev=10000)
                self.fit_parameters.append(popt)
                self.fit_covariance.append(pcov)
                self.max.append( np.max(self.func(time[fit], *popt)) )
            except:
                self.fit_parameters.append([0,0,0,0,0])
                self.fit_covariance.append([])
    
    def func(self,x,base,V0,sigma,tau,mu):
        return base + V0/2.0 * np.exp(0.5 * (sigma/tau)**2 - (x-mu)/tau) * erfc(1.0/np.sqrt(2.0) * (sigma/tau - (x-mu)/sigma))
    
    def gauss_conv(self, x, mu=0, sigma=0.1):
        x = x-mu
        return np.exp(-np.square((x-mu)/sigma)/2)/(sigma*np.sqrt(2*np.pi))



    
    def get_sampling(self):
        '''
        Sets up sampling time and frequency that gets called only once
        '''
        self.start = self.Ch[0].Time[0]
        self.length = self.Ch[0].Time[-1] - self.Ch[0].Time[0]
        # print(self.Ch[0].Time[1])
        self.sampling_freq=(1/(self.Ch[0].Time[1] - self.Ch[0].Time[0]))*self.Ch[0].TScale


    def get_convolution_filter(self):
        xdata = self.Ch[0].Time
        desire_x = np.arange(self.start,self.start+self.length, xdata[1]-xdata[0])
        desire_x = np.arange(-500,500, xdata[1]-xdata[0])
        desire_y = self.gauss_conv(desire_x)
        resp_x = np.arange(0,1000, xdata[1]-xdata[0])
        resp_y = self.func(resp_x, 0, 1, 1.83, 46.93, 0)
        resp_f = np.fft.fft(resp_y)
        desire_f = np.fft.fft(desire_y)
        filter_f = desire_f/resp_f
        filter_y = np.real(np.fft.ifft(filter_f))
        return filter_y
    
    def get_deconvolved_waveform(self, data):
        if self.deconv_filter is None:
            # print('Getting deconvolution filter...')
            self.deconv_filter = self.get_convolution_filter()
        return np.convolve(data, self.deconv_filter, 'same')  

    def run_filter(self, data, lowfreq=100.0, highfreq=100000.0, type='band', order=3):  
        x_filt = [] 
        for x in data:
            x_filt.append(self.get_filtered_waveform(self.Ch[0].Time, x, lowfreq, highfreq, order, type))
        return np.array(x_filt)
    
    def run_deconvolution(self, data, window=1000):
        x_deconv = []
        # count = int(self.length/window)
        for i,x in enumerate(data):
            # for j in range(count):
            #     print(i,j)
            #     x_cut = x[window*j:window*(j+1)]
            x_deconv.append(self.get_deconvolved_waveform(x))
        return np.array(x_deconv)
    
    def get_peaks(self, data, height=20, distance=1):
        # self.wvf_num = 0
        if self.peak_height:
            pass
        else:
            self.peak_height = []
            self.peak_pos = []
            self.peak_num = []
            self.peak_wvf_num = []
            self.wvf_num = 0
            
        for x in data:
            peaks,pdict = find_peaks(x, height=height, distance=distance)

            if len(peaks)>0:
                self.peak_height.extend(pdict['peak_heights'])
                self.peak_pos.extend(self.Ch[0].Time[peaks])
                self.peak_num.extend([len(peaks)])
                self.peak_wvf_num.extend([self.wvf_num]*len(peaks))
            else:
                self.peak_height.extend([-99999.99])
                self.peak_pos.extend([-99999.99])
                self.peak_num.extend([0])
                self.peak_wvf_num.extend([self.wvf_num])
            self.wvf_num += 1
    
    def get_peak_integral(self, data, window=[100]):
        if len(self.integral.keys())>0:
            pass
        else: 
            print('recreate')
            for wind in window: 
                self.integral[wind] = []
        for i,x in enumerate(self.peak_pos):
            if self.peak_height[i] == -99999.99:
                print(i,'skipping')
                continue
            else:
                print(i)
                for wind in window:
                    cut = np.where((self.Ch[0].Time>x-10) & (self.Ch[0].Time<x+wind))
                    self.integral[wind].append(np.sum(self.Ch[0].Amp[self.peak_wvf_num[i]][cut]))
        print(i)

    def clear(self):
        self.Ch[0].Amp = []
        self.Ch[0].Deconv = []
        self.Ch[0].Peak_height=[]