from sipm.analysis.AdvancedAnalyzer import AdvancedAnalyzer
from typing import Dict, List
import numpy as np
import scipy
from scipy.optimize import curve_fit
import sipm.util.functions as func
import os
import csv

class SipmCalibration(AdvancedAnalyzer):
    """The class for SiPM calibration analysis
    """
    def __init__(self, positions:List[str], channels:List[int], voltages:List[float], directory:str, metadata_dict:Dict, wf:bool, merge:bool, verbose:bool):
        """SipmCalibration constructor.

        Args:
            positions (List[str]): A list of positions (e.g. ['top','bottom'])
            channels (List[int]): A list of channels (e.g. [0,1,2,3])
            voltages (List[float]): A list of voltages (e.g. [63,65,67,69,71])
            directory (str): Directory containing processed HDF5 files
            metadata_dict (Dict): Metadata used to specify the file names. Arranged into a nested dictionary.
            wf (bool): Whether the file name ends with '_wf'
            merge (bool): Whether to merge different runs
            verbose (bool): Whether to print out more information
        """
        super().__init__(directory, metadata_dict, wf, merge, verbose)
        self.positions = positions
        self.channels = channels
        self.voltages = voltages
        self.amp_hist = {} # filtered amplitude histogram
        self.crosstalk = {} # DiCT fits (Probability-PE)
        self.charge_hist = {} # integrated charge histogram
        self.charge_fits = {} # Gaussian fits to PE peaks
        self.gain_peak_fits = {} # linear fits to peak position vs PE
        self.gain_avg_fits = {} # linear fits to averaged charge vs PE
        self.ap_charge = {} # average afterpulsing charge
        self.ap_prob = {} # afterpulsing probability
        self.ap_prob_fits = {} # fits for ap probability
        self.vbd_fits = {} # fits for breakdown voltage
        self.results = {'vbd': {}, 'dict': {}, 'ap_charge': {}, 'ap_prob': {}, 'gain': {}}
        for pos in self.positions:
            self.amp_hist[pos] = {}
            self.crosstalk[pos] = {}
            self.charge_hist[pos] = {}
            self.charge_fits[pos] = {} 
            self.gain_peak_fits[pos] = {}
            self.gain_avg_fits[pos] = {} 
            self.ap_charge[pos] = {}
            self.ap_prob[pos] = {}
            self.ap_prob_fits[pos] = {}
            self.vbd_fits[pos] = {}
            self.results['dict'][pos] = {}
            self.results['ap_charge'][pos] = {}
            self.results['ap_prob'][pos] = {}
            self.results['gain'][pos] = {}
            self.results['vbd'][pos] = {}
            for ch in channels:
                self.amp_hist[pos][ch] = {}
                self.crosstalk[pos][ch] = {}
                self.charge_hist[pos][ch] = {}
                self.charge_fits[pos][ch] = {} 
                self.gain_peak_fits[pos][ch] = {}
                self.gain_avg_fits[pos][ch] = {} 
                self.ap_charge[pos][ch] = {}
                self.ap_prob[pos][ch] = {}
                self.ap_prob_fits[pos][ch] = {}
                self.vbd_fits[pos][ch] = {}
                self.results['dict'][pos][ch] = {}
                self.results['ap_charge'][pos][ch] = {}
                self.results['ap_prob'][pos][ch] = {}
                self.results['gain'][pos][ch] = {}
                self.results['vbd'][pos][ch] = {}
                for volt in voltages:
                    self.amp_hist[pos][ch][volt] = {}
                    self.crosstalk[pos][ch][volt] = {}
                    self.charge_hist[pos][ch][volt] = {}
                    self.charge_fits[pos][ch][volt] = {} 
                    self.gain_peak_fits[pos][ch][volt] = {}
                    self.gain_avg_fits[pos][ch][volt] = {} 
                    self.ap_charge[pos][ch][volt] = {}
                    self.ap_prob[pos][ch][volt] = {}
                    self.ap_prob_fits[pos][ch][volt] = {}
                
    def amplitude_analysis(self, boundary_par_dict, prom=70, wid=10, dist=15, nbins=1500, hist_range=(-1e2, 1.6e3)):
        """Analyze filtered amplitude histograms.

        Args:
            boundary_par_dict (Dict): Parameter to set the boundary between different PEs. 0.5=middle point between two peaks. 0=the left peak. 1=the right peak. Arranged into a nested dictionary with the same structure as self.metadata.
            prom (int, optional): Prominence for the peak finder. Defaults to 70.
            wid (int, optional): Width for the peak finder. Defaults to 10.
            dist (int, optional): Distance for the peak finder. Defaults to 15.
            nbins (int, optional): Number of bins for the histograms. Defaults to 1500.
            hist_range (tuple, optional): Range for the histogram. Defaults to (-1e2, 1.6e3).
        """
        # Generate histograms
        bin_width = (hist_range[1]-hist_range[0])/nbins
        for pos in self.positions:
            for ch in self.channels:
                for volt in self.voltages:
                    nevents = np.sum(self.data[pos][ch][volt]['data']['bsl_cut'])
                    self.amp_hist[pos][ch][volt]['hist'], self.amp_hist[pos][ch][volt]['bins'] = np.histogram(
                        self.data[pos][ch][volt]['data']['amplitude_trig'].loc[self.data[pos][ch][volt]['data']['bsl_cut']], 
                        bins=nbins, range=hist_range
                    )
                    # find PE peaks in histogram
                    p, pdict = scipy.signal.find_peaks(
                        self.amp_hist[pos][ch][volt]['hist'], prominence=prom, width=wid, distance=dist)
                    # discriminate different PE counts and calculate probability distribution P_k
                    self.data[pos][ch][volt]['data']['pe'] = np.zeros(self.data[pos][ch][volt]['data'].shape[0]).astype(int)
                    P_k = []
                    npe = len(p)
                    pe_cuts_in_bins = []
                    bound_par = boundary_par_dict[pos][ch][volt]
                    for ipe in range(npe):
                        if ipe == 0:
                            pe_cuts_in_bins.append(int(1.5*p[0]-0.5*p[1]))
                        else:
                            pe_cuts_in_bins.append(int(bound_par*p[ipe]+(1-bound_par)*p[ipe-1]))
                            self.data[pos][ch][volt]['data']['pe'] += (self.data[pos][ch][volt]['data']['amplitude_trig']>pe_cuts_in_bins[-1]*bin_width+hist_range[0]).astype(int)
                            P_k.append([np.sum(self.amp_hist[pos][ch][volt]['hist'][pe_cuts_in_bins[ipe-1]:pe_cuts_in_bins[ipe]])/nevents,
                                    np.sqrt(np.sum(self.amp_hist[pos][ch][volt]['hist'][pe_cuts_in_bins[ipe-1]:pe_cuts_in_bins[ipe]]))/nevents])
                    self.amp_hist[pos][ch][volt]['boundaries'] = list(
                        np.array(pe_cuts_in_bins)*bin_width+hist_range[0])
                    # Save P_k for Vinogradov fit
                    self.crosstalk[pos][ch][volt]['y'] = np.array(P_k)[:, 0]
                    self.crosstalk[pos][ch][volt]['yerr'] = np.array(P_k)[:, 1]
                    self.crosstalk[pos][ch][volt]['x'] = np.arange(len(P_k))

    def crosstalk_analysis(self):
        """Perform Vinogradov fit to obtain the direct crosstalk probability.
        """
        for pos in self.positions:
            for ch in self.channels:
                for volt in self.voltages:
                    # Do Vinogradov fit
                    self.crosstalk[pos][ch][volt]['par'], self.crosstalk[pos][ch][volt]['cov'] = curve_fit(
                        func.compound_poisson,
                        self.crosstalk[pos][ch][volt]['x'],
                        self.crosstalk[pos][ch][volt]['y'],
                        p0=[2, 0.2], sigma=self.crosstalk[pos][ch][volt]['yerr'], maxfev=10000)
                    # Save fit results
                    self.crosstalk[pos][ch][volt]['dict'] = self.crosstalk[pos][ch][volt]['par'][1]
                    self.crosstalk[pos][ch][volt]['dict_err'] = func.error_distance(df=2, sigma=1)*np.sqrt(self.crosstalk[pos][ch][volt]['cov'][1, 1])
                    print(f'{pos} ch{ch} {volt}V P_dict = {self.crosstalk[pos][ch][volt]["dict"]:.4f} +/- {self.crosstalk[pos][ch][volt]["dict_err"]:.4f}')
                self.results['dict'][pos][ch]['bias'] = self.voltages
                self.results['dict'][pos][ch]['dict'] = [self.crosstalk[pos][ch][volt]['dict'] for volt in self.voltages]
                self.results['dict'][pos][ch]['dict_err'] = [self.crosstalk[pos][ch][volt]['dict_err'] for volt in self.voltages]

    def charge_analysis(self,nbins=1000, hist_range=(-2e3, 6e3), fit_range_thre = (0.02, 0.3)):
        for pos in self.positions:
            for ch in self.channels:
                for volt in self.voltages:
                    range_min = hist_range[0]
                    range_max = hist_range[1]+(volt-63)*1e3
                    bin_width = (range_max-range_min)/nbins
                    for pe in np.arange(1,len(self.amp_hist[pos][ch][volt]['boundaries'])-1):
                        self.charge_hist[pos][ch][volt][pe] = {}
                        self.charge_fits[pos][ch][volt][pe] = {}
                        self.ap_prob[pos][ch][volt][pe] = {}
                        # Generate histograms
                        selected_charges = self.data[pos][ch][volt]['data']['integral_5p00us'].loc[(self.data[pos][ch][volt]['data']['bsl_cut'])&(self.data[pos][ch][volt]['data']['pe']==pe)&(self.data[pos][ch][volt]['data']['integral_5p00us']<range_max)&(self.data[pos][ch][volt]['data']['integral_5p00us']>range_min)]
                        self.charge_hist[pos][ch][volt][pe]['hist'], self.charge_hist[pos][ch][volt][pe]['bins'] = np.histogram(
                            selected_charges, bins=nbins, range=(range_min, range_max))
                        # find appropriate fit range
                        peak_bin = np.argmax(self.charge_hist[pos][ch][volt][pe]['hist'])
                        fit_min = peak_bin
                        while self.charge_hist[pos][ch][volt][pe]['hist'][fit_min] > fit_range_thre[0]*self.charge_hist[pos][ch][volt][pe]['hist'][peak_bin]:
                            fit_min -= 1
                        if self.charge_hist[pos][ch][volt][pe]['hist'][fit_min] == 0:
                            fit_min += 1
                        self.charge_fits[pos][ch][volt][pe]['min_bin'] = fit_min
                        fit_max = peak_bin
                        while self.charge_hist[pos][ch][volt][pe]['hist'][fit_max] > fit_range_thre[1]*self.charge_hist[pos][ch][volt][pe]['hist'][peak_bin]:
                            fit_max += 1
                        if self.charge_hist[pos][ch][volt][pe]['hist'][fit_max] == 0:
                            fit_max -= 1
                        self.charge_fits[pos][ch][volt][pe]['max_bin'] = fit_max
                        # Gaussian fits
                        self.charge_fits[pos][ch][volt][pe]['par'], self.charge_fits[pos][ch][volt][pe]['cov'] = curve_fit(
                            func.gauss_normalized,
                            0.5*(self.charge_hist[pos][ch][volt][pe]['bins'][fit_min:fit_max]+self.charge_hist[pos][ch][volt][pe]['bins'][fit_min+1:fit_max+1]),
                            self.charge_hist[pos][ch][volt][pe]['hist'][fit_min:fit_max],
                            p0=[self.charge_hist[pos][ch][volt][pe]['hist'][peak_bin]*(fit_max-fit_min)*bin_width/3,
                                peak_bin*bin_width+range_min,
                                (fit_max-fit_min)*bin_width/3],
                            sigma=np.sqrt(self.charge_hist[pos][ch][volt][pe]['hist'][fit_min:fit_max]),
                            maxfev=100000)
                        # Store peak positions and averages
                        self.charge_fits[pos][ch][volt][pe]['Ipeak'] = self.charge_fits[pos][ch][volt][pe]['par'][1]
                        self.charge_fits[pos][ch][volt][pe]['Ipeak_err'] = func.error_distance(
                            df=3, sigma=1)*np.sqrt(self.charge_fits[pos][ch][volt][pe]['cov'][1, 1])
                        self.charge_fits[pos][ch][volt][pe]['Iavg'] = np.mean(selected_charges)
                        self.charge_fits[pos][ch][volt][pe]['Iavg_err'] = np.std(selected_charges)/np.sqrt(len(selected_charges))
                        # Store afterpulse probability
                        self.ap_prob[pos][ch][volt][pe]['prob'] = self.charge_fits[pos][ch][volt][pe]['par'][0]/bin_width/len(selected_charges)
                        self.ap_prob[pos][ch][volt][pe]['prob_err'] = np.sqrt(self.ap_prob[pos][ch][volt][pe]['prob']*(1-self.ap_prob[pos][ch][volt][pe]['prob'])/len(selected_charges))

    def gain_analysis(self):
        for pos in self.positions:
            for ch in self.channels:
                for volt in self.voltages:
                    # Qpeak fit
                    self.gain_peak_fits[pos][ch][volt]['x'] = np.array(list(self.charge_fits[pos][ch][volt].keys()))
                    self.gain_peak_fits[pos][ch][volt]['y'] = np.array([self.charge_fits[pos][ch][volt][pe]['Ipeak'] for pe in self.charge_fits[pos][ch][volt].keys()])
                    self.gain_peak_fits[pos][ch][volt]['yerr'] = np.array([self.charge_fits[pos][ch][volt][pe]['Ipeak_err'] for pe in self.charge_fits[pos][ch][volt].keys()])
                    self.gain_peak_fits[pos][ch][volt]['par'], self.gain_peak_fits[pos][ch][volt]['cov'] = curve_fit(
                        func.line_simple,
                        self.gain_peak_fits[pos][ch][volt]['x'],
                        self.gain_peak_fits[pos][ch][volt]['y'],
                        sigma=self.gain_peak_fits[pos][ch][volt]['yerr'],
                        p0=[500, 0],
                        maxfev=10000)
                    self.gain_peak_fits[pos][ch][volt]['Qpeak'] = self.gain_peak_fits[pos][ch][volt]['par'][0]
                    self.gain_peak_fits[pos][ch][volt]['Qpeak_err'] = func.error_distance(df=2, sigma=1)*np.sqrt(self.gain_peak_fits[pos][ch][volt]['cov'][0, 0])
                    # Qavg fit
                    self.gain_avg_fits[pos][ch][volt]['x'] = np.array(list(self.charge_fits[pos][ch][volt].keys()))
                    self.gain_avg_fits[pos][ch][volt]['y'] = np.array([self.charge_fits[pos][ch][volt][pe]['Iavg'] for pe in self.charge_fits[pos][ch][volt].keys()])
                    self.gain_avg_fits[pos][ch][volt]['yerr'] = np.array([self.charge_fits[pos][ch][volt][pe]['Iavg_err'] for pe in self.charge_fits[pos][ch][volt].keys()])
                    self.gain_avg_fits[pos][ch][volt]['par'], self.gain_avg_fits[pos][ch][volt]['cov'] = curve_fit(
                        func.line_simple,
                        self.gain_avg_fits[pos][ch][volt]['x'],
                        self.gain_avg_fits[pos][ch][volt]['y'],
                        sigma=self.gain_avg_fits[pos][ch][volt]['yerr'],
                        p0=[500, 0],
                        maxfev=10000)
                    self.gain_avg_fits[pos][ch][volt]['Qavg'] = self.gain_avg_fits[pos][ch][volt]['par'][0]
                    self.gain_avg_fits[pos][ch][volt]['Qavg_err'] = func.error_distance(df=2, sigma=1)*np.sqrt(self.gain_avg_fits[pos][ch][volt]['cov'][0, 0])
                self.results['gain'][pos][ch]['bias'] = self.voltages
                self.results['gain'][pos][ch]['gain'] = [self.gain_peak_fits[pos][ch][volt]['Qpeak'] for volt in self.voltages]
                self.results['gain'][pos][ch]['gain_err'] = [self.gain_peak_fits[pos][ch][volt]['Qpeak_err'] for volt in self.voltages]
    
    def afterpulse_analysis(self):
        for pos in self.positions:
            for ch in self.channels:
                for volt in self.voltages:
                    # Afterpulse charge
                    self.ap_charge[pos][ch][volt]['Qap'] = self.gain_avg_fits[pos][ch][volt]['Qavg']/self.gain_peak_fits[pos][ch][volt]['Qpeak'] - 1
                    self.ap_charge[pos][ch][volt]['Qap_err'] = (1+self.ap_charge[pos][ch][volt]['Qap'])*np.sqrt((self.gain_avg_fits[pos][ch][volt]['Qavg_err']/self.gain_avg_fits[pos][ch][volt]['Qavg'])**2+(self.gain_peak_fits[pos][ch][volt]['Qpeak_err']/self.gain_peak_fits[pos][ch][volt]['Qpeak'])**2)
                    self.ap_charge[pos][ch][volt]['Qap'] *= 1-self.crosstalk[pos][ch][volt]['dict']  # correct for APxDiCT
                    self.ap_charge[pos][ch][volt]['Qap_err'] *= 1-self.crosstalk[pos][ch][volt]['dict']  # correct for APxDiCT
                    # Afterpulse probability fit
                    self.ap_prob_fits[pos][ch][volt]['x'] = np.array(list(self.ap_prob[pos][ch][volt].keys()))
                    self.ap_prob_fits[pos][ch][volt]['y'] = [self.ap_prob[pos][ch][volt][pe]['prob'] for pe in self.ap_prob[pos][ch][volt].keys()]
                    self.ap_prob_fits[pos][ch][volt]['yerr'] = [self.ap_prob[pos][ch][volt][pe]['prob_err'] for pe in self.ap_prob[pos][ch][volt].keys()]
                    self.ap_prob_fits[pos][ch][volt]['par'], self.ap_prob_fits[pos][ch][volt]['cov'] = curve_fit(
                        func.power,
                        self.ap_prob_fits[pos][ch][volt]['x'],
                        self.ap_prob_fits[pos][ch][volt]['y'],
                        sigma=self.ap_prob_fits[pos][ch][volt]['yerr'],
                        p0=[0.9],
                        maxfev=10000)
                    self.ap_prob_fits[pos][ch][volt]['Pap'] = 1-self.ap_prob_fits[pos][ch][volt]['par'][0]
                    self.ap_prob_fits[pos][ch][volt]['Pap_err'] = func.error_distance(df=1, sigma=1)*np.sqrt(self.ap_prob_fits[pos][ch][volt]['cov'][0, 0])
                # Qap-Vbias
                self.results['ap_charge'][pos][ch]['bias'] = self.voltages
                self.results['ap_charge'][pos][ch]['ap_charge'] = [self.ap_charge[pos][ch][volt]['Qap'] for volt in self.voltages]
                self.results['ap_charge'][pos][ch]['ap_charge_err'] = [self.ap_charge[pos][ch][volt]['Qap_err'] for volt in self.voltages]
                # Pap-Vbias
                self.results['ap_prob'][pos][ch]['bias'] = self.voltages
                self.results['ap_prob'][pos][ch]['ap_prob'] = [self.ap_prob_fits[pos][ch][volt]['Pap'] for volt in self.voltages]
                self.results['ap_prob'][pos][ch]['ap_prob_err'] = [self.ap_prob_fits[pos][ch][volt]['Pap_err'] for volt in self.voltages]
    
    def breakdown_analysis(self, init_pars=[100,55]):
        # Fitting for breakdown voltage
        for pos in self.positions:
            for ch in self.channels:
                self.vbd_fits[pos][ch]['x'] = self.voltages
                self.vbd_fits[pos][ch]['y'] = [self.gain_peak_fits[pos][ch][volt]['Qpeak'] for volt in self.voltages]
                self.vbd_fits[pos][ch]['yerr'] = [self.gain_peak_fits[pos][ch][volt]['Qpeak_err'] for volt in self.voltages]
                self.vbd_fits[pos][ch]['par'], self.vbd_fits[pos][ch]['cov'] = curve_fit(
                    func.line,
                    self.vbd_fits[pos][ch]['x'],
                    self.vbd_fits[pos][ch]['y'],
                    sigma=self.vbd_fits[pos][ch]['yerr'],
                    p0=init_pars,
                    maxfev=10000)
                self.vbd_fits[pos][ch]['vbd'] = self.vbd_fits[pos][ch]['par'][1]
                self.vbd_fits[pos][ch]['vbd_err'] = func.error_distance(df=2, sigma=1)*np.sqrt(self.vbd_fits[pos][ch]['cov'][1, 1])
                print(f'{pos} ch{ch} Vbd = {self.vbd_fits[pos][ch]["vbd"]:.2f} +/- {self.vbd_fits[pos][ch]["vbd_err"]:.2f} V')
                self.results['vbd'][pos][ch]['vbd_sipm'] = self.vbd_fits[pos][ch]['vbd']/2
                self.results['vbd'][pos][ch]['vbd_sipm_err'] = self.vbd_fits[pos][ch]['vbd_err']/2
                self.results['dict'][pos][ch]['ov'] = np.array(self.results['dict'][pos][ch]['bias'])/2-self.results['vbd'][pos][ch]['vbd_sipm']
                self.results['dict'][pos][ch]['ov_err'] = np.ones(self.results['dict'][pos][ch]['ov'].shape[0])*self.results['vbd'][pos][ch]['vbd_sipm_err']
                self.results['ap_charge'][pos][ch]['ov'] = np.array(self.results['ap_charge'][pos][ch]['bias'])/2-self.results['vbd'][pos][ch]['vbd_sipm']
                self.results['ap_charge'][pos][ch]['ov_err'] = np.ones(self.results['ap_charge'][pos][ch]['ov'].shape[0])*self.results['vbd'][pos][ch]['vbd_sipm_err']
                self.results['ap_prob'][pos][ch]['ov'] = np.array(self.results['ap_prob'][pos][ch]['bias'])/2-self.results['vbd'][pos][ch]['vbd_sipm']
                self.results['ap_prob'][pos][ch]['ov_err'] = np.ones(self.results['ap_prob'][pos][ch]['ov'].shape[0])*self.results['vbd'][pos][ch]['vbd_sipm_err']
                self.results['gain'][pos][ch]['ov'] = np.array(self.results['gain'][pos][ch]['bias'])/2-self.results['vbd'][pos][ch]['vbd_sipm']
                self.results['gain'][pos][ch]['ov_err'] = np.ones(self.results['gain'][pos][ch]['ov'].shape[0])*self.results['vbd'][pos][ch]['vbd_sipm_err']
                
    def write_to_csv(self, name):
        if not os.path.exists(f'data/{name}'):
            os.makedirs(f'data/{name}')
        for volt in self.voltages:
            with open(f'data/{name}/{name}_{volt}V.csv', 'w') as f:
                w = csv.writer(f)
                w.writerow(['CH', 'A1min', 'A1max', 'p', 'p_err', 'Qavg', 'Qavg_err', 'Qpeak', 'Qpeak_err', 'Qap', 'Qap_err', 'bsl_rms'])
                for pos in self.positions:
                    for ch in self.channels:
                        row = [f'{pos[0].upper()}{ch}']
                        row += [str(self.amp_hist[pos][ch][volt]['boundaries'][1]),
                                str(self.amp_hist[pos][ch][volt]['boundaries'][2])]
                        row += [str(self.crosstalk[pos][ch][volt]['dict']),
                                str(self.crosstalk[pos][ch][volt]['dict_err'])]
                        row += [str(self.gain_avg_fits[pos][ch][volt]['Qavg']),
                                str(self.gain_avg_fits[pos][ch][volt]['Qavg_err'])]
                        row += [str(self.gain_peak_fits[pos][ch][volt]['Qpeak']),
                                str(self.gain_peak_fits[pos][ch][volt]['Qpeak_err'])]
                        row += [str(self.ap_charge[pos][ch][volt]['Qap']),
                                str(self.ap_charge[pos][ch][volt]['Qap_err'])]
                        row += [str(self.bsl_rms_thre[pos][ch][volt])]
                        w.writerow(row)