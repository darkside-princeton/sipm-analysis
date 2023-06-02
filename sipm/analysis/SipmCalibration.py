from sipm.analysis.AdvancedAnalyzer import AdvancedAnalyzer
from typing import Dict, List
import numpy as np
import scipy
from scipy.optimize import curve_fit
import sipm.util.functions as func

class SipmCalibration(AdvancedAnalyzer):
    def __init__(self, positions:List[str], channels:List[int], voltages:List[float], directory:str, metadata_dict:Dict, wf:bool, merge:bool, verbose:bool):
        super().__init__(directory, metadata_dict, wf, merge, verbose)
        self.positions = positions
        self.channels = channels
        self.voltages = voltages
        self.amp_hist = {}
        self.crosstalk = {}
        self.results = {'vbd': {}, 'dict': {}, 'ap_charge': {}, 'ap_prob': {}, 'gain': {}}
        for pos in self.positions:
            self.amp_hist[pos] = {}
            self.crosstalk[pos] = {}
            self.results['dict'][pos] = {}
            self.results['ap_charge'][pos] = {}
            self.results['ap_prob'][pos] = {}
            self.results['gain'][pos] = {}
            for ch in channels:
                self.amp_hist[pos][ch] = {}
                self.crosstalk[pos][ch] = {}
                self.results['dict'][pos][ch] = {}
                self.results['ap_charge'][pos][ch] = {}
                self.results['ap_prob'][pos][ch] = {}
                self.results['gain'][pos][ch] = {}
                for volt in voltages:
                    self.amp_hist[pos][ch][volt] = {}
                    self.crosstalk[pos][ch][volt] = {}
                
    def amplitude_analysis(self, boundary_par_dict, prom=70, wid=10, dist=15, nbins=1500, hist_range=(-1e2, 1.6e3)):
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
                    P_k = []
                    npe = len(p)
                    pe_cuts_in_bins = []
                    bound_par = boundary_par_dict[pos][ch][volt]
                    for ipe in range(npe):
                        if ipe == 0:
                            pe_cuts_in_bins.append(int(1.5*p[0]-0.5*p[1]))
                        else:
                            pe_cuts_in_bins.append(int(bound_par*p[ipe]+(1-bound_par)*p[ipe-1]))
                            # if pos == 'bottom' and ch == 3:
                            #     pe_cuts_in_bins.append(int(0.35*p[ipe]+0.65*p[ipe-1]))
                            # elif pos == 'bottom' and ch == 0:
                            #     pe_cuts_in_bins.append(int(0.45*p[ipe]+0.55*p[ipe-1]))
                            # else:
                            #     pe_cuts_in_bins.append(int(0.5*(p[ipe]+p[ipe-1])))
                        # if ipe > 0:
                            P_k.append([np.sum(self.amp_hist[pos][ch][volt]['hist'][pe_cuts_in_bins[ipe-1]:pe_cuts_in_bins[ipe]])/nevents,
                                    np.sqrt(np.sum(self.amp_hist[pos][ch][volt]['hist'][pe_cuts_in_bins[ipe-1]:pe_cuts_in_bins[ipe]]))/nevents])
                    self.amp_hist[pos][ch][volt]['boundaries'] = list(
                        np.array(pe_cuts_in_bins)*bin_width+hist_range[0])
                    # Save P_k for Vinogradov fit
                    self.crosstalk[pos][ch][volt]['y'] = np.array(P_k)[:, 0]
                    self.crosstalk[pos][ch][volt]['yerr'] = np.array(P_k)[:, 1]
                    self.crosstalk[pos][ch][volt]['x'] = np.arange(len(P_k))

    def crosstalk_analysis(self):
        for pos in self.positions:
            for ch in self.channels:
                for volt in self.voltages:
                    # Do Vinogradov fit
                    popt, pcov = curve_fit(func.compound_poisson,
                                        self.crosstalk[pos][ch][volt]['x'],
                                        self.crosstalk[pos][ch][volt]['y'],
                                        p0=[2, 0.2], sigma=self.crosstalk[pos][ch][volt]['yerr'], maxfev=10000)
                    # Save fit results
                    self.crosstalk[pos][ch][volt]['par'] = popt
                    self.crosstalk[pos][ch][volt]['cov'] = pcov
                    self.crosstalk[pos][ch][volt]['dict'] = popt[1]
                    self.crosstalk[pos][ch][volt]['dict_err'] = func.error_distance(df=2, sigma=1)*np.sqrt(pcov[1, 1])
                    print(f'{pos} ch{ch} {volt}V P_dict = {self.crosstalk[pos][ch][volt]["dict"]:.4f} +/- {self.crosstalk[pos][ch][volt]["dict_err"]:.4f}')
                self.results['dict'][pos][ch]['x'] = self.voltages
                self.results['dict'][pos][ch]['y'] = [self.crosstalk[pos][ch][volt]['dict'] for volt in self.voltages]
                self.results['dict'][pos][ch]['yerr'] = [self.crosstalk[pos][ch][volt]['dict_err'] for volt in self.voltages]
