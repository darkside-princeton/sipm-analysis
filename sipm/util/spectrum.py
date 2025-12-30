from scipy.stats import norm
import numpy as np
import pandas as pd
from datetime import datetime
from scipy.optimize import least_squares
from scipy.stats import norm

def load_data(file_list, channel_list, volt):
    df_ch = {}
    dt = 0
    for ch in channel_list:
        df_ch[ch] = []
        for f in file_list:
            df = pd.read_hdf(f, key=f'{volt}/{ch}')
            df_ch[ch].append(df)
            if ch==-1:
                dt = datetime(*np.array(df['start_datetime'][:6]).astype(int))
            df = None
        df_ch[ch] = pd.concat(df_ch[ch]).sort_index()
        if ch!=-1:
            df_ch[ch].rename(columns=lambda x: x+f'_{ch}',inplace=True)            
    df_ch = pd.concat(list(df_ch.values()),axis=1)
    return df_ch, dt

def convert_to_pe(df_data,df_pe_scale,channels):
    pe_ch = df_data[[f'integral_10p00us_{ch}' for ch in range(8)]]
    pe_ch.columns = np.arange(8)
    pe_ch = pe_ch/df_pe_scale
    df_pe = pe_ch.iloc[:,channels].sum(axis=1)
    df_pe.name = 'pe_' + ''.join([str(c) for c in channels])
    pe_ch = None
    return df_pe

def get_bsl_filt(df_data,channels,rms_thre_ch):
    df_bsl_filt = {}
    for ch in channels:
        df_bsl_filt[f'bsl_filt_{ch}'] = df_data[f'baseline_rms_{ch}']<rms_thre_ch[ch]
        # print(df_bsl_filt[f'bsl_filt_{ch}'].sum())
    df_bsl_filt = pd.DataFrame(df_bsl_filt)
    df_bsl_filt['bsl_filt'] =  df_bsl_filt[[f'bsl_filt_{ch}' for ch in channels]].all(axis=1)
    return df_bsl_filt

def get_bsl_rms_hist(df_data,channels,bins,range):
    hist_ch = {}
    for ch in channels:
        hist_ch[ch] = np.histogram(df_data[f'baseline_rms_{ch}'], bins=bins, range=range)
    return hist_ch

def get_bsl_mean_hist(df_data,channels,bsl_filt,bins,range):
    hist_ch = {}
    for ch in channels:
        if bsl_filt:
            hist_ch[ch] = np.histogram(df_data.loc[df_data[f'bsl_filt_{ch}'],f'baseline_mean_{ch}'], bins=bins, range=range)
        else:
            hist_ch[ch] = np.histogram(df_data[f'baseline_mean_{ch}'], bins=bins, range=range)
    return hist_ch

def get_sat_filt(df_data,channels,sat_list):
    df_sat_filt = {}
    for ch in channels:
        df_sat_filt[f'sat_filt_{ch}'] = df_data[f'amplitude_{ch}']<sat_list[ch]
    df_sat_filt = pd.DataFrame(df_sat_filt)
    df_sat_filt['sat_filt'] =  df_sat_filt[[f'sat_filt_{ch}' for ch in range(8)]].all(axis=1)
    return df_sat_filt

def get_fp_filt(df_data,channel_groups,fp_dict,length):
    df_fp_filt = {}
    for g in channel_groups:
        df_fp_filt[f'fp_filt_{g}'] = (df_data[f'fprompt_{length}_{g}']>fp_dict[g][0]) & (df_data[f'fprompt_{length}_{g}']<fp_dict[g][1])
    df_fp_filt = pd.DataFrame(df_fp_filt)
    df_fp_filt['fp_filt'] =  df_fp_filt[[f'fp_filt_{g}' for g in channel_groups]].all(axis=1)
    return df_fp_filt

def get_fp_hist2d(df_data,channel_groups,length,bins_xy,range_xy,log_scale):
    hist2d_g = {}
    for g in channel_groups:
        filt = df_data['bsl_filt']
        hist2d_g[g] = list(np.histogram2d(
            df_data[f'pe_{g}'][filt], 
            df_data[f'fprompt_{length}_{g}'][filt], 
            bins=bins_xy, range=range_xy
        ))
        temp = hist2d_g[g][0].T
        hist2d_g[g][0], hist2d_g[g][1] = np.meshgrid(
            hist2d_g[g][1], hist2d_g[g][2]
        )
        if log_scale:
            hist2d_g[g][2] = np.log10(temp)
        else:
            hist2d_g[g][2] = temp
    return hist2d_g

def get_fs_nofs_hist2d(df_data,bsl_filt,fp_filt,sat_filt,bins_xy,range_xy,log_scale):
    filt = np.ones(df_data.shape[0])
    if bsl_filt:
        filt = filt & df_data['bsl_filt']
    if fp_filt:
        filt = filt & df_data['fp_filt_01234567']
    if sat_filt:
        filt = filt & df_data['sat_filt']
    hist2d = list(np.histogram2d(
        df_data['nofs_pe'][filt], 
        df_data['fs_pe'][filt], 
        bins=bins_xy, range=range_xy
    ))
    temp = hist2d[0].T
    hist2d[0], hist2d[1] = np.meshgrid(
        hist2d[1], hist2d[2]
    )
    if log_scale:
        hist2d[2] = np.log10(temp)
    else:
        hist2d[2] = temp
    return hist2d

class DetectorSmearing():
    def __init__(self, fano, hist_in):
        self.hist_in_bin = hist_in[1]
        self.hist_in_width = self.hist_in_bin[1]-self.hist_in_bin[0]
        self.hist_in_count = hist_in[0]/np.sum(hist_in[0])/self.hist_in_width # normalized
        self.fano = fano

    def response_matrix(self, pe_in, pe_out, pde):
        pe_det = pde*pe_in
        return norm.pdf(pe_out, loc=pe_det, scale=np.sqrt((self.fano-pde)*pe_det))

    def get_spectrum(self, bins, n, pde):
        self.hist_out_bin = bins
        bin_cen_in = (self.hist_in_bin[1:]+self.hist_in_bin[:-1])/2
        bin_cen_out = (self.hist_out_bin[1:]+self.hist_out_bin[:-1])/2
        pe_out_mesh, pe_in_mesh = np.meshgrid(bin_cen_out, bin_cen_in)
        self.hist_out_count = self.hist_in_count @ self.response_matrix(pe_in_mesh, pe_out_mesh, pde)
        return n*self.hist_out_count, self.hist_out_bin
        
    def fit_to_data(self, data, sigma, x0, fit_range):
        def residual(x, *args, **kwargs):
            n, pde = x[0], x[1]
            data_y, data_bin = args[0], args[1]
            bincen = (data_bin[1:] + data_bin[:-1])/2
            ans = (data_y-self.get_spectrum(data_bin, n, pde)[0])
            sigma[sigma==0] = 1
            ans = ans/sigma
            mask = (bincen<fit_range[1]) & (bincen>fit_range[0])
            return ans[mask]
        res = least_squares(residual,x0=x0,args=data,bounds=(0,np.inf))
        cov = np.linalg.inv(res.jac.T@res.jac)
        return res.x, cov
        