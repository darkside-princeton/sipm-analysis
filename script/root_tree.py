import numpy as np
import sipm.sipm as sipm
import sipm.dataset as ds
import sipm.beta as beta
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import sipm.constants as const
import ROOT
from numpy.random import normal

dataDir = '/scratch/gpfs/GALBIATI/data/sipm/reflector_studies/'
dir1118top = ['2022-11-20/2022-11-20_volt_63_pos_top_light_scintillation_coinc_000_cond_gamma/', '2022-11-18/2022-11-18_volt_63_pos_top_light_scintillation_coinc_000_cond_no_gamma/',
              '2022-11-20/2022-11-20_volt_65_pos_top_light_scintillation_coinc_000_cond_gamma/', '2022-11-18/2022-11-18_volt_65_pos_top_light_scintillation_coinc_000_cond_no_gamma/',
              '2022-11-20/2022-11-20_volt_67_pos_top_light_scintillation_coinc_000_cond_gamma/', '2022-11-18/2022-11-18_volt_67_pos_top_light_scintillation_coinc_000_cond_no_gamma/',
              '2022-11-20/2022-11-20_volt_69_pos_top_light_scintillation_coinc_000_cond_gamma/', '2022-11-18/2022-11-18_volt_69_pos_top_light_scintillation_coinc_000_cond_no_gamma/',
              '2022-11-20/2022-11-20_volt_71_pos_top_light_scintillation_coinc_000_cond_gamma/', '2022-11-18/2022-11-18_volt_71_pos_top_light_scintillation_coinc_000_cond_no_gamma/']
dir1118bot = ['2022-11-20/2022-11-20_volt_63_pos_bottom_light_scintillation_coinc_000_cond_gamma/', '2022-11-18/2022-11-18_volt_63_pos_bottom_light_scintillation_coinc_000_cond_no_gamma/',
              '2022-11-20/2022-11-20_volt_65_pos_bottom_light_scintillation_coinc_000_cond_gamma/', '2022-11-18/2022-11-18_volt_65_pos_bottom_light_scintillation_coinc_000_cond_no_gamma/',
              '2022-11-20/2022-11-20_volt_67_pos_bottom_light_scintillation_coinc_000_cond_gamma/', '2022-11-18/2022-11-18_volt_67_pos_bottom_light_scintillation_coinc_000_cond_no_gamma/',
              '2022-11-20/2022-11-20_volt_69_pos_bottom_light_scintillation_coinc_000_cond_gamma/', '2022-11-18/2022-11-18_volt_69_pos_bottom_light_scintillation_coinc_000_cond_no_gamma/',
              '2022-11-20/2022-11-20_volt_71_pos_bottom_light_scintillation_coinc_000_cond_gamma/', '2022-11-18/2022-11-18_volt_71_pos_bottom_light_scintillation_coinc_000_cond_no_gamma/']
bias = [63, 65, 67, 69, 71]
gain = [[556.758,548.693,548.862,541.270], [556.997, 503.579, 549.769, 558.519]] # gain = [[T0,T1,T2,T3],[B0,B1,B2,B3]]

ds1118top = []
for i,dir in enumerate(dir1118top):
    data = ds.Dataset('', pol=-1, channels=range(4), spe=gain[0])
    for j in range(20):
        for ch in range(4):
            subdir = '{}{}{}/'.format(dataDir, dir, j)
            print(subdir)
            data.ch[ch].path = subdir
            data.ch[ch].read_data(root_file=)
            data.ch[ch].baseline_subtraction()
            data.ch[ch].get_integral(prompt=0.5, long=5)
            data.ch[ch].clear()
    data.get_summed_integral_pe()
    ds1118top.append(data)

ds1118bot = []
for i,dir in enumerate(dir1118bot):
    data = ds.Dataset('', pol=-1, channels=range(4), spe=gain[1])
    for j in range(20):
        for ch in range(4):
            subdir = '{}{}{}/'.format(dataDir, dir, j)
            print(subdir)
            data.ch[ch].path = subdir
            data.ch[ch].read_data(root_file=)
            data.ch[ch].baseline_subtraction()
            data.ch[ch].get_integral(prompt=0.5, long=5)
            data.ch[ch].clear()
    data.get_summed_integral_pe()
    ds1118bot.append(data)