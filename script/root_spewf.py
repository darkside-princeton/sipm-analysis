import numpy as np
import sipm.sipm as sipm
import sipm.dataset as ds
import sipm.constants as const
import ROOT
import sys
import os
from array import array
import csv

if __name__ == '__main__':
    #usage: python root_calibration.py index
    index = int(sys.argv[1])
    data_dir = '/scratch/gpfs/GALBIATI/data/sipm/reflector_studies/'
    direc = [   '2022-11-22/2022-11-22_volt_63_pos_top_light_laser_coinc_laser_cond_calibration/',
                '2022-11-22/2022-11-22_volt_65_pos_top_light_laser_coinc_laser_cond_calibration/',
                '2022-11-22/2022-11-22_volt_67_pos_top_light_laser_coinc_laser_cond_calibration/',
                '2022-11-22/2022-11-22_volt_69_pos_top_light_laser_coinc_laser_cond_calibration/',
                '2022-11-22/2022-11-22_volt_71_pos_top_light_laser_coinc_laser_cond_calibration/',
                '2022-11-22/2022-11-22_volt_63_pos_bottom_light_laser_coinc_laser_cond_calibration/',
                '2022-11-22/2022-11-22_volt_65_pos_bottom_light_laser_coinc_laser_cond_calibration/',
                '2022-11-22/2022-11-22_volt_67_pos_bottom_light_laser_coinc_laser_cond_calibration/',
                '2022-11-22/2022-11-22_volt_69_pos_bottom_light_laser_coinc_laser_cond_calibration/',
                '2022-11-22/2022-11-22_volt_71_pos_bottom_light_laser_coinc_laser_cond_calibration/']
    nsubdir = [20]*10
    outfile = [ 'SPE_waveform_1122_63V_top',
                'SPE_waveform_1122_65V_top',
                'SPE_waveform_1122_67V_top',
                'SPE_waveform_1122_69V_top',
                'SPE_waveform_1122_71V_top',
                'SPE_waveform_1122_63V_bottom',
                'SPE_waveform_1122_65V_bottom',
                'SPE_waveform_1122_67V_bottom',
                'SPE_waveform_1122_69V_bottom',
                'SPE_waveform_1122_71V_bottom']
    #####
    # Read A1 sigma1 from calibration_xxV.csv
    #####
    bias = [63,65,67,69,71]
    volt = bias[index%5]
    A1min = []
    A1max = []
    with open('../calibration_1122_{}V.csv'.format(volt)) as f:
        r = csv.reader(f)
        line_count = 0
        for row in r:
            if line_count>0:
                if index//5==0:
                    if line_count<=4:
                        A1min.append(float(row[1]))
                        A1max.append(float(row[2]))
                else:
                    if line_count>4:
                        A1min.append(float(row[1]))
                        A1max.append(float(row[2]))
            line_count += 1
    for ch in range(4):
        print('Ch{} A1min={:.3f} A1max={:.3f}'.format(ch, A1min[ch], A1max[ch]))
    
    file = ROOT.TFile('../root/{}.root'.format(outfile[index]), 'recreate')
    data = ds.Dataset('', mode='spewf', pol=-1, channels=range(4))
    for i in range(nsubdir[index]):
        subdir = '{}{}{}/'.format(data_dir, direc[index],i)
        print(subdir)
        for ch in range(4):    
            data.ch[ch].path = subdir
            data.ch[ch].read_data(header=True)
            data.ch[ch].baseline_subtraction(analysis=True)
            data.ch[ch].ar_filter(tau=20)
            data.ch[ch].get_famp()
            data.ch[ch].get_spe_sumwf(famp_range=(A1min[ch], A1max[ch]))
            data.ch[ch].clear_all()
    for ch in range(4):
        data.ch[ch].get_spe_avgwf()
        print('{} waveforms averaged'.format(data.ch[ch].spe_sumwf_count))
        gSPE = ROOT.TGraph(int(data.ch[ch].samples), array('f', list(data.ch[ch].time)), array('f', list(data.ch[ch].spe_avgwf)))
        gSPE.Write('gSPE_ch{}'.format(ch))
    file.Close()

