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
    direc = [   '2022-11-01/2022-11-01_volt_61_pos_top_light_laser_coinc_none_cond_calibration',
                '2022-11-01/2022-11-01_volt_63_pos_top_light_laser_coinc_none_cond_calibration',
                '2022-11-01/2022-11-01_volt_65_pos_top_light_laser_coinc_none_cond_calibration',
                '2022-11-01/2022-11-01_volt_67_pos_top_light_laser_coinc_none_cond_calibration',
                '2022-11-01/2022-11-01_volt_69_pos_top_light_laser_coinc_none_cond_calibration',
                '2022-11-01/2022-11-01_volt_61_pos_bottom_light_laser_coinc_none_cond_calibration',
                '2022-11-01/2022-11-01_volt_63_pos_bottom_light_laser_coinc_none_cond_calibration',
                '2022-11-01/2022-11-01_volt_65_pos_bottom_light_laser_coinc_none_cond_calibration',
                '2022-11-01/2022-11-01_volt_67_pos_bottom_light_laser_coinc_none_cond_calibration',
                '2022-11-01/2022-11-01_volt_69_pos_bottom_light_laser_coinc_none_cond_calibration']
    nsubdir = [1,1,1,1,1,1,1,1,1,1]
    outfile = [ 'SPE_waveform_61V_top',
                'SPE_waveform_63V_top',
                'SPE_waveform_65V_top',
                'SPE_waveform_67V_top',
                'SPE_waveform_69V_top',
                'SPE_waveform_61V_bottom',
                'SPE_waveform_63V_bottom',
                'SPE_waveform_65V_bottom',
                'SPE_waveform_67V_bottom',
                'SPE_waveform_69V_bottom']
    #####
    # Read A1 sigma1 from calibration_xxV.csv
    #####
    bias = [61,63,65,67,69]
    volt = bias[index%5]
    A1pe = []
    sigma1pe = []
    with open('../calibration_{}V.csv'.format(volt)) as f:
        r = csv.reader(f)
        line_count = 0
        for row in r:
            if line_count>0:
                if index//5==0:
                    if line_count<=4:
                        A1pe.append(float(row[1]))
                        sigma1pe.append(float(row[2]))
                else:
                    if line_count>4:
                        A1pe.append(float(row[1]))
                        sigma1pe.append(float(row[2]))
            line_count += 1
    for ch in range(4):
        print('Ch{} A1pe={:.3f} sigma1pe={:.3f}'.format(ch, A1pe[ch], sigma1pe[ch]))
    
    file = ROOT.TFile('{}.root'.format(outfile[index]), 'recreate')
    data = ds.Dataset('', mode='spewf', pol=-1, channels=range(4))
    for i in range(nsubdir[index]):
        # subdir = '{}{}/'.format(direc, i)
        subdir = '{}{}/'.format(data_dir, direc[index])
        print(subdir)
        for ch in range(4):    
            data.ch[ch].path = subdir
            data.ch[ch].read_data(header=True)
            data.ch[ch].baseline_subtraction(analysis=True)
            data.ch[ch].ar_filter(tau=20)
            data.ch[ch].get_famp()
            data.ch[ch].get_spe_sumwf(famp_range=(A1pe[ch]-3*sigma1pe[ch], A1pe[ch]+3*sigma1pe[ch]))
            data.ch[ch].clear_all()
    for ch in range(4):
        data.ch[ch].get_spe_avgwf()
        print('{} waveforms averaged'.format(data.ch[ch].spe_avgwf_count))
        gSPE = ROOT.TGraph(int(data.ch[ch].samples), array('f', list(data.ch[ch].time)), array('f', list(data.ch[ch].spe_avgwf)))
        gSPE.Write('gSPE_ch{}'.format(ch))
    file.Close()

