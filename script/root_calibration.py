import numpy as np
import sipm.sipm as sipm
import sipm.dataset as ds
import sipm.constants as const
import ROOT
import sys
import os

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
    outfile = [ 'calibration_61V_top',
                'calibration_63V_top',
                'calibration_65V_top',
                'calibration_67V_top',
                'calibration_69V_top',
                'calibration_61V_bottom',
                'calibration_63V_bottom',
                'calibration_65V_bottom',
                'calibration_67V_bottom',
                'calibration_69V_bottom']

    data = ds.Dataset('', mode='calibration', pol=-1, channels=range(4), root_file_name='{}.root'.format(outfile[index]))
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
            data.ch[ch].get_integral(long=5.0)
            data.ch[ch].clear()
        data.fill_tree()
        for ch in range(4):
            data.ch[ch].clear_all()
    data.write_tree()
