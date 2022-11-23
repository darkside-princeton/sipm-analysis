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
    direc = [   '2022-11-20/2022-11-20_volt_63_pos_top_light_scintillation_coinc_000_cond_gamma/',
                '2022-11-20/2022-11-20_volt_65_pos_top_light_scintillation_coinc_000_cond_gamma/',
                '2022-11-20/2022-11-20_volt_67_pos_top_light_scintillation_coinc_000_cond_gamma/',
                '2022-11-20/2022-11-20_volt_69_pos_top_light_scintillation_coinc_000_cond_gamma/',
                '2022-11-20/2022-11-20_volt_71_pos_top_light_scintillation_coinc_000_cond_gamma/',
                '2022-11-20/2022-11-20_volt_63_pos_bottom_light_scintillation_coinc_000_cond_gamma/',
                '2022-11-20/2022-11-20_volt_65_pos_bottom_light_scintillation_coinc_000_cond_gamma/',
                '2022-11-20/2022-11-20_volt_67_pos_bottom_light_scintillation_coinc_000_cond_gamma/',
                '2022-11-20/2022-11-20_volt_69_pos_bottom_light_scintillation_coinc_000_cond_gamma/',
                '2022-11-20/2022-11-20_volt_71_pos_bottom_light_scintillation_coinc_000_cond_gamma/',
                '2022-11-18/2022-11-18_volt_63_pos_top_light_scintillation_coinc_000_cond_no_gamma/',
                '2022-11-18/2022-11-18_volt_65_pos_top_light_scintillation_coinc_000_cond_no_gamma/',
                '2022-11-18/2022-11-18_volt_67_pos_top_light_scintillation_coinc_000_cond_no_gamma/',
                '2022-11-18/2022-11-18_volt_69_pos_top_light_scintillation_coinc_000_cond_no_gamma/',
                '2022-11-18/2022-11-18_volt_71_pos_top_light_scintillation_coinc_000_cond_no_gamma/',
                '2022-11-18/2022-11-18_volt_63_pos_bottom_light_scintillation_coinc_000_cond_no_gamma/',
                '2022-11-18/2022-11-18_volt_65_pos_bottom_light_scintillation_coinc_000_cond_no_gamma/',
                '2022-11-18/2022-11-18_volt_67_pos_bottom_light_scintillation_coinc_000_cond_no_gamma/',
                '2022-11-18/2022-11-18_volt_69_pos_bottom_light_scintillation_coinc_000_cond_no_gamma/',
                '2022-11-18/2022-11-18_volt_71_pos_bottom_light_scintillation_coinc_000_cond_no_gamma/'
            ]
    nsubdir = [20]*20
    outfile = [ 'gamma_63V_top',
                'gamma_65V_top',
                'gamma_67V_top',
                'gamma_69V_top',
                'gamma_71V_top',
                'gamma_63V_bottom',
                'gamma_65V_bottom',
                'gamma_67V_bottom',
                'gamma_69V_bottom',
                'gamma_71V_bottom',
                'backgrounds_63V_top',
                'backgrounds_65V_top',
                'backgrounds_67V_top',
                'backgrounds_69V_top',
                'backgrounds_71V_top',
                'backgrounds_63V_bottom',
                'backgrounds_65V_bottom',
                'backgrounds_67V_bottom',
                'backgrounds_69V_bottom',
                'backgrounds_71V_bottom'
            ]
    ##########
    # SPE gain is preliminary. Should use Qavg/(1-p) obtained from calibration analysis!
    ##########
    ##########
    # Read Qavg/(1-p) from calibration file
    ##########
    bias = [63,65,67,69,71]
    volt = bias[index%5]
    spe_gain = []
    try:
        with open('../calibration_1122_{}V.csv'.format(volt)) as f:
            r = csv.reader(f)
            line_count = 0
            for row in r:
                if line_count>0:
                    if index%10<5:
                        if line_count<=4:
                            spe_gain.append(float(row[5])/(1-float(row[3])))
                    else:
                        if line_count>4:
                            spe_gain.append(float(row[5])/(1-float(row[3])))
                line_count += 1
    except:
        spe_gain = [550]*4
        print('No calibration csv file. Use default spe_gain of 550.')

    for ch in range(4):
        print('Ch{} spe_gain=[{:.3f},{:.3f},{:.3f},{:.3f}]'.format(ch, *spe_gain))

    data = ds.Dataset('', mode='scintillation', spe=spe_gain, pol=-1, channels=range(4), root_file_name='../root/{}.root'.format(outfile[index]))
    for i in range(nsubdir[index]):
        subdir = '{}{}{}/'.format(data_dir, direc[index], i)
        print(subdir)
        for ch in range(4):    
            data.ch[ch].path = subdir
            data.ch[ch].read_data(header=True)
            data.ch[ch].baseline_subtraction(analysis=True)
            data.ch[ch].get_integral(prompt=0.5, long=5.0)
            data.ch[ch].clear()
        data.get_summed_integral_pe()
        data.get_fprompt()
        data.get_avgwf_all(pe_range=(50,1e4), fprompt_range=(0.2,0.4))
        data.fill_tree()
        for ch in range(4):
            data.ch[ch].clear_all()
        data.clear()
    for ch in range(4):
        data.ch[ch].get_scint_avgwf()
        gScintWf = ROOT.TGraph(int(data.ch[ch].samples), array('f', list(data.ch[ch].time)), array('f', list(data.ch[ch].scint_avgwf)))
        gScintWf.Write('gScintWf_ch{}'.format(ch))
    data.write_tree()
