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
    direc = [   '2022-11-08/2022-11-08_volt_65_pos_top_light_scintillation_coinc_111_cond_with_gamma_overnight/0/',
                '2022-11-08/2022-11-08_volt_65_pos_top_light_scintillation_coinc_111_cond_with_gamma_overnight/1/',
                '2022-11-08/2022-11-08_volt_65_pos_top_light_scintillation_coinc_111_cond_with_gamma_overnight/2/',
                '2022-11-08/2022-11-08_volt_65_pos_top_light_scintillation_coinc_111_cond_with_gamma_overnight/3/',
                '2022-11-08/2022-11-08_volt_65_pos_top_light_scintillation_coinc_111_cond_with_gamma_overnight/4/',
                '2022-11-08/2022-11-08_volt_65_pos_top_light_scintillation_coinc_111_cond_with_gamma_overnight/5/',
                '2022-11-08/2022-11-08_volt_65_pos_top_light_scintillation_coinc_111_cond_with_gamma_overnight/6/',
                '2022-11-08/2022-11-08_volt_65_pos_top_light_scintillation_coinc_111_cond_with_gamma_overnight/7/',
                '2022-11-08/2022-11-08_volt_65_pos_top_light_scintillation_coinc_111_cond_with_gamma_overnight/8/',
                '2022-11-08/2022-11-08_volt_65_pos_top_light_scintillation_coinc_111_cond_with_gamma_overnight/9/',
                '2022-11-08/2022-11-08_volt_65_pos_top_light_scintillation_coinc_111_cond_with_gamma_overnight/10/',
                '2022-11-08/2022-11-08_volt_65_pos_top_light_scintillation_coinc_111_cond_with_gamma_overnight/11/',
                '2022-11-08/2022-11-08_volt_65_pos_top_light_scintillation_coinc_111_cond_with_gamma_overnight/12/',
                '2022-11-08/2022-11-08_volt_65_pos_top_light_scintillation_coinc_111_cond_with_gamma_overnight/13/',
                '2022-11-08/2022-11-08_volt_65_pos_top_light_scintillation_coinc_111_cond_with_gamma_overnight/14/',
                '2022-11-08/2022-11-08_volt_65_pos_top_light_scintillation_coinc_111_cond_with_gamma_overnight/15/',

                '2022-11-09/2022-11-09_volt_65_pos_top_light_scintillation_coinc_111_cond_with_gamma_high_stat/',
                '2022-11-10/2022-11-10_volt_65_pos_top_light_scintillation_coinc_111_cond_with_gamma_high_stat/',
                '2022-11-14/2022-11-14_volt_65_pos_top_light_scintillation_coinc_111_cond_gamma_high_stat/',
                '2022-11-15/2022-11-15_volt_65_pos_top_light_scintillation_coinc_111_cond_gamma_high_stat/',
                '2022-11-17/2022-11-17_volt_65_pos_top_light_scintillation_coinc_000_cond_gamma/',
                '2022-11-18/2022-11-18_volt_65_pos_top_light_scintillation_coinc_000_cond_gamma/',
                '2022-11-20/2022-11-20_volt_65_pos_top_light_scintillation_coinc_000_cond_gamma/',

                '2022-11-07/2022-11-07_volt_65_pos_bottom_light_scintillation_coinc_111_cond_with_gamma_purification_1/',
                '2022-11-07/2022-11-07_volt_65_pos_bottom_light_scintillation_coinc_111_cond_with_gamma_purification_2/',
                '2022-11-07/2022-11-07_volt_65_pos_bottom_light_scintillation_coinc_111_cond_with_gamma_purification_3/',
                '2022-11-07/2022-11-07_volt_65_pos_bottom_light_scintillation_coinc_111_cond_with_gamma_purification_4/',
                '2022-11-07/2022-11-07_volt_65_pos_bottom_light_scintillation_coinc_111_cond_with_gamma_purification_5/',
                '2022-11-07/2022-11-07_volt_65_pos_bottom_light_scintillation_coinc_111_cond_with_gamma_purification_6/',

                '2022-11-07/2022-11-07_volt_65_pos_bottom_light_scintillation_coinc_111_cond_with_gamma_purification_overnight/0/',
                '2022-11-07/2022-11-07_volt_65_pos_bottom_light_scintillation_coinc_111_cond_with_gamma_purification_overnight/1/',
                '2022-11-07/2022-11-07_volt_65_pos_bottom_light_scintillation_coinc_111_cond_with_gamma_purification_overnight/2/',
                '2022-11-07/2022-11-07_volt_65_pos_bottom_light_scintillation_coinc_111_cond_with_gamma_purification_overnight/3/',
                '2022-11-07/2022-11-07_volt_65_pos_bottom_light_scintillation_coinc_111_cond_with_gamma_purification_overnight/4/',
                '2022-11-07/2022-11-07_volt_65_pos_bottom_light_scintillation_coinc_111_cond_with_gamma_purification_overnight/5/',
                '2022-11-07/2022-11-07_volt_65_pos_bottom_light_scintillation_coinc_111_cond_with_gamma_purification_overnight/6/',
                '2022-11-07/2022-11-07_volt_65_pos_bottom_light_scintillation_coinc_111_cond_with_gamma_purification_overnight/7/',
                '2022-11-07/2022-11-07_volt_65_pos_bottom_light_scintillation_coinc_111_cond_with_gamma_purification_overnight/8/',
                '2022-11-07/2022-11-07_volt_65_pos_bottom_light_scintillation_coinc_111_cond_with_gamma_purification_overnight/9/',
                '2022-11-07/2022-11-07_volt_65_pos_bottom_light_scintillation_coinc_111_cond_with_gamma_purification_overnight/10/',
                '2022-11-07/2022-11-07_volt_65_pos_bottom_light_scintillation_coinc_111_cond_with_gamma_purification_overnight/11/',
                '2022-11-07/2022-11-07_volt_65_pos_bottom_light_scintillation_coinc_111_cond_with_gamma_purification_overnight/12/',
                '2022-11-07/2022-11-07_volt_65_pos_bottom_light_scintillation_coinc_111_cond_with_gamma_purification_overnight/13/',

                '2022-11-09/2022-11-09_volt_65_pos_bottom_light_scintillation_coinc_111_cond_with_gamma_high_stat/',
                '2022-11-10/2022-11-10_volt_65_pos_bottom_light_scintillation_coinc_111_cond_gamma_high_stat/',

                '2022-11-10/2022-11-10_volt_65_pos_bottom_light_scintillation_coinc_111_cond_with_gamma_overnight/0/',
                '2022-11-10/2022-11-10_volt_65_pos_bottom_light_scintillation_coinc_111_cond_with_gamma_overnight/1/',
                '2022-11-10/2022-11-10_volt_65_pos_bottom_light_scintillation_coinc_111_cond_with_gamma_overnight/2/',
                '2022-11-10/2022-11-10_volt_65_pos_bottom_light_scintillation_coinc_111_cond_with_gamma_overnight/3/',
                '2022-11-10/2022-11-10_volt_65_pos_bottom_light_scintillation_coinc_111_cond_with_gamma_overnight/4/',
                '2022-11-10/2022-11-10_volt_65_pos_bottom_light_scintillation_coinc_111_cond_with_gamma_overnight/5/',
                '2022-11-10/2022-11-10_volt_65_pos_bottom_light_scintillation_coinc_111_cond_with_gamma_overnight/6/',
                '2022-11-10/2022-11-10_volt_65_pos_bottom_light_scintillation_coinc_111_cond_with_gamma_overnight/7/',
                '2022-11-10/2022-11-10_volt_65_pos_bottom_light_scintillation_coinc_111_cond_with_gamma_overnight/8/',
                '2022-11-10/2022-11-10_volt_65_pos_bottom_light_scintillation_coinc_111_cond_with_gamma_overnight/9/',
                '2022-11-10/2022-11-10_volt_65_pos_bottom_light_scintillation_coinc_111_cond_with_gamma_overnight/10/',
                '2022-11-10/2022-11-10_volt_65_pos_bottom_light_scintillation_coinc_111_cond_with_gamma_overnight/11/',
                '2022-11-10/2022-11-10_volt_65_pos_bottom_light_scintillation_coinc_111_cond_with_gamma_overnight/12/',
                '2022-11-10/2022-11-10_volt_65_pos_bottom_light_scintillation_coinc_111_cond_with_gamma_overnight/13/',
                '2022-11-10/2022-11-10_volt_65_pos_bottom_light_scintillation_coinc_111_cond_with_gamma_overnight/14/',
                '2022-11-10/2022-11-10_volt_65_pos_bottom_light_scintillation_coinc_111_cond_with_gamma_overnight/15/',
                '2022-11-10/2022-11-10_volt_65_pos_bottom_light_scintillation_coinc_111_cond_with_gamma_overnight/16/',
                '2022-11-10/2022-11-10_volt_65_pos_bottom_light_scintillation_coinc_111_cond_with_gamma_overnight/17/',
                '2022-11-10/2022-11-10_volt_65_pos_bottom_light_scintillation_coinc_111_cond_with_gamma_overnight/18/',
                '2022-11-10/2022-11-10_volt_65_pos_bottom_light_scintillation_coinc_111_cond_with_gamma_overnight/19/',

                '2022-11-15/2022-11-15_volt_65_pos_bottom_light_scintillation_coinc_111_cond_gamma_high_stat/',
                '2022-11-16/2022-11-16_volt_65_pos_bottom_light_scintillation_coinc_111_cond_gamma/',
                '2022-11-17/2022-11-17_volt_65_pos_bottom_light_scintillation_coinc_000_cond_gamma/',
                '2022-11-18/2022-11-18_volt_65_pos_bottom_light_scintillation_coinc_000_cond_gamma/',
                '2022-11-20/2022-11-20_volt_65_pos_bottom_light_scintillation_coinc_000_cond_gamma/'
            ]
    nsubdir = [ 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, 20,20,20,20,20,20,20,  10,10,10,10,10,10, 1,1,1,1,1,1,1,1,1,1,1,1,1,1, 20,20, 20,20,20,20,20]
    outfile = [ 'gamma_1108_0_65V_top',
                'gamma_1108_1_65V_top',
                'gamma_1108_2_65V_top',
                'gamma_1108_3_65V_top',
                'gamma_1108_4_65V_top',
                'gamma_1108_5_65V_top',
                'gamma_1108_6_65V_top',
                'gamma_1108_7_65V_top',
                'gamma_1108_8_65V_top',
                'gamma_1108_9_65V_top',
                'gamma_1108_10_65V_top',
                'gamma_1108_11_65V_top',
                'gamma_1108_12_65V_top',
                'gamma_1108_13_65V_top',
                'gamma_1108_14_65V_top',
                'gamma_1108_15_65V_top',

                'gamma_1109_65V_top',
                'gamma_1110_65V_top',
                'gamma_1114_65V_top',
                'gamma_1115_65V_top',
                'gamma_1117_65V_top',
                'gamma_1118_65V_top',
                'gamma_1120_65V_top',

                'gamma_1107_65V_bottom',
                'gamma_1109_65V_bottom',
                'gamma_1110_65V_bottom',
                'gamma_1115_65V_bottom',
                'gamma_1116_65V_bottom',
                'gamma_1117_65V_bottom',
                'gamma_1118_65V_bottom',
                'gamma_1120_65V_bottom'
            ]
    ##########
    # Read Qavg/(1-p) from calibration file
    ##########
    spe_gain = []
    try:
        with open('../calibration_1122_65V.csv') as f:
            r = csv.reader(f)
            line_count = 0
            for row in r:
                if line_count>0:
                    if index<8:
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
        if nsubdir[index]==1:
            subdir = '{}{}'.format(data_dir, direc[index])
        else:
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
        data.get_avgwf_all(pe_range=(50,5000), fprompt_range=(0.1,1))
        data.fill_tree()
        for ch in range(4):
            data.ch[ch].clear_all()
        data.clear()
    for ch in range(4):
        data.ch[ch].get_scint_avgwf()
        print('Ch{} {} waveforms averaged'.format(ch, data.ch[ch].scint_sumwf_count))
        gScintWf = ROOT.TGraph(int(data.ch[ch].samples), array('f', list(data.ch[ch].time)), array('f', list(data.ch[ch].scint_avgwf)))
        gScintWf.Write('gScintWf_ch{}'.format(ch))
    data.write_tree()
