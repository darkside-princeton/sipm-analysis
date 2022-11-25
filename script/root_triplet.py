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

    files = []
    #Top: 23 datasets
    for i in range(16):
        files.append(['2022-11-08/2022-11-08_volt_65_pos_top_light_scintillation_coinc_111_cond_with_gamma_overnight/{}/'.format(i),1,'gamma_1108_{}_65V_top'.format(i)])
    files.append(['2022-11-09/2022-11-09_volt_65_pos_top_light_scintillation_coinc_111_cond_with_gamma_high_stat/',20,'gamma_1109_65V_top'])
    files.append(['2022-11-10/2022-11-10_volt_65_pos_top_light_scintillation_coinc_111_cond_with_gamma_high_stat/',20,'gamma_1110_65V_top'])
    files.append(['2022-11-14/2022-11-14_volt_65_pos_top_light_scintillation_coinc_111_cond_gamma_high_stat/',20,'gamma_1114_65V_top'])
    files.append(['2022-11-15/2022-11-15_volt_65_pos_top_light_scintillation_coinc_111_cond_gamma_high_stat/',20,'gamma_1115_65V_top'])
    files.append(['2022-11-17/2022-11-17_volt_65_pos_top_light_scintillation_coinc_000_cond_gamma/',20,'gamma_1117_65V_top'])
    files.append(['2022-11-18/2022-11-18_volt_65_pos_top_light_scintillation_coinc_000_cond_gamma/',20,'gamma_1118_65V_top'])
    files.append(['2022-11-20/2022-11-20_volt_65_pos_top_light_scintillation_coinc_000_cond_gamma/',20,'gamma_1120_65V_top'])
    #Bottom: 70 datasets
    for i in range(6):
        files.append(['2022-11-07/2022-11-07_volt_65_pos_bottom_light_scintillation_coinc_111_cond_with_gamma_purification_{}/'.format(i+1),10,'gamma_1107_p{}_65V_bottom'.format(i+1)])
    for i in range(14):
        files.append(['2022-11-07/2022-11-07_volt_65_pos_bottom_light_scintillation_coinc_111_cond_with_gamma_purification_overnight/{}/'.format(i),1,'gamma_1107_{}_65V_bottom'.format(i)])
    files.append(['2022-11-09/2022-11-09_volt_65_pos_bottom_light_scintillation_coinc_111_cond_with_gamma_high_stat/',20,'gamma_1109_65V_bottom'])
    files.append(['2022-11-10/2022-11-10_volt_65_pos_bottom_light_scintillation_coinc_111_cond_gamma_high_stat/',20,'gamma_1110_65V_bottom'])
    for i in range(20):
        files.append(['2022-11-10/2022-11-10_volt_65_pos_bottom_light_scintillation_coinc_111_cond_with_gamma_overnight/{}/'.format(i),1,'gamma_1110_{}_65V_bottom'.format(i)])
    files.append(['2022-11-15/2022-11-15_volt_65_pos_bottom_light_scintillation_coinc_111_cond_gamma_high_stat/',20,'gamma_1115_65V_bottom'])
    for i in range(6):
        files.append(['2022-11-15/2022-11-15_volt_65_pos_bottom_light_scintillation_coinc_111_cond_gamma_overnight/{}/'.format(i),1,'gamma_1115_{}_65V_bottom'.format(i)])
    files.append(['2022-11-16/2022-11-16_volt_65_pos_bottom_light_scintillation_coinc_111_cond_gamma/',20,'gamma_1116_65V_bottom'])
    for i in range(6):
        files.append(['2022-11-16/2022-11-16_volt_65_pos_bottom_light_scintillation_coinc_111_cond_gamma_overnight/{}/'.format(i),1,'gamma_1116_{}_65V_bottom'.format(i)])
    files.append(['2022-11-17/2022-11-17_volt_65_pos_bottom_light_scintillation_coinc_000_cond_gamma/',20,'gamma_1117_65V_bottom'])
    files.append(['2022-11-18/2022-11-18_volt_65_pos_bottom_light_scintillation_coinc_000_cond_gamma/',20,'gamma_1118_65V_bottom'])
    for i in range(11):
        files.append(['2022-11-18/2022-11-18_volt_65_pos_bottom_light_scintillation_coinc_000_cond_gamma_overnight/{}/'.format(i),1,'gamma_1118_{}_65V_bottom'.format(i)])
    files.append(['2022-11-20/2022-11-20_volt_65_pos_bottom_light_scintillation_coinc_000_cond_gamma/',20,'gamma_1120_65V_bottom'])

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
                    if index<23:
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

    pe_cut = 0
    if index<23:
        pe_cut = 50
    else:
        pe_cut = 300

    data = ds.Dataset('', mode='scintillation', spe=spe_gain, pol=-1, channels=range(4), root_file_name='../root/{}.root'.format(files[index][2]))
    for i in range(files[index][1]):
        if files[index][1]==1:
            subdir = '{}{}'.format(data_dir, files[index][0])
        else:
            subdir = '{}{}{}/'.format(data_dir, files[index][0], i)
        print(subdir)
        for ch in range(4):    
            data.ch[ch].path = subdir
            data.ch[ch].read_data(header=True)
            data.ch[ch].baseline_subtraction(analysis=True)
            data.ch[ch].get_integral(prompt=0.5, long=5.0)
            data.ch[ch].clear()
        data.get_summed_integral_pe()
        data.get_fprompt()
        data.get_avgwf_all(pe_range=(pe_cut,5000), fprompt_range=(0.1,1))
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
