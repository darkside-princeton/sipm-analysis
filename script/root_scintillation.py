import numpy as np
import sipm.sipm as sipm
import sipm.dataset as ds
import sipm.constants as const
import ROOT
import sys
import os
from array import array

if __name__ == '__main__':
    #usage: python root_scintillation.py <data directory> <number of subdirectories> <output file name> <calibration file>
    direc = sys.argv[1]
    nsubdir = int(sys.argv[2])
    outfile = sys.argv[3]
    calib_file = sys.argv[4]
    ##########
    # SPE gain is preliminary. Should use Qavg/(1-p) obtained from calibration analysis!
    ##########
    ##########
    # Read Qavg/(1-p) from calibration file
    ##########
    spe_gain = [550,550,550,550] 
    data = ds.Dataset('', mode='scintillation', spe=spe_gain, pol=-1, channels=range(4), root_file_name='{}.root'.format(outfile))
    for i in range(nsubdir):
        subdir = '{}{}/'.format(direc, i)
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
