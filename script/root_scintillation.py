import numpy as np
import sipm.sipm as sipm
import sipm.dataset as ds
import sipm.constants as const
import ROOT
import sys
import os

if __name__ == '__main__':
    #usage: python root_scintillation.py <data directory> <number of subdirectories> <output file name>
    direc = sys.argv[1]
    nsubdir = int(sys.argv[2])
    outfile = sys.argv[3]
    ##########
    # SPE gain is preliminary. Should use Qavg/(1-p) obtained from calibration analysis!
    ##########
    data = ds.Dataset('', mode='scintillation', spe=[550,550,550,550], pol=-1, channels=range(4), root_file_name='{}.root'.format(outfile))
    for i in range(nsubdir):
        subdir = '{}{}/'.format(direc, i)
        print(subdir)
        for ch in range(4):    
            data.ch[ch].path = subdir
            data.ch[ch].read_data(header=True)
            data.ch[ch].baseline_subtraction(analysis=True)
            data.ch[ch].get_integral(prompt=0.5, long=5.0)
            data.ch[ch].clear()
        data.fill_tree()
        for ch in range(4):
            data.ch[ch].clear_all()
        data.clear()
    data.write_tree()
