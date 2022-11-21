import numpy as np
import sipm.sipm as sipm
import sipm.dataset as ds
import sipm.constants as const
import ROOT
import sys
import os

if __name__ == '__main__':
    #usage: python root_calibration.py <data directory> <number of subdirectories> <output file name>
    direc = sys.argv[1]
    nsubdir = int(sys.argv[2])
    outfile = sys.argv[3]

    data = ds.Dataset('', mode='calibration', pol=-1, channels=range(4), root_file_name='{}.root'.format(outfile))
    for i in range(nsubdir):
        # subdir = '{}{}/'.format(direc, i)
        subdir = direc
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
