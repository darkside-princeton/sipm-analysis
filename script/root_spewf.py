import numpy as np
import sipm.sipm as sipm
import sipm.dataset as ds
import sipm.constants as const
import ROOT
import sys
import os
from array import array

if __name__ == '__main__':
    #usage: python root_spewf.py <data directory> <number of subdirectories> <output file name> <A1> <sigma1>
    direc = sys.argv[1]
    nsubdir = int(sys.argv[2])
    outfile = sys.argv[3]
    A1pe = float(sys.argv[4])
    sigma1pe = float(sys.argv[5])

    file = ROOT.TFile('{}.root'.format(outfile), 'recreate')
    data = ds.Dataset('', mode='spewf', pol=-1, channels=range(4), root_file_name='{}.root'.format(outfile))
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
            data.ch[ch].get_spe_sumwf(famp_range=(A1pe-3*sigma1pe, A1pe+3*sigma1pe))
            data.ch[ch].clear_all()
    for ch in range(4):
        data.ch[ch].get_spe_avgwf()
        print('{} waveforms averaged'.format(data.ch[ch].spe_avgwf_count))
        gSPE = ROOT.TGraph(int(data.ch[ch].samples), array('f', list(data.ch[ch].time)), array('f', list(data.ch[ch].spe_avgwf)))
        gSPE.Write('gSPE_ch{}'.format(ch))
    file.Close()

