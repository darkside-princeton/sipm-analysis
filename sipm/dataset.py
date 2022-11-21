import numpy as np
import yaml 
import ROOT
from array import array
import sipm.sipm as sipm

class Dataset: 
    def __init__(self,  path, pol=1, channels=range(4), spe=[], root_file_name=None):
        self.path = path
        self.channels = channels
        self.pol = pol
        self.ch = self.InitializeChannels()
        self.summed_integral_pe = []
        self.fprompt = []
        self.gain=spe
        # output root file
        self.root_file_name = root_file_name
        self.root_file = None
        self.root_tree = None
        # variables
        #quality
        self.bsl_avg = array('f', [0]*4)
        self.bsl_med = array('f', [0]*4)
        self.bsl_std = array('f', [0]*4)
        self.acq_max = array('f', [0]*4)
        self.acq_min = array('f', [0]*4)
        #amplitude
        self.fil_amp = array('f', [0]*4)
        #integral
        self.int_long = array('f', [0]*4)
        self.int_prompt = array('f', [0]*4)
        self.initialize_tree()
        
    def InitializeChannels(self):
        channels = []
        for i in self.channels:
            new_channel = sipm.SiPM(id=i, pol=self.pol, path=self.path)
            channels.append(new_channel)
        return channels

    def initialize_tree(self):
        if self.root_file_name==None:
            print('No root file name provided. Use default name tree.root')
            self.root_file_name = 'tree.root'
        self.root_file = ROOT.TFile(self.root_file_name, 'recreate')
        self.root_tree = ROOT.TTree('tree', 'tree')
        # Set branches
        self.root_tree.Branch('bsl_avg', self.bsl_avg, 'bsl_avg[4]/F')
        self.root_tree.Branch('bsl_med', self.bsl_med, 'bsl_med[4]/F')
        self.root_tree.Branch('bsl_std', self.bsl_std, 'bsl_std[4]/F')
        self.root_tree.Branch('acq_max', self.acq_max, 'acq_max[4]/F')
        self.root_tree.Branch('acq_min', self.acq_max, 'acq_min[4]/F')
        self.root_tree.Branch('fil_amp', self.fil_amp, 'fil_amp[4]/F')
        self.root_tree.Branch('int_long', self.int_long, 'int_long[4]/F')
        self.root_tree.Branch('int_prompt', self.int_prompt, 'int_prompt[4]/F')

    def fill_tree(self):
        for ch in self.channels:
            print('ch {} nevents {}'.format(ch, self.ch[ch].nevents))
        for ev in range(self.ch[0].nevents):
            for ch in self.channels:
                self.bsl_avg[ch] = self.ch[ch].baseline_avg[ev]
                self.bsl_med[ch] = self.ch[ch].baseline_med[ev]
                self.bsl_std[ch] = self.ch[ch].baseline_std[ev]
                self.acq_min[ch] = self.ch[ch].acquisition_min[ev]
                self.acq_max[ch] = self.ch[ch].acquisition_max[ev]
                self.fil_amp[ch] = self.ch[ch].famp[ev]
                self.int_long[ch] = self.ch[ch].integral_long[ev]
                self.int_prompt[ch] = self.ch[ch].integral_prompt[ev]
            self.root_tree.Fill()
            if ev>0 and ev%1000==0:
                print('{} events filled'.format(ev))

    def write_tree(self):
        self.root_file.cd()
        self.root_tree.Write("tree")
        self.root_file.Close()

    def get_summed_integral_pe(self):
        self.summed_integral_pe = np.zeros(self.ch[0].cumulative_nevents)
        for i in self.channels:
#             print('ch {} nevents {} charge entries {}'.format(i,self.ch[i].cumulative_nevents, np.shape(self.ch[i].integral_long)[0]))
#             print(np.shape(self.ch[i].integral_long))
#             print(np.shape(self.gain[i]))
            self.summed_integral_pe += np.array(self.ch[i].integral_long)/self.gain[i]

    def get_fprompt(self):
        summed_prompt_integral_pe = np.zeros(self.ch[0].cumulative_nevents)
        for i in self.channels:
            summed_prompt_integral_pe += np.array(self.ch[i].integral_prompt)/self.gain[i]
        self.fprompt = summed_prompt_integral_pe/self.summed_integral_pe

    def get_waveforms_id(self, count=-1, integral_range=(0,1e4), fprompt_range=(0,1)):
        count_ = 0
        event_id = []
        ev = self.ch[0].cumulative_nevents-self.ch[0].nevents
        while (count==-1 or count_<count) and ev < self.ch[0].cumulative_nevents:
            if self.summed_integral_pe[ev]<integral_range[1] and self.summed_integral_pe[ev]>integral_range[0]:
                if self.fprompt[ev]<fprompt_range[1] and self.fprompt[ev]>fprompt_range[0]:
                    event_id.append(ev-(self.ch[0].cumulative_nevents-self.ch[0].nevents))
                    count_ += 1
            ev += 1
        return event_id

    def get_avgwf_all(self, count=-1, integral_range=(0,1e4), fprompt_range=(0,1)):
        indices = self.get_waveforms_id(count, integral_range, fprompt_range)
        for ch in self.channels:
            self.ch[ch].clear()
            self.ch[ch].read_data(simple=True)
            self.ch[ch].baseline_subtraction(analysis=False)
            self.ch[ch].get_avgwf(indices)
            self.ch[ch].clear()