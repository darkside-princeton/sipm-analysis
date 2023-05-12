from sipm.analysis.AdvancedAnalyzer import AdvancedAnalyzer

class SipmTemplate(AdvancedAnalyzer):
    def __init__(self):
        super().__init__()
        self.results = {'bsl_rms_hist':{},
                        'bsl_mean_hist':{},
                        'fil_amp_hist':{},
                        'dict_fit':{},
                        'charge_hist':{},
                        'charge_fit':{},
                        'ap_prob_fit':{},
                        'vbd_fit':{},
                        'dict':{},
                        'ap_prof':{},
                        'ap_charge':{},
                        'vbd_ch':{}}
    
    