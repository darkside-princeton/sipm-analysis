from sipm.analysis.AdvancedAnalyzer import AdvancedAnalyzer
from typing import Dict, List


class SipmCalibration(AdvancedAnalyzer):
    def __init__(self, positions:List[str], channels:List[int], voltages:List[float], directory:str, metadata_dict:Dict, wf:bool, merge:bool, verbose:bool):
        super().__init__(directory, metadata_dict, wf, merge, verbose)
        self.positions = positions
        self.channels = channels
        self.voltages = voltages