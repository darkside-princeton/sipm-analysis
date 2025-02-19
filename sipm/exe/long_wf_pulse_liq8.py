import argparse
import sipm.recon.WaveformDataset as wfd
import sipm.recon.h5_io as h5_io
import numpy as np

parser = argparse.ArgumentParser("Princeton SiPM Analysis")
parser.add_argument("-n", "--num_events", type=int, default=100000)
parser.add_argument("-f", "--file_dir", type=str, default="")
parser.add_argument("-s", '--calib_file', type=str, default="")
args = parser.parse_args()

def main():
    # Create new dataset object
    d = wfd.WaveformDataset(path=args.file_dir,channels=[0,1,2,3,4,7],trig=9960)

    # Run pulse analysis on scintillation data

    d.read_calibration_h5(args.calib_file)
    
    for i in d.channels:
        d.ch[i].read_data(header=True, num_events=args.num_events)
        d.ch[i].baseline_subtraction(samples=d.ch[i].trigger_position)
        d.ch[i].get_integral(length_us=[0.3,10]) # <=0.5us for Fprompt analysis
        d.ch[i].get_max()
    d.get_total_pe(length_us=10,channels=[1,2,4,7])
    d.get_fprompt(tprompt=[0.3],channels=[1,2,4,7],t_all=10)
    d.clear()

    # Create an IO objects to save the high level information
    io = h5_io.IO(dataset=d)

    # Save data to HDF5
    io.save(script='long_wf_pulse_liq8')

if __name__ == "__main__":
    main()