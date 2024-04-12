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
    d = wfd.WaveformDataset(path=args.file_dir, samples=4000, channels=np.arange(8))

    # Run pulse analysis on scintillation data

    d.read_calibration_h5(args.calib_file)
    
    for i in d.channels:
        d.ch[i].read_data(header=True, num_events=args.num_events)
        d.ch[i].baseline_subtraction(samples=d.ch[i].trigger_position-int(0.5/d.ch[i].sample_step))
        d.ch[i].get_integral(length_us=[0.3,5]) # <=0.5us for Fprompt analysis
        d.ch[i].output['fired'] = d.ch[i].output['integral_5p00us']>0.5*d.calib_df['cn_corrected_gain'][i]
        d.ch[i].get_max()
    d.output['nch'] = np.zeros_like(d.ch[0].output['fired'])
    d.output['nch_top'] = np.zeros_like(d.ch[0].output['fired'])
    d.output['nch_bot'] = np.zeros_like(d.ch[0].output['fired'])
    for i in d.channels:
        d.output['nch'] = d.output['nch'] + d.ch[i].output['fired'].astype(int)
    for i in [0,1,2,3]:
        d.output['nch_bot'] = d.output['nch_bot'] + d.ch[i].output['fired'].astype(int)
    for i in [4,5,6,7]:
        d.output['nch_top'] = d.output['nch_top'] + d.ch[i].output['fired'].astype(int)
    d.get_total_pe()
    d.get_fprompt(tprompt=[0.3],channels=np.arange(8))
    d.get_fprompt(tprompt=[0.3],channels=np.arange(4))
    d.get_fprompt(tprompt=[0.3],channels=np.arange(4,8))
    d.clear()

    # Create a IO objects to save the high level information
    io = h5_io.IO(dataset=d)

    # Save data to HDF5
    io.save(script='scintillation_pulse_liq5')

if __name__ == "__main__":
    main()