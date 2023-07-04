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
    d = wfd.WaveformDataset(path=args.file_dir, samples=4000)

    d.read_calibration_h5(args.calib_file)
    for i in d.channels:
        d.ch[i].read_data(header=True, num_events=args.num_events)
        d.ch[i].baseline_subtraction(samples=d.ch[i].trigger_position-int(0.5/d.ch[i].sample_step))

        
        d.ch[i].get_fft2()
        # Store 1PE power spectrum (|fft|^2)
        
        d.ch[i].output['frequency_MHz'] = d.ch[i].time/d.ch[i].sample_step**2/d.ch[i].samples
        # Store 0PE power spectrum (|fft|^2)
        d.ch[i].output['n_0pe_wfs'] = np.sum(cut_0pe)
        d.ch[i].output['avg_0pe_psd'] = np.dot(d.ch[i].fft2.T,cut_0pe)/d.ch[i].output['n_0pe_wfs']
        # Clean up unnecessary variables
        d.ch[i].output.pop('baseline_mean')
        d.ch[i].output.pop('baseline_rms')
        d.ch[i].clear()

    # Create a IO objects to save the high level information
    io = h5_io.IO(dataset=d)

    # Save data to HDF5
    io.save(script='fft_below_breakdown')

if __name__ == "__main__":
    main()