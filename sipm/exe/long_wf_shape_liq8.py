import argparse
import sipm.recon.WaveformDataset as wfd
import sipm.recon.h5_io as h5_io
import numpy as np

parser = argparse.ArgumentParser("Princeton SiPM Analysis")
parser.add_argument("-n", "--num_events", type=int, default=100000)
parser.add_argument("-f", "--file_dir", type=str, default="")
parser.add_argument("-s", '--calib_file', type=str, default="")
parser.add_argument("-p", '--fprompt', type=float, nargs=2, default=[0.1,0.6])
parser.add_argument("-e", '--pe', type=float, nargs=2, default=[300,700])
args = parser.parse_args()

def main():
    # Create new dataset object
    d = wfd.WaveformDataset(path=args.file_dir,channels=[1,2,4,7],trig=9960)

    # Run waveform shape analysis on scintillation data
    d.read_calibration_h5(args.calib_file)
    for i in d.channels:
        d.ch[i].read_data(header=True, num_events=args.num_events)
        d.ch[i].baseline_subtraction(samples=d.ch[i].trigger_position)
        d.ch[i].get_max()
    # Make cut on baseline rms of all the channels
    for i in d.channels:
        cut = (np.array(d.ch[i].output['amplitude'])<d.calib_df['max_amp'][i]) & (np.array(d.ch[i].output['baseline_rms'])<2.0)
    print(f'pre-trigger cut fraction: {1-np.sum(cut)/cut.shape[0]}')
    # Store average LAr scintillation waveform and number of selected waveforms
    for i in d.channels:
        d.ch[i].output['n_scint_wfs'] = np.sum(cut)
        d.ch[i].output['avg_scint_wf'] = np.dot(d.ch[i].traces.T,cut)/d.ch[i].output['n_scint_wfs']
        d.ch[i].output['time'] = d.ch[i].time
    # Clean up unnecessary variables
    for i in d.channels:
        d.ch[i].output.pop('baseline_mean')
        d.ch[i].output.pop('baseline_rms')
        d.ch[i].output.pop('amplitude')
        d.ch[i].output.pop('peakpos')
    d.clear()

    # Create a IO objects to save the high level information
    io = h5_io.IO(dataset=d)

    # Save data to HDF5
    io.save(script='long_wf_shape_liq8')

if __name__ == "__main__":
    main()