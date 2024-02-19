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
    d = wfd.WaveformDataset(path=args.file_dir, samples=4000, channels=np.arange(8))

    # Run waveform shape analysis on scintillation data
    d.read_calibration_h5(args.calib_file)
    for i in d.channels:
        d.ch[i].read_data(header=True, num_events=args.num_events)
        d.ch[i].baseline_subtraction(samples=d.ch[i].trigger_position-int(0.5/d.ch[i].sample_step))
        d.ch[i].get_integral(length_us=[0.3,5]) # 0.3us for Fprompt analysis
        d.ch[i].get_max()
    d.get_total_pe()
    d.get_fprompt(tprompt=[0.3],channels=np.array([1,2,4,7]))
    # Make cut on total pe, fprompt, baseline rms of all the channels
    cut = (np.array(d.output['total_pe'])<args.pe[1]) & \
        (np.array(d.output['total_pe'])>args.pe[0]) & \
        (np.array(d.output['fprompt_0p30us_1247'])<args.fprompt[1]) & \
        (np.array(d.output['fprompt_0p30us_1247'])>args.fprompt[0])
    for i in d.channels:
        cut = cut & (np.array(d.ch[i].output['baseline_rms'])<d.calib_df['bsl_rms'][i]) & (np.array(d.ch[i].output['amplitude'])<d.calib_df['max_amp'][i])
    # Store average LAr scintillation waveform and number of selected waveforms
    for i in d.channels:
        d.ch[i].output['n_scint_wfs'] = np.sum(cut)
        d.ch[i].output['avg_scint_wf'] = np.dot(d.ch[i].traces.T,cut)/d.ch[i].output['n_scint_wfs']
        d.ch[i].output['time'] = d.ch[i].time
    # Clean up unnecessary variables
    d.output.pop('total_pe')
    d.output.pop('fprompt_0p30us_1247')
    for i in d.channels:
        d.ch[i].output.pop('baseline_mean')
        d.ch[i].output.pop('baseline_rms')
        d.ch[i].output.pop('integral_0p30us')
        d.ch[i].output.pop('integral_5p00us')
        d.ch[i].output.pop('amplitude')
    d.clear()

    # Create a IO objects to save the high level information
    io = h5_io.IO(dataset=d)

    # Save data to HDF5
    io.save(script='scintillation_waveform_liq6')

if __name__ == "__main__":
    main()