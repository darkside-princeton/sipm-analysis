import argparse
import sipm.recon.WaveformDataset as wfd
import sipm.recon.h5_io as h5_io
import numpy as np

parser = argparse.ArgumentParser("Princeton SiPM Analysis")
parser.add_argument("-n", "--num_events", type=int, default=100000)
parser.add_argument("-f", "--file_dir", type=str, default="")
parser.add_argument("-s", '--sipm_calib_dir', type=str, default="")
args = parser.parse_args()

def main():
    # Create new dataset object
    d = wfd.WaveformDataset(path=args.file_dir, samples=4000)

    d.read_calibration(args.sipm_calib_dir)
    for i in d.channels:
        d.ch[i].read_data(header=True, num_events=args.num_events)
        d.ch[i].baseline_subtraction(samples=d.ch[i].trigger_position-int(0.5/d.ch[i].sample_step))
        d.ch[i].ar_filter(tau=20) # 20 samples = 80us = fast component
        d.ch[i].get_max(ar=True, trig=True) # AR matched filter, maximum near trigger position
        d.ch[i].get_integral() # full integral (from trigger-10 samples to end)
        # Make cut on filtered amplitude->SPE, baseline rms->pre-trigger pulses, and total integral->post-trigger scintillation pulses
        cut_1pe = (np.array(d.ch[i].output['baseline_rms'])<d.bslrms[i]) & \
            (np.array(d.ch[i].output['amplitude_trig'])<d.a1max[i]) & \
            (np.array(d.ch[i].output['amplitude_trig'])>d.a1min[i]) & \
            (np.array(d.ch[i].output['integral'])<6*d.gain[i])
        cut_0pe = (np.array(d.ch[i].output['baseline_rms'])<d.bslrms[i]) & \
            (np.array(d.ch[i].output['amplitude_trig'])<d.a1min[i]) & \
            (np.array(d.ch[i].output['integral'])<6*d.gain[i])
        d.ch[i].get_fft2()
        # Store 1PE power spectrum (|fft|^2)
        d.ch[i].output['n_1pe_wfs'] = np.sum(cut_1pe)
        d.ch[i].output['avg_1pe_psd'] = np.dot(d.ch[i].fft2.T,cut_1pe)/d.ch[i].output['n_1pe_wfs']
        d.ch[i].output['frequency_MHz'] = d.ch[i].time/d.ch[i].sample_step**2/d.ch[i].samples
        # Store 0PE power spectrum (|fft|^2)
        d.ch[i].output['n_0pe_wfs'] = np.sum(cut_0pe)
        d.ch[i].output['avg_0pe_psd'] = np.dot(d.ch[i].fft2.T,cut_0pe)/d.ch[i].output['n_0pe_wfs']
        # Clean up unnecessary variables
        d.ch[i].output.pop('baseline_mean')
        d.ch[i].output.pop('baseline_rms')
        d.ch[i].output.pop('amplitude_trig')
        d.ch[i].output.pop('peakpos_trig')
        d.ch[i].output.pop('integral')
        d.ch[i].clear()

    # Create a IO objects to save the high level information
    io = h5_io.IO(dataset=d)

    # Save data to HDF5
    io.save(script='fft')

if __name__ == "__main__":
    main()