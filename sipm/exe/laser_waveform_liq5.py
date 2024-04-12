import argparse
import sipm.recon.WaveformDataset as wfd
import sipm.recon.h5_io as h5_io
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser("Princeton SiPM Analysis")
parser.add_argument("-n", "--num_events", type=int, default=100000)
parser.add_argument("-f", "--file_dir", type=str, default="")
parser.add_argument("-s", '--calib_file', type=str, default="")
args = parser.parse_args()

def main():
    # Create new dataset object
    d = wfd.WaveformDataset(path=args.file_dir, samples=4000, channels=np.arange(8))

    # Run waveform shape analysis on laser data
    d.calib_df = pd.read_hdf(args.calib_file, key=f'{d.volt}V')
    d.calib_df['cn_corrected_gain'] = d.calib_df['Qpeak']*(1+d.calib_df['Qap'])/(1-d.calib_df['DiCT']) # effective SPE gain corrected for correlated noises (DiCT and afterpulsing)
    for i in d.channels:
        if i in d.calib_df['channel']:
            d.ch[i].read_data(header=True, num_events=args.num_events)
            d.ch[i].baseline_subtraction(samples=d.ch[i].trigger_position-int(0.5/d.ch[i].sample_step))
            d.ch[i].ar_filter(tau=20) # 20 samples = 80us = fast component
            d.ch[i].get_max(ar=True, trig=True) # AR matched filter, maximum near trigger position
            d.ch[i].get_integral() # full integral (from trigger-10 samples to end)
            # Make cut on filtered amplitude->SPE, baseline rms->pre-trigger pulses, and total integral->post-trigger scintillation pulses
            cut = (np.array(d.ch[i].output['baseline_rms'])<d.calib_df['bsl_rms'][i]) & \
                (np.array(d.ch[i].output['amplitude_trig'])<d.calib_df['A1max'][i]) & \
                (np.array(d.ch[i].output['amplitude_trig'])>d.calib_df['A1min'][i]) & \
                (np.array(d.ch[i].output['integral'])<2*d.calib_df['Qpeak'][i])
            # Store SPE average waveform and number of selected waveforms
            d.ch[i].output['n_spe_wfs'] = np.sum(cut)
            d.ch[i].output['avg_spe_wf'] = np.dot(d.ch[i].traces.T,cut)/d.ch[i].output['n_spe_wfs']
            d.ch[i].output['time'] = d.ch[i].time
            # Clean up unnecessary variables
            d.ch[i].output.pop('baseline_mean')
            d.ch[i].output.pop('baseline_rms')
            d.ch[i].output.pop('amplitude_trig')
            d.ch[i].output.pop('peakpos_trig')
            d.ch[i].output.pop('integral')
    d.clear()

    # Create a IO objects to save the high level information
    io = h5_io.IO(dataset=d)

    # Save data to HDF5
    io.save(script='laser_waveform_liq5')

if __name__ == "__main__":
    main()