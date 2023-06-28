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
        cut = (np.array(d.ch[i].output['baseline_rms'])<d.bslrms[i]) & \
            (np.array(d.ch[i].output['amplitude_trig'])<d.a1max[i]) & \
            (np.array(d.ch[i].output['amplitude_trig'])>d.a1min[i])
        # Store SPE average waveform and number of selected waveforms
        d.ch[i].output['n_spe_wfs'] = np.sum(cut)
        if np.sum(cut) > 0:
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
    io.save(wf=True)

if __name__ == "__main__":
    main()