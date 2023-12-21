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
    d.process_scintillation_waveforms(num_events=args.num_events, calib=args.calib_file, fprompt=args.fprompt, pe=args.pe)

    # Create a IO objects to save the high level information
    io = h5_io.IO(dataset=d)

    # Save data to HDF5
    io.save(script='scintillation_waveform_liq5')

if __name__ == "__main__":
    main()