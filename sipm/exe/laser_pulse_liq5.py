import argparse
import sipm.recon.WaveformDataset as wfd
import sipm.recon.h5_io as h5_io
import numpy as np

parser = argparse.ArgumentParser("Princeton SiPM Analysis")
parser.add_argument("-n", "--num_events", type=int, default=100000)
parser.add_argument("-f", "--file_dir", type=str, default="")
args = parser.parse_args()

def main():
    # Create new dataset object
    d = wfd.WaveformDataset(path=args.file_dir, channels=np.arange(8), samples=4000)

    # Run pulse analysis on laser data
    d.process_laser_pulses(num_events=args.num_events)

    # Create a IO objects to save the high level information
    io = h5_io.IO(dataset=d)

    # Save data to HDF5
    io.save(script='laser_pulse_liq5')

if __name__ == "__main__":
    main()