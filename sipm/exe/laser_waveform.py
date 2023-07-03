import argparse
import sipm.recon.WaveformDataset as wfd
import sipm.recon.h5_io as h5_io

parser = argparse.ArgumentParser("Princeton SiPM Analysis")
parser.add_argument("-n", "--num_events", type=int, default=100000)
parser.add_argument("-f", "--file_dir", type=str, default="")
parser.add_argument("-s", '--calib_file', type=str, default="")
args = parser.parse_args()

def main():
    # Create new dataset object
    d = wfd.WaveformDataset(path=args.file_dir, samples=4000)

    # Run waveform shape analysis on laser data
    d.process_laser_waveforms(num_events=args.num_events, calib=args.calib_file)

    # Create a IO objects to save the high level information
    io = h5_io.IO(dataset=d)

    # Save data to HDF5
    io.save(script='laser_waveform')

if __name__ == "__main__":
    main()