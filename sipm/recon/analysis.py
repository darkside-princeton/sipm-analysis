import argparse
import sipm.io.dataset as ds
import sipm.io.h5_io as h5_io

parser = argparse.ArgumentParser("Princeton SiPM Analysis")
parser.add_argument("-n", "--num_events", type=int, default=100000)
parser.add_argument("-f", "--file_dir", type=str, default="")
parser.add_argument("-s", '--sum', type=bool, default=True)
parser.add_argument("-c", '--clear', type=bool, default=True)
args = parser.parse_args()

def main():
    d = ds.Dataset(path=args.file_dir)
    d.analyze(num_events=args.num_events, clear=args.clear, sum=args.sum)
    io = h5_io.IO(dataset=d)
    io.save()
    print('Done...')

if __name__ == "__main__":
    main()