import argparse

import sipm.io.sipm as sipm
import sipm.io.dataset as ds
import sipm.util.functions as func

parser = argparse.ArgumentParser("Princeton SiPM Analysis")
parser.add_argument("-n", "--num_events", type=int, default=100000)
parser.add_argument("-f", "--file_dir", type=str, default="")
parser.add_argument("-s", '--sum', type=bool, default=True)
parser.add_argument("-c", '--clear', type=bool, default=True)
args = parser.parse_args()

def main():
    d = ds.Dataset(path=args.file_dir)
    d.analyze(num_events=args.num_events, clear=args.clear, sum=args.sum)
    print('Done...')
    # d.save()

if __name__ == "__main__":
    main()