import glob
import sipm.util.submit as submit
import yaml
import argparse

parser = argparse.ArgumentParser("Princeton SiPM Analysis")
parser.add_argument("-c", "--config", type=str)
args = parser.parse_args()

with open(args.config, 'r') as f:
    conf = yaml.safe_load(f)
print(conf['description'])

dirs = []
for _ in conf['files']:
    dirs += glob.glob(_)

scheduler = submit.Scheduler(dirs=dirs)

if conf['script'] == 'sipm/exe/analysis.py':
    scheduler.submit(script=conf['script'],args=f"-n {conf['num_events']}")
elif conf['script'] == 'sipm/exe/laser_pulse.py':
    scheduler.submit(script=conf['script'],args=f"-n {conf['num_events']}")
elif conf['script'] == 'sipm/exe/scintillation_pulse.py':
    scheduler.submit(script=conf['script'],args=f"-n {conf['num_events']} -s {conf['sipm_calib_dir']}")
elif conf['script'] == 'sipm/exe/laser_waveform.py':
    scheduler.submit(script=conf['script'],args=f"-n {conf['num_events']} -s {conf['sipm_calib_dir']}")
elif conf['script'] == 'sipm/exe/scintillation_waveform.py':
    if conf.has_key('cpu_memory'):
        scheduler.cpu_memory=conf['cpu_memory']
    scheduler.submit(script=conf['script'],args=f"-n {conf['num_events']} -s {conf['sipm_calib_dir']} -p {conf['fprompt_range'][0]} {conf['fprompt_range'][1]} -e {conf['pe_range'][0]} {conf['pe_range'][1]}")