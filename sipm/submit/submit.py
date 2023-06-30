import glob
import sipm.util.scheduler as scheduler
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

scheduler = scheduler.Scheduler(dirs=dirs)

if 'cpu_memory' in conf:
    scheduler.cpu_memory=conf['cpu_memory']

submit_string = ''
for k,v in conf['arguments'].items():
    submit_string += f' --{k}'
    if isinstance(v,list):
        for v_ in v:
            submit_string += f' {v_}'
    else:
        submit_string += f' {v}'
submit_string = submit_string[1:]

scheduler.submit(script=conf['script'],args=submit_string)