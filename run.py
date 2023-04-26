import glob
import sipm.util.scheduler as scheduler

dirs = glob.glob('/scratch/gpfs/GALBIATI/data/sipm/reflector_studies/2022-11-22/**/*/')

scheduler = scheduler.Scheduler(dirs=dirs)
scheduler.submit()